import { useEffect, useState } from "react";
import { Condiciones, Informe, Oligo, Proyecto, fmt, post } from "./lib";

// ---------- IDT: comparar la tabla de control con el OligoAnalyzer oficial ----------
const CAMPOS_IDT: [string, keyof Informe, number][] = [
  ["Tm fábrica", "tm_fabrica", 1], ["Tm PCR", "tm_pcr", 1], ["Horquilla ΔG", "horquilla_dg", 2], ["Horquilla Tm", "horquilla_tm", 1],
  ["Autodímero ΔG", "autodimero_dg_fabrica", 2], ["Peso molecular", "pm", 1], ["nmol/OD", "nmol_od", 2], ["µg/OD", "ug_od", 2],
];
type FilaIdt = Record<string, number | null> & { campos_idt?: string[] };
const GUARDA_CLAVE = "primerflow-clave-idt";

export function ComprobarIdt({ res, condiciones }: { res: Informe[]; condiciones: Condiciones }) {
  const [clave, setClave] = useState(() => { try { return localStorage.getItem(GUARDA_CLAVE) ?? ""; } catch { return ""; } });
  const [idt, setIdt] = useState<FilaIdt[] | null>(null);
  const [cargando, setCargando] = useState(false);
  const [error, setError] = useState("");

  async function comprobar() {
    setCargando(true); setError("");
    try { localStorage.setItem(GUARDA_CLAVE, clave); } catch { /* sin almacenamiento: se pide cada vez */ }
    try {
      const r = await post<{ oligos: FilaIdt[] }>("idt", { oligos: res.map((o) => ({ nombre: o.nombre, seq: o.seq })), condiciones, clave });
      setIdt(r.oligos);
    } catch (e) { setError((e as Error).message); }
    setCargando(false);
  }
  const faltan = idt ? CAMPOS_IDT.filter(([, k]) => idt.every((f) => f[k] == null)).map(([t]) => t) : [];

  return (
    <>
      <div className="fila">
        <h3>Comprobar con IDT (OligoAnalyzer oficial)</h3>
        <label>Clave de acceso <input type="password" value={clave} onChange={(e) => setClave(e.target.value)} /></label>
        <button disabled={!clave || cargando} onClick={comprobar}>{cargando ? "Consultando IDT…" : `Comprobar ${res.length} primers con IDT`}</button>
      </div>
      {error && <p className="error">{error}</p>}
      {idt && (
        <div className="desliza">
          <table>
            <thead>
              <tr><th rowSpan={2}>oligo</th>{CAMPOS_IDT.map(([t]) => <th key={t} colSpan={2}>{t}</th>)}</tr>
              <tr>{CAMPOS_IDT.flatMap(([t]) => [<th key={t + "a"}>app</th>, <th key={t + "i"}>IDT</th>])}</tr>
            </thead>
            <tbody>
              {res.map((o, i) => (
                <tr key={o.nombre + o.seq}>
                  <td>{o.nombre}</td>
                  {CAMPOS_IDT.flatMap(([t, k, d]) => {
                    const a = o[k] as number | null, b = idt[i]?.[k] ?? null;
                    const dif = a != null && b != null ? Math.abs(a - b) : null;
                    return [<td key={t + "a"}>{fmt(a, d)}</td>,
                      <td key={t + "i"} className={dif == null ? "" : dif <= 0.5 ? "verde" : dif <= 1.5 ? "amarillo" : "rojo"}>{fmt(b, d)}</td>];
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
      <p className="nota">
        Manda tus primers al OligoAnalyzer de IDT con tu cuenta (sólo quien tenga la clave). Verde: IDT y la app difieren ≤ 0.5; amarillo ≤ 1.5; rojo más.
        {faltan.length > 0 && ` IDT no devolvió: ${faltan.join(", ")} (campos recibidos: ${idt?.[0]?.campos_idt?.join(", ")}).`}
      </p>
    </>
  );
}

// ---------- NCBI: ¿los pares amplifican algo fuera del panel? (BLAST) ----------
type Sitio = { primer: string; accesion: string; especie: string; desajustes: number };
type Producto = { par: string; primers: string; tamano: number; accesion: string; especie: string; titulo: string; desajustes: number };
export type TrabajoBlast = { rid: string; enviado: number; excluir: string; estado: "esperando" | "listo" | "error"; sitios?: Sitio[]; productos?: Producto[]; error?: string };

/** Pares por nombre: Saureus_nuc_cF + Saureus_nuc_cR → «Saureus_nuc». */
function pares(oligos: Oligo[]) {
  const m = new Map<string, { nombre: string; F?: string; R?: string }>();
  for (const o of oligos) {
    const x = o.nombre.match(/^(.*?)_c?([FR])$/);
    if (x) m.set(x[1], { ...(m.get(x[1]) ?? { nombre: x[1] }), [x[2]]: o.nombre });
  }
  return [...m.values()].filter((p) => p.F && p.R);
}

export function Especificidad({ p, setP, especies }: { p: Proyecto; setP: (f: (q: Proyecto) => Proyecto) => void; especies: string[] }) {
  const t = p.blast;
  const oligos = p.oligos.filter((o) => /^[ACGT]{8,}$/i.test(o.seq));
  const [excluir, setExcluir] = useState(especies.join(", "));
  const [enviando, setEnviando] = useState(false);
  const [vuelta, setVuelta] = useState(0);
  const guardar = (b: TrabajoBlast) => setP((q) => ({ ...q, blast: b }));

  async function enviar() {
    setEnviando(true);
    const lista = excluir.split(",").map((s) => s.trim()).filter(Boolean);
    try {
      const r = await post<{ rid: string }>("ncbi/blast", { oligos: oligos.map((o) => ({ nombre: o.nombre, seq: o.seq })), excluir: lista });
      guardar({ rid: r.rid, enviado: Date.now(), excluir: lista.join(", "), estado: "esperando" });
    } catch (e) { guardar({ rid: "", enviado: Date.now(), excluir: "", estado: "error", error: (e as Error).message }); }
    setEnviando(false);
  }

  // NCBI pide no preguntar por la misma búsqueda más de una vez por minuto.
  useEffect(() => {
    if (t?.estado !== "esperando") return;
    const id = setTimeout(async () => {
      try {
        const r = await post<{ estado: string; sitios: Sitio[]; productos: Producto[] }>("ncbi/blast/resultado", { rid: t.rid, pares: pares(p.oligos) });
        if (r.estado === "listo") guardar({ ...t, estado: "listo", sitios: r.sitios, productos: r.productos });
      } catch (e) { guardar({ ...t, estado: "error", error: (e as Error).message }); }
      setVuelta((v) => v + 1);
    }, 60000);
    return () => clearTimeout(id);
  }, [t?.rid, t?.estado, vuelta]); // eslint-disable-line react-hooks/exhaustive-deps

  const porPrimer = (nombre: string) => {
    const s = (t?.sitios ?? []).filter((x) => x.primer === nombre);
    const especiesHit = [...new Set(s.map((x) => x.especie))];
    return { n: s.length, especies: especiesHit };
  };
  const minutos = t ? Math.round((Date.now() - t.enviado) / 60000) : 0;

  return (
    <>
      <h3>6 · Especificidad en NCBI (BLAST contra toda la base nucleotídica)</h3>
      <div className="fila">
        <label style={{ flex: 1 }}>Especies a excluir (las del panel, que ya revisa la PCR virtual){" "}
          <input style={{ width: "100%" }} value={excluir} onChange={(e) => setExcluir(e.target.value)} />
        </label>
        <button className="primario" disabled={!oligos.length || enviando || t?.estado === "esperando"} onClick={enviar}>
          {enviando ? "Enviando…" : `Buscar ${oligos.length} primers en NCBI`}
        </button>
      </div>
      {especies.some((e) => /^\w\. /.test(e)) && (
        <p className="caja amarillo">Escribe los nombres de especie completos (p. ej. «Staphylococcus aureus», no «S. aureus»); NCBI no reconoce las abreviaturas.</p>
      )}
      {t?.estado === "esperando" && (
        <p className="caja">Búsqueda enviada hace {minutos} min (RID <code>{t.rid}</code>). NCBI suele tardar entre 5 y 20 min con primers cortos;
          se revisa sola cada minuto. Puedes cerrar la página y volver: NCBI guarda el resultado unas 24 h.</p>
      )}
      {t?.estado === "error" && <p className="error">{t.error}</p>}
      {t?.estado === "listo" && (
        <>
          {t.productos?.length ? (
            <div className="desliza">
              <table>
                <thead><tr><th>Par</th><th>Primers</th><th>Especie</th><th>Registro</th><th>Tamaño (pb)</th><th>Desajustes</th></tr></thead>
                <tbody>
                  {t.productos.map((x, i) => (
                    <tr key={i} className={x.desajustes <= 2 ? "rojo" : "amarillo"}>
                      <td>{x.par}</td><td>{x.primers}</td><td><i>{x.especie}</i></td>
                      <td><a href={`https://www.ncbi.nlm.nih.gov/nuccore/${x.accesion}`} target="_blank" rel="noreferrer" title={x.titulo}>{x.accesion}</a></td>
                      <td>{x.tamano}</td><td>{x.desajustes}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          ) : (
            <p className="caja verde">Ningún par amplifica fuera de las especies excluidas (≤ 3 desajustes, extremo 3′ completo, ≤ 4000 pb).</p>
          )}
          <div className="desliza">
            <table>
              <thead><tr><th>Primer</th><th>Sitios fuera del panel</th><th>Especies donde pega</th></tr></thead>
              <tbody>
                {oligos.map((o) => {
                  const r = porPrimer(o.nombre);
                  return (
                    <tr key={o.nombre}>
                      <td>{o.nombre}</td><td className={r.n ? "amarillo" : "verde"}>{r.n}</td>
                      <td>{r.especies.slice(0, 5).join(", ")}{r.especies.length > 5 && ` y ${r.especies.length - 5} más`}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
          <p className="nota">Excluidas: {t.excluir || "ninguna"}. Un sitio = el primer pega con su extremo 3′ completo y ≤ 3 desajustes; un producto = dos sitios
            del mismo par mirándose a ≤ 4000 pb en el mismo registro. Es la misma idea que Primer-BLAST.</p>
        </>
      )}
    </>
  );
}
