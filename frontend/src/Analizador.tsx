import { useMemo, useState } from "react";
import type { Ctx } from "./App";
import { Esquema } from "./Alineamiento";
import {
  Condicion, Condiciones, Diana, Estructura, Informe, Nivel, Oligo, PrimerUbicado, csv, descargar, especieDe, fmt, genDe,
  especiesPegadas, pegarPrimers, post, ubicarPrimer,
} from "./lib";
import { escenaEstructura } from "./Vista3D";

type Matriz = { nombres: string[]; matriz: number[][]; pares: { a: string; b: string; dg: number; anclaje_3p: number; nivel: Nivel; estructura: Estructura }[] };
type Fila = Informe & { tubo: string };
type Ubicado = { o: Oligo; d: Diana; u: PrimerUbicado } | { o: Oligo; d?: undefined; u?: undefined };

// Tus dos tablas, con las mismas columnas y en el mismo orden.
const COLUMNAS: [string, keyof Fila, number][] = [
  ["tubo", "tubo", 0], ["oligo", "nombre", 0], ["secuencia 5a3", "seq", 0], ["largo", "largo", 0], ["gc pct", "gc_pct", 1],
  ["Tm IDT valores de fabrica", "tm_fabrica", 1], ["Tm IDT condiciones PCR", "tm_pcr", 1], ["peso molecular g mol", "pm", 1],
  ["coef extincion L mol cm", "eps260", 0], ["nmol por OD260", "nmol_od", 2], ["ug por OD260", "ug_od", 2],
  ["dH kcal mol Entalpía", "dh", 1], ["dS cal K mol Entropía", "ds", 1], ["horquilla dG fabrica 25C", "horquilla_dg", 2],
  ["horquilla Tm", "horquilla_tm", 1], ["autodimero dG fabrica", "autodimero_dg_fabrica", 2], ["autodimero dG PCR", "autodimero_dg", 2],
  ["dG maximo duplex perfecto", "dg_duplex_37", 1],
];
const COLUMNAS_POS = ["Archivo", "Especie", "Gen", "Primer", "Hebra", "Secuencia 5' 3'", "pegar en Find Motif", "inicio", "fin", "largo",
  "coincidencias esperadas", "posiciones variables dentro del primer", "codones", "aminoacidos en esa region"];
const PALETA = ["#84cc16", "#06b6d4", "#7e22ce", "#f97316", "#e11d48", "#2563eb"];

/** Busca el primer en todos los genes alineados y se queda con el de menos desajustes. */
function ubicarEnPanel(o: Oligo, dianas: Diana[]): Ubicado {
  let mejor: Ubicado = { o };
  for (const d of dianas) {
    if (!d.alin) continue;
    const u = ubicarPrimer(o.nombre, o.conjunto, o.seq, d.alin);
    if (u && (!mejor.u || u.desajustes.length < mejor.u.desajustes.length)) mejor = { o, d, u };
  }
  return mejor;
}
const filaPosicion = (x: Ubicado) => x.u
  ? [x.d.archivo, especieDe(x.d), genDe(x.d), x.u.nombre, x.u.hebra, x.u.seq, x.u.motivo, x.u.inicio, x.u.final, x.u.motivo.length,
    x.u.coincidencias, x.u.variables, x.u.codones, x.u.aminoacidos]
  : ["—", "—", "—", x.o.nombre, "no encontrado en los FASTA cargados", x.o.seq, "", "", "", x.o.seq.length, "", "", "", ""];
const tsv = (enc: string[], filas: unknown[][]) => [enc, ...filas].map((f) => f.map((v) => (v == null ? "" : String(v))).join("\t")).join("\n");

export function Analizador({ p, setP, ver3d }: Ctx) {
  const [res, setRes] = useState<Fila[] | null>(null);
  const [matrices, setMatrices] = useState<Record<string, Matriz>>({});
  const [presetDim, setPresetDim] = useState<"pcr" | "multiplex">("pcr");
  const [celda, setCelda] = useState<Matriz["pares"][number] | null>(null);
  const [error, setError] = useState("");
  const [aviso, setAviso] = useState("");
  const [cargando, setCargando] = useState(false);
  const [nFilas, setNFilas] = useState(3);
  const [ancho, setAncho] = useState(100);
  const [colorPor, setColorPor] = useState<"conjunto" | "hebra" | "tubo">("conjunto");
  const [completo, setCompleto] = useState(false);
  const oligos = p.oligos;

  /** Agrega lo nuevo; si un primer ya estaba y lo pegado trae su tubo, se lo actualiza. */
  function agregar(nuevos: Oligo[], especies = new Map<string, string>()) {
    const igual = (a: Oligo, b: Oligo) => a.nombre === b.nombre && a.seq === b.seq;
    let movidos = 0;
    const lista = oligos.map((o) => {
      const n = nuevos.find((x) => igual(x, o));
      if (!n?.tubo || n.tubo === o.tubo) return o;
      movidos++;
      return { ...o, tubo: n.tubo };
    });
    const altas = nuevos.filter((n) => !oligos.some((o) => igual(o, n))).map((n) => ({ ...n, tubo: n.tubo || p.tubos[0] }));
    const dianas = p.dianas.map((d) => (especies.get(d.archivo) && !d.especie ? { ...d, especie: especies.get(d.archivo) } : d));
    setP({ ...p, dianas, oligos: [...lista, ...altas], tubos: [...new Set([...p.tubos, ...altas.map((n) => n.tubo)])] });
    setAviso(`Se agregaron ${altas.length} primers${movidos ? ` y se actualizó el tubo de ${movidos}` : ""}. Total: ${lista.length + altas.length}.`);
  }
  const cambiar = (i: number, c: Partial<Oligo>) => setP({ ...p, oligos: oligos.map((o, j) => (j === i ? { ...o, ...c } : o)) });

  async function analizar() {
    const validos = oligos.filter((o) => o.seq.trim());
    const malos = validos.filter((o) => !/^[ACGT]{8,120}$/i.test(o.seq.replace(/\s/g, "")));
    if (!validos.length) return setError("Agrega al menos un primer.");
    if (malos.length) return setError(`Secuencia no válida (sólo A, C, G, T; 8–120 nt): ${malos.map((o) => o.nombre || o.seq).join(", ")}`);
    setCargando(true); setError(""); setCelda(null);
    try {
      const lista = validos.map((o, i) => ({ nombre: o.nombre.trim() || `primer${i + 1}`, seq: o.seq.replace(/\s/g, "").toUpperCase(), tubo: o.tubo }));
      const r = await post<{ oligos: Informe[] }>("analizar", { oligos: lista, condiciones: p.condiciones });
      setRes(r.oligos.map((o, i) => ({ ...o, tubo: lista[i].tubo })));
      const porTubo: Record<string, Matriz> = {};
      for (const t of [...new Set(lista.map((o) => o.tubo))]) {
        porTubo[t] = await post<Matriz>("dimeros", { oligos: lista.filter((o) => o.tubo === t), condiciones: p.condiciones, preset: presetDim });
      }
      setMatrices(porTubo);
    } catch (e) { setError((e as Error).message); }
    setCargando(false);
  }

  // Posiciones y esquemas: cada primer buscado en todos los genes alineados (pasos 1–2).
  const ubicados = useMemo(() => (res ? res.map((r) => ubicarEnPanel({ nombre: r.nombre, seq: r.seq, tubo: r.tubo,
    conjunto: oligos.find((o) => o.nombre === r.nombre)?.conjunto ?? "" }, p.dianas)) : []), [res, p.dianas]); // oligos sólo aporta el conjunto (color)
  const alineadas = p.dianas.some((d) => d.alin);
  const porGen = p.dianas.map((d) => ({ d, us: ubicados.filter((x) => x.d === d).map((x) => ({ ...x.u!, tubo: x.o.tubo })) })).filter((g) => g.us.length);
  const filasControl = res?.map((o) => COLUMNAS.map(([, k]) => o[k] as string | number | null)) ?? [];
  const filasPos = ubicados.map(filaPosicion);

  async function descargarExcel() {
    const r = await fetch("/api/excel", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ hojas: [
        { nombre: "Tabla de control", encabezados: COLUMNAS.map((c) => c[0]), filas: filasControl },
        ...(alineadas ? [{ nombre: "Posiciones", encabezados: COLUMNAS_POS, filas: filasPos }] : []),
      ] }),
    });
    if (!r.ok) return setError("No se pudo generar el Excel.");
    descargar(`${p.nombre}_tablas.xlsx`, await r.blob(), "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");
  }
  async function copiar(texto: string, que: string) {
    try { await navigator.clipboard.writeText(texto); setAviso(`${que} copiada: pégala en Excel con Ctrl+V.`); }
    catch { setError("El navegador no dejó copiar; usa «Descargar Excel»."); }
  }

  const c = p.condiciones;
  const metodos = `Las temperaturas de fusión (Tm) se calcularon por el método de vecino más cercano (Allawi y SantaLucia, 1997; SantaLucia, 1998) con la corrección de sales de Owczarzy et al. (2004) en ausencia de Mg²⁺ y de Owczarzy et al. (2008) en presencia de Mg²⁺, para ${c.pcr.oligo_uM} µM de oligo, ${c.pcr.na_mM} mM de sal monovalente, ${c.pcr.mg_mM} mM de Mg²⁺ y ${c.pcr.dntp_mM} mM de dNTPs (condiciones de PCR) y para ${c.fabrica.oligo_uM} µM de oligo y ${c.fabrica.na_mM} mM de Na⁺ sin Mg²⁺ (valores de referencia de IDT). Las horquillas (25 °C), autodímeros y heterodímeros (37 °C) se evaluaron con el alineamiento termodinámico de Primer3 (primer3-py 2.3.1). Cálculos hechos con PrimerFlow, motor v0.1.0.`;

  return (
    <section>
      <h2>Analizador: llena tus tablas y esquemas</h2>

      <h3>1 · Tus primers</h3>
      <textarea className="pegar" rows={3} value=""
        placeholder={"📋 Pega aquí (Ctrl+V) filas copiadas de Excel —cualquiera de tus tablas, con o sin encabezados— o líneas «nombre secuencia». Se agregan solas a la lista de abajo."}
        onChange={() => undefined}
        onPaste={(e) => {
          e.preventDefault();
          const texto = e.clipboardData.getData("text");
          const nuevos = pegarPrimers(texto, "Tesis");
          nuevos.length ? agregar(nuevos, especiesPegadas(texto)) : setError("No encontré secuencias (A, C, G, T de 10 nt o más) en lo que pegaste.");
        }} />
      <div className="fila">
        <button onClick={() => setP({ ...p, oligos: [...oligos, { nombre: "", seq: "", tubo: p.tubos[0], conjunto: "Tesis" }] })}>+ Agregar primer a mano</button>
        {p.dianas.some((d) => d.elegido) && (
          <button onClick={() => agregar(p.dianas.filter((d) => d.elegido).flatMap((d) => [d.elegido!.F, d.elegido!.R]
            .map((o) => ({ nombre: o.nombre, seq: o.seq, tubo: d.tubo, conjunto: "Diseño" }))))}>Traer los pares del diseño (paso 3)</button>
        )}
        {oligos.length > 0 && <button onClick={() => confirm("¿Quitar todos los primers de la lista?") && setP({ ...p, oligos: [] })}>Vaciar lista</button>}
        {aviso && <span className="nota">{aviso}</span>}
      </div>
      {oligos.length > 0 && (
        <table>
          <thead><tr><th>#</th><th>Nombre</th><th>Secuencia 5′→3′</th><th>nt</th><th>Tubo</th><th>Conjunto (color del esquema)</th><th /></tr></thead>
          <tbody>
            {oligos.map((o, i) => {
              const s = o.seq.replace(/\s/g, "");
              return (
                <tr key={i}>
                  <td>{i + 1}</td>
                  <td><input value={o.nombre} placeholder="Saureus_nuc_cF" onChange={(e) => cambiar(i, { nombre: e.target.value })} /></td>
                  <td className={s && !/^[ACGT]+$/i.test(s) ? "rojo" : ""}>
                    <input className="secuencia" value={o.seq} placeholder="CTGAAGCAAGTGCATTTACGA" onChange={(e) => cambiar(i, { seq: e.target.value.toUpperCase() })} />
                  </td>
                  <td>{s.length || ""}</td>
                  <td><select value={o.tubo} onChange={(e) => cambiar(i, { tubo: e.target.value })}>{p.tubos.map((t) => <option key={t}>{t}</option>)}</select></td>
                  <td><input value={o.conjunto} onChange={(e) => cambiar(i, { conjunto: e.target.value })} /></td>
                  <td><button onClick={() => setP({ ...p, oligos: oligos.filter((_, j) => j !== i) })} aria-label="Quitar">✕</button></td>
                </tr>
              );
            })}
          </tbody>
        </table>
      )}
      <div className="dos">
        <div className="fila">
          <button className="primario" disabled={cargando || !oligos.length} onClick={analizar}>{cargando ? "Calculando…" : `Analizar ${oligos.length} primers y llenar las tablas`}</button>
          <label>Dímeros con preset{" "}
            <select value={presetDim} onChange={(e) => setPresetDim(e.target.value as "pcr" | "multiplex")}>
              <option value="pcr">PCR</option><option value="multiplex">multiplex</option>
            </select>
          </label>
        </div>
        <details><summary>Condiciones de cálculo</summary><EditorCondiciones c={p.condiciones} onChange={(condiciones) => setP({ ...p, condiciones })} /></details>
      </div>
      {error && <p className="error">{error}</p>}

      {res && (
        <>
          <div className="fila">
            <h3>2 · Tabla de control (estilo IDT), llena automáticamente</h3>
            <button onClick={() => copiar(tsv(COLUMNAS.map((x) => x[0]), filasControl), "Tabla de control")}>Copiar para Excel</button>
            <button className="primario" onClick={descargarExcel}>Descargar Excel (las dos tablas)</button>
            <button onClick={() => descargar(`${p.nombre}_pedido_IDT.txt`, res.map((o) => `${o.nombre}\t${o.seq}\t25nm\tSTD`).join("\n"))}>Pedido IDT</button>
          </div>
          <div className="desliza">
            <table>
              <thead><tr>{COLUMNAS.map(([t]) => <th key={t}>{t}</th>)}<th>3D</th></tr></thead>
              <tbody>
                {res.map((o) => (
                  <tr key={o.nombre + o.seq} title={o.banderas.map((b) => `${b.regla}: ${b.valor} (${b.nivel})`).join("\n")}>
                    {COLUMNAS.map(([, k, dec]) => (
                      <td key={k} className={o.niveles[k as string] ?? ""}>
                        {k === "seq" ? <code>{o.seq}</code> : typeof o[k] === "number" ? fmt(o[k] as number, dec) : String(o[k] ?? "—")}
                      </td>
                    ))}
                    <td className="nowrap">
                      <button onClick={() => ver3d(escenaEstructura(o.estructuras.horquilla, `Horquilla de ${o.nombre} (ΔG ${o.horquilla_dg})`))}>horquilla</button>
                      <button onClick={() => ver3d(escenaEstructura(o.estructuras.autodimero, `Autodímero de ${o.nombre} (ΔG ${o.autodimero_dg})`))}>autodímero</button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <p className="nota">
            Color = criterios del perfil «Tesis» (verde bien, amarillo revisar, rojo fuera de rango). Las horquillas las calcula Primer3;
            IDT usa UNAFold y puede dar hasta ±1.5 kcal/mol distinto.
          </p>

          <div className="fila">
            <h3>3 · Posiciones en el alineamiento (estilo MEGA / Find Motif)</h3>
            {alineadas && <button onClick={() => copiar(tsv(COLUMNAS_POS, filasPos), "Tabla de posiciones")}>Copiar para Excel</button>}
          </div>
          {!alineadas ? <p className="caja amarillo">Para llenar esta tabla y los esquemas, carga tus FASTA (paso 1) y alinéalos (paso 2). Cada primer se busca solo en todos los genes.</p> : (
            <div className="desliza">
              <table>
                <thead><tr>{COLUMNAS_POS.map((t) => <th key={t}>{t}</th>)}</tr></thead>
                <tbody>
                  {ubicados.map((x, i) => (
                    <tr key={i} className={x.u ? "" : "rojo"}>
                      {filaPosicion(x).map((v, j) => (
                        <td key={j} className={x.u && j === 10 ? (x.u.coincidencias.split("/")[0] === x.u.coincidencias.split("/")[1] ? "verde" : "amarillo")
                          : x.u && j === 11 ? (x.u.variables ? "amarillo" : "verde") : ""}>
                          {j === 5 || j === 6 ? <code>{v}</code> : String(v)}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {porGen.length > 0 && (
            <>
              <div className="fila">
                <h3>4 · Esquemas por gen (estilo MultAlin)</h3>
                <label>Cepas a mostrar <input type="number" min={1} max={30} value={nFilas} onChange={(e) => setNFilas(Math.max(1, +e.target.value))} /></label>
                <label>Columnas por línea{" "}
                  <select value={ancho} onChange={(e) => setAncho(+e.target.value)}>{[60, 80, 100, 130].map((n) => <option key={n}>{n}</option>)}</select>
                </label>
                <label>Color de las líneas según{" "}
                  <select value={colorPor} onChange={(e) => setColorPor(e.target.value as typeof colorPor)}>
                    <option value="conjunto">conjunto (artículo, convencional, qPCR…)</option>
                    <option value="hebra">hebra (Forward / Reverse)</option>
                    <option value="tubo">tubo</option>
                  </select>
                </label>
                <label><input type="checkbox" checked={completo} onChange={(e) => setCompleto(e.target.checked)} /> gen completo</label>
              </div>
              <p className="nota">
                Cada primer lleva una línea de color hasta donde empieza y, al final, su <b>inicio–fin</b> en el alineamiento (los mismos números de la tabla de posiciones).
                El conjunto sale solo del nombre (<code>Pri_artic</code> → artic, <code>Mpr_QPCR</code> → QPCR) o lo escribes en la lista de arriba; si un gen tiene un solo
                conjunto, se colorea por hebra: Forward <b style={{ color: "#16a34a" }}>verde</b>, Reverse <b style={{ color: "#dc2626" }}>rojo</b>. Los desajustes con el consenso salen en rojo.
              </p>
              {porGen.map(({ d, us }) => {
                const hebra = (u: (typeof us)[number]) => (u.hebra === "sentido" ? "Forward" : "Reverse");
                const porConjunto = colorPor === "conjunto" && new Set(us.map((u) => u.conjunto)).size > 1;
                const primers = us.map((u) => ({ ...u, conjunto: colorPor === "tubo" ? u.tubo : porConjunto ? u.conjunto : hebra(u) }));
                const valores = [...new Set(primers.map((u) => u.conjunto))];
                const colorDe = (v: string) => (v === "Forward" ? "#16a34a" : v === "Reverse" ? "#dc2626" : PALETA[valores.indexOf(v) % PALETA.length]);
                return (
                  <div key={d.id}>
                    <h4><i>{especieDe(d)}</i> · {genDe(d)} <small className="nota">({d.archivo})</small>{" "}
                      {valores.map((v) => <small key={v} className="nota"><span className="punto" style={{ background: colorDe(v) }} />{v} </small>)}
                    </h4>
                    <Esquema alin={d.alin!} primers={primers} ancho={ancho} soloPrimers={!completo} nFilas={nFilas} nombre={d.prefijo} colorDe={colorDe} />
                  </div>
                );
              })}
            </>
          )}

          <h3>5 · Dímeros por tubo</h3>
          {Object.entries(matrices).map(([t, m]) => <MatrizDimeros key={t} tubo={t} m={m} onCelda={setCelda} />)}
          {celda && (
            <div className="tarjeta-grande">
              <b>{celda.a} × {celda.b}</b>: ΔG {celda.dg} kcal/mol · anclaje 3′ {celda.anclaje_3p}
              <button onClick={() => ver3d(escenaEstructura(celda.estructura, `${celda.a} × ${celda.b} (ΔG ${celda.dg})`))}>Ver en 3D</button>
              <pre className="ascii">{celda.estructura.ascii.join("\n") || "Sin estructura estable"}</pre>
            </div>
          )}

          <h3>Párrafo para Métodos</h3>
          <p className="metodos">{metodos}</p>
        </>
      )}
    </section>
  );
}

function MatrizDimeros({ tubo, m, onCelda }: { tubo: string; m: Matriz; onCelda: (c: Matriz["pares"][number]) => void }) {
  const par = (a: string, b: string) => m.pares.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a))!;
  const peores = m.pares.filter((x) => x.a !== x.b).slice(0, 5);
  const fondo = (dg: number) => (dg <= -9 ? "#fca5a5" : dg <= -6 ? "#fde68a" : dg <= -4 ? "#fef9c3" : "#dcfce7");
  return (
    <div>
      <h4>{tubo}: ΔG kcal/mol a 37 °C · {m.nombres.length * (m.nombres.length - 1) / 2} combinaciones</h4>
      <div className="desliza">
        <table className="matriz">
          <thead><tr><th />{m.nombres.map((n) => <th key={n} className="vertical"><span>{n}</span></th>)}</tr></thead>
          <tbody>
            {m.nombres.map((a, i) => (
              <tr key={a}><th>{a}</th>
                {m.nombres.map((b, j) => {
                  const x = par(a, b);
                  return (
                    <td key={b} style={{ background: fondo(m.matriz[i][j]), outline: x.anclaje_3p <= -5 ? "2px solid #dc2626" : undefined }}
                      title={`${a} × ${b}: ${m.matriz[i][j]} · anclaje 3′ ${x.anclaje_3p}`} onClick={() => onCelda(x)}>
                      {m.matriz[i][j].toFixed(1)}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p className="nota">
        Las 5 peores: {peores.map((x) => `${x.a} × ${x.b} (${x.dg})`).join(" · ") || "—"}. Diagonal = autodímero.
        Borde rojo = extremo 3′ anclado (≤ −5). Clic en una celda para ver el apareamiento.
      </p>
    </div>
  );
}

function EditorCondiciones({ c, onChange }: { c: Condiciones; onChange: (c: Condiciones) => void }) {
  const campos: [keyof Condicion, string][] = [["oligo_uM", "Oligo µM"], ["na_mM", "Na⁺ mM"], ["mg_mM", "Mg²⁺ mM"], ["dntp_mM", "dNTPs mM"]];
  return (
    <table className="condiciones">
      <thead><tr><th>Preset</th>{campos.map(([, t]) => <th key={t}>{t}</th>)}</tr></thead>
      <tbody>
        {(Object.keys(c) as (keyof Condiciones)[]).map((k) => (
          <tr key={k}><td>{k === "fabrica" ? "fábrica (IDT)" : k}</td>
            {campos.map(([campo]) => (
              <td key={campo}>
                <input type="number" step="0.05" min={campo === "oligo_uM" ? 0.01 : 0} value={c[k][campo]}
                  onChange={(e) => onChange({ ...c, [k]: { ...c[k], [campo]: Math.max(0, +e.target.value) } })} />
              </td>
            ))}
          </tr>
        ))}
      </tbody>
    </table>
  );
}
