import { useEffect, useMemo, useState } from "react";
import SeqViz from "seqviz";
import type { Ctx } from "./App";
import { Condiciones, Diana, Informe, Nivel, Par, Proyecto, fmt, post } from "./lib";
import { escenaAmplicon, escenaEstructura } from "./Vista3D";

type Respuesta = { candidatos: Par[]; descartados: Par[]; rangos: [number, number][]; explicacion: Record<string, string> };
type Elegido = { d: Diana; par: Par };

/** Tamaños ocupados en el tubo: ± máx(40 pb, 10 %) para que las bandas se separen en el gel. */
const excluirDe = (otros: Elegido[]) => otros.map(({ par }) => {
  const t = par.producto.tamano, g = Math.max(40, Math.round(t * 0.1));
  return [t - g, t + g];
});

function disenarDiana(d: Diana, otros: Elegido[], condiciones: Condiciones, minCons: number) {
  return post<Respuesta>("disenar", {
    plantilla: d.alin!.plantilla, conservacion: d.alin!.conservacion, prefijo: d.prefijo, region: d.region ?? null,
    tamanos: [d.tamanos], excluir_tamanos: excluirDe(otros), condiciones, min_conservacion: minCons,
    oligos_tubo: otros.flatMap(({ par }) => [par.F, par.R].map((o) => ({ nombre: o.nombre, seq: o.seq }))),
    plantillas_tubo: otros.map(({ d: x }) => ({ nombre: x.prefijo, seq: x.alin!.plantilla })),
  });
}

/** Diseña todos los tubos: en cada uno, la diana más difícil primero (menos nt en bloques 100 % conservados). */
async function disenarPanel(p: Proyecto, actualizar: Ctx["actualizar"], minCons: number) {
  const sinPar: string[] = [];
  const dificultad = (d: Diana) => d.alin!.bloques.reduce((s, b) => s + b.largo, 0);
  for (const tubo of p.tubos) {
    const elegidos: Elegido[] = [];
    for (const d of p.dianas.filter((x) => x.tubo === tubo && x.alin).sort((a, b) => dificultad(a) - dificultad(b))) {
      const r = await disenarDiana(d, elegidos, p.condiciones, minCons);
      const par = r.candidatos[0] ?? null;
      actualizar(d.id, { candidatos: r.candidatos, descartados: r.descartados, explicacion: r.explicacion, elegido: par });
      if (par) elegidos.push({ d, par }); else sinPar.push(d.prefijo);
    }
  }
  return sinPar;
}

export function Diseno({ p, actualizar, ver3d, d }: Ctx & { d: Diana }) {
  const [sel, setSel] = useState<Par | null>(d.elegido ?? null);
  const [minCons, setMinCons] = useState(100);
  const [cargando, setCargando] = useState(false);
  const [error, setError] = useState("");
  const [aviso, setAviso] = useState("");
  const a = d.alin;
  useEffect(() => { if (d.elegido) setSel((s) => s ?? d.elegido!); }, [d.elegido]);

  // Lo ya elegido en el mismo tubo: sus tamaños se evitan y sus oligos entran en dímeros y productos cruzados.
  const otros: Elegido[] = p.dianas.filter((x) => x.id !== d.id && x.tubo === d.tubo && x.elegido && x.alin).map((x) => ({ d: x, par: x.elegido! }));
  const excluir = excluirDe(otros);

  const highlights = useMemo(() => (a ? a.conservacion.flatMap((c, i) => (c < 1 ? [{ start: i, end: i + 1, color: "#fb923c" }] : [])) : []), [a]);
  const primers = useMemo(() => (sel ? [
    { name: sel.F.nombre, start: sel.F.inicio!, end: sel.F.fin!, direction: 1 as const, color: "#16a34a" },
    { name: sel.R.nombre, start: sel.R.inicio!, end: sel.R.fin!, direction: -1 as const, color: "#dc2626" },
  ] : []), [sel]);
  const annotations = useMemo(() => [
    ...(d.region ? [{ name: "región blanco", start: d.region[0], end: d.region[1], direction: 0, color: "#c7d2fe" }] : []),
    ...(sel ? [{ name: `amplicón ${sel.producto.tamano} pb`, start: sel.producto.inicio, end: sel.producto.fin, direction: 1, color: "#fde68a" }] : []),
  ], [d.region, sel]);

  if (!a) return <p className="nota">Alinea esta diana primero (paso 2).</p>;

  async function disenar() {
    setCargando(true); setError("");
    try {
      const r = await disenarDiana(d, otros, p.condiciones, minCons);
      actualizar(d.id, { candidatos: r.candidatos, descartados: r.descartados, explicacion: r.explicacion });
      setSel(r.candidatos[0] ?? null);
    } catch (e) { setError((e as Error).message); }
    setCargando(false);
  }

  async function panel() {
    if (p.dianas.some((x) => !x.alin)) return setError("Alinea todas las dianas primero (paso 2).");
    setCargando(true); setError(""); setAviso("");
    try {
      const sinPar = await disenarPanel(p, actualizar, minCons);
      setAviso(sinPar.length ? `Sin par válido: ${sinPar.join(", ")}. Ábrelas y prueba bajar la conservación o ampliar el rango de producto.`
        : `Panel completo: ${p.dianas.length} pares elegidos, sin dímeros en rojo ni productos cruzados en los tubos.`);
    } catch (e) { setError((e as Error).message); }
    setCargando(false);
  }

  return (
    <section>
      <h2>Diseño · {d.prefijo} <small>({d.tubo})</small></h2>
      <div className="fila">
        <span>Región: <b>{d.region ? `${d.region[0] + 1}–${d.region[1]}` : `gen completo (1–${a.plantilla.length})`}</b></span>
        {d.region && <button onClick={() => actualizar(d.id, { region: null })}>Usar gen completo</button>}
        <label>Producto{" "}
          <input type="number" value={d.tamanos[0]} onChange={(e) => actualizar(d.id, { tamanos: [+e.target.value, d.tamanos[1]] })} />–
          <input type="number" value={d.tamanos[1]} onChange={(e) => actualizar(d.id, { tamanos: [d.tamanos[0], +e.target.value] })} /> pb
        </label>
        <label>Conservación mínima{" "}
          <select value={minCons} onChange={(e) => setMinCons(+e.target.value)}>
            <option value={100}>100 %</option><option value={93}>≥ 93 % (1 cepa distinta de 15)</option><option value={86}>≥ 86 %</option>
          </select>
        </label>
        <button className="primario" disabled={cargando} onClick={disenar}>{cargando ? "Diseñando…" : "Diseñar primers"}</button>
        <button disabled={cargando} onClick={panel} title="Cada tubo: la diana más difícil primero, eligiendo el mejor par de cada una">
          Diseñar panel completo (todas las dianas)
        </button>
      </div>
      {aviso && <p className="caja verde">{aviso}</p>}
      {otros.length > 0 && (
        <p className="nota">
          Multiplex: ya hay {otros.length} par(es) elegidos en {d.tubo}. Se evitan productos de {excluir.map(([x, y]) => `${x}–${y}`).join(", ")} pb
          y se calculan dímeros y productos cruzados contra sus {otros.length * 2} oligos.
        </p>
      )}
      {error && <p className="error">{error}</p>}
      <p className="nota">Arrastra sobre la secuencia para elegir la región blanco. Naranja = posición variable entre cepas.</p>
      <div className="mapa">
        <SeqViz
          name={d.prefijo} seq={a.plantilla} viewer="linear" primers={primers} highlights={highlights} annotations={annotations}
          zoom={{ linear: 45 }} style={{ height: "100%", width: "100%" }}
          onSelection={(s) => {
            if (s.type === "SEQ" && s.start != null && s.end != null && Math.abs(s.end - s.start) >= 60)
              actualizar(d.id, { region: [Math.min(s.start, s.end), Math.max(s.start, s.end)] });
          }}
        />
      </div>

      {d.candidatos && (
        <div className="dos">
          <div>
            <h3>Candidatos ({d.candidatos.length})</h3>
            {d.candidatos.length === 0 && (
              <div className="error">
                {d.descartados?.length ? (
                  <>
                    Primer3 encontró pares, pero todos tienen algún criterio en rojo. Motivos más frecuentes:
                    <ul>{motivos(d.descartados).map(([m, n]) => <li key={m}>{m}: {n} pares</li>)}</ul>
                    {otros.length > 0 && "En multiplex conviene diseñar primero la diana más difícil (la de menos bloques conservados) y luego las demás. "}
                  </>
                ) : (
                  <>
                    Primer3 no encontró pares:
                    <ul>{Object.entries(d.explicacion ?? {}).map(([k, v]) => <li key={k}><code>{k}</code>: {v}</li>)}</ul>
                  </>
                )}
                Prueba ampliar la región o el rango de producto, o bajar la conservación mínima.
              </div>
            )}
            <table className="clic">
              <thead><tr><th>#</th><th>Puntaje</th><th>Forward 5′→3′</th><th>Reverse 5′→3′</th><th>Tm F/R</th><th>pb</th><th>Hetero ΔG</th><th>Avisos</th></tr></thead>
              <tbody>
                {d.candidatos.map((c, i) => {
                  const avisos = c.F.banderas.length + c.R.banderas.length + c.banderas.length;
                  return (
                    <tr key={c.id} className={sel?.id === c.id ? "sel" : ""} onClick={() => setSel(c)}>
                      <td>{i + 1}{d.elegido?.id === c.id && " ✓"}</td>
                      <td><b>{c.puntaje}</b></td>
                      <td><code>{c.F.seq}</code></td>
                      <td><code>{c.R.seq}</code></td>
                      <td>{fmt(c.F.tm_pcr)} / {fmt(c.R.tm_pcr)}</td>
                      <td>{c.producto.tamano}</td>
                      <td className={c.niveles.hetero_dg}>{fmt(c.hetero_dg, 2)}</td>
                      <td className={avisos ? "amarillo" : "verde"}>{avisos || "—"}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
            {!!d.descartados?.length && (
              <details>
                <summary>{d.descartados.length} pares descartados (algún criterio en rojo)</summary>
                <ul>{d.descartados.map((c) => <li key={c.id}><code>{c.F.seq} / {c.R.seq}</code>: {c.rojos.join(" · ")}</li>)}</ul>
              </details>
            )}
          </div>
          {sel && (
            <div>
              <h3>Par seleccionado</h3>
              <div className="fila">
                <button className="primario" onClick={() => actualizar(d.id, { elegido: sel })}>{d.elegido?.id === sel.id ? "Elegido ✓" : "Elegir este par"}</button>
                <button onClick={() => ver3d(escenaAmplicon(sel))}>Amplicón en 3D</button>
                <button onClick={() => ver3d(escenaEstructura(sel.estructuras.heterodimero, `Heterodímero ${sel.F.nombre} × ${sel.R.nombre}`))}>Heterodímero 3D</button>
              </div>
              <FichaPar par={sel} />
            </div>
          )}
        </div>
      )}
    </section>
  );
}

function motivos(descartados: Par[]) {
  const n = new Map<string, number>();
  for (const r of descartados.flatMap((c) => c.rojos.map((x) => x.split(":")[0]))) n.set(r, (n.get(r) ?? 0) + 1);
  return [...n].sort((a, b) => b[1] - a[1]).slice(0, 4);
}

const FILAS: [string, keyof Informe, number][] = [
  ["Posición", "inicio", 0], ["Largo (nt)", "largo", 0], ["%GC", "gc_pct", 1], ["Tm fábrica (°C)", "tm_fabrica", 1], ["Tm PCR (°C)", "tm_pcr", 1],
  ["Pinza GC", "pinza_gc", 0], ["Repetición máx.", "poly_x", 0], ["Horquilla ΔG", "horquilla_dg", 2], ["Horquilla Tm", "horquilla_tm", 1],
  ["Autodímero ΔG PCR", "autodimero_dg", 2], ["Estabilidad 3′", "estabilidad_3p", 2],
];

function FichaPar({ par }: { par: Par }) {
  const celda = (o: Informe, campo: keyof Informe, dec: number) =>
    campo === "inicio" ? `${o.inicio! + 1}–${o.fin}` : fmt(o[campo] as number | null, dec);
  const nivel = (n?: Nivel) => n ?? "";
  return (
    <>
      <table>
        <thead><tr><th /><th>{par.F.nombre}</th><th>{par.R.nombre}</th></tr></thead>
        <tbody>
          {FILAS.map(([t, campo, dec]) => (
            <tr key={campo}><td>{t}</td>
              <td className={nivel(par.F.niveles[campo])}>{celda(par.F, campo, dec)}</td>
              <td className={nivel(par.R.niveles[campo])}>{celda(par.R, campo, dec)}</td>
            </tr>
          ))}
          <tr><td>ΔTm F–R</td><td colSpan={2} className={nivel(par.niveles.dtm)}>{fmt(par.dtm)} °C</td></tr>
          <tr><td>Heterodímero F×R</td><td colSpan={2} className={nivel(par.niveles.hetero_dg)}>{fmt(par.hetero_dg, 2)} kcal/mol · anclaje 3′ {fmt(par.anclaje_3p, 2)}</td></tr>
          {par.tubo_con && (
            <tr><td>Productos cruzados en el tubo</td>
              <td colSpan={2} className={nivel(par.niveles.cruzados)}>{par.cruzados ? par.cruzados_detalle.join(" · ") : "ninguno"}</td>
            </tr>
          )}
          {par.tubo_con && <tr><td>Peor dímero con el tubo</td><td colSpan={2} className={nivel(par.niveles.tubo_dg)}>{fmt(par.tubo_dg, 2)} con {par.tubo_con}</td></tr>}
          <tr><td>Producto</td><td colSpan={2}>{par.producto.tamano} pb · {fmt(par.producto.gc_pct)} % GC · posiciones variables en primers: {par.variables}</td></tr>
        </tbody>
      </table>
      <h4>Puntaje {par.puntaje} / 100</h4>
      <ul className="desglose">{par.desglose.map((x) => <li key={x.regla}>{x.regla}: {x.puntos}</li>)}</ul>
      {[...par.F.banderas.map((b) => ["F", b] as const), ...par.R.banderas.map((b) => ["R", b] as const), ...par.banderas.map((b) => ["Par", b] as const)].map(([q, b], i) => (
        <div key={i} className={`aviso ${b.nivel}`}>{q} · {b.regla}: {b.valor}</div>
      ))}
      <pre className="ascii">{par.estructuras.heterodimero.ascii.join("\n") || "Sin heterodímero estable"}</pre>
    </>
  );
}
