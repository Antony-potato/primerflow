import { useState } from "react";
import type { Ctx } from "./App";
import { emparejar, especie, post } from "./lib";

type Resp = {
  resumen: { par: string; especie: string; tubo: string | null; propia: string; otras: number; probadas_otras: number; tamanos: number[] }[];
  resultados: { par: string; secuencia: string; especie: string; esperado: boolean; producto: boolean; tamano: number | null; desajustes: string | null; combinacion: string | null }[];
  cruzados: { tubo: string; secuencia: string; especie: string; primers: string; tamano: number; desajustes: string }[];
};

export function PcrVirtual({ p }: Ctx) {
  const [mm, setMm] = useState(2);
  const [estricto, setEstricto] = useState(false);
  const [res, setRes] = useState<Resp | null>(null);
  const [error, setError] = useState("");
  const [cargando, setCargando] = useState(false);

  // Pares elegidos en el diseño + pares que se puedan formar con los oligos del analizador (por nombre _cF/_cR).
  const pares = new Map<string, { nombre: string; especie: string; tubo: string; F: string; R: string }>();
  for (const x of emparejar(p.oligos)) pares.set(x.nombre, { ...x, especie: especie(x.nombre) });
  for (const d of p.dianas.filter((d) => d.elegido))
    pares.set(d.prefijo, { nombre: d.prefijo, especie: especie(d.prefijo), tubo: d.tubo, F: d.elegido!.F.seq, R: d.elegido!.R.seq });
  const secuencias = p.dianas.flatMap((d) => d.secuencias.map((s) => ({ ...s, especie: especie(d.prefijo) })));

  async function correr() {
    setCargando(true); setError("");
    try {
      setRes(await post<Resp>("pcr-virtual", { pares: [...pares.values()], secuencias, max_desajustes: mm, estricto_3p: estricto, cruzados: true }));
    } catch (e) { setError((e as Error).message); }
    setCargando(false);
  }

  const raros = res?.resultados.filter((r) => r.producto !== r.esperado) ?? [];
  // Un renglón por (tubo, primers, especie, tamaño): cuántas secuencias lo dan.
  const cruces = new Map<string, Resp["cruzados"][number] & { n: number }>();
  for (const c of res?.cruzados ?? []) {
    const k = [c.tubo, c.primers, c.especie, c.tamano].join("|");
    cruces.set(k, { ...c, n: (cruces.get(k)?.n ?? 0) + 1 });
  }
  return (
    <section>
      <h2>PCR virtual contra el panel</h2>
      <p className="nota">
        {pares.size} pares × {secuencias.length} secuencias, ambas hebras. Los pares salen del paso 3 (elegidos) y del analizador
        (oligos con nombre <code>…_cF</code> / <code>…_cR</code>). La especie es la primera parte del nombre (<code>Saureus</code>).
      </p>
      <div className="fila">
        <label>Desajustes por primer <input type="number" min={0} max={6} value={mm} onChange={(e) => setMm(+e.target.value)} /></label>
        <label><input type="checkbox" checked={estricto} onChange={(e) => setEstricto(e.target.checked)} /> 3′ estricto (0 desajustes en los últimos 3 nt)</label>
        <button className="primario" disabled={cargando || !pares.size || !secuencias.length} onClick={correr}>{cargando ? "Calculando…" : "Correr PCR virtual"}</button>
      </div>
      {error && <p className="error">{error}</p>}
      {res && (
        <>
          <table>
            <thead><tr><th>Par</th><th>Tubo</th><th>Productos en su especie</th><th>Productos en otras especies</th><th>Tamaños (pb)</th></tr></thead>
            <tbody>
              {res.resumen.map((r) => {
                const [a, b] = r.propia.split("/").map(Number);
                return (
                  <tr key={r.par}>
                    <td>{r.par}</td><td>{r.tubo ?? "—"}</td>
                    <td className={a === b && b > 0 ? "verde" : a / Math.max(b, 1) >= 0.9 ? "amarillo" : "rojo"}>{r.propia}</td>
                    <td className={r.otras ? "rojo" : "verde"}>{r.otras} de {r.probadas_otras}</td>
                    <td>{r.tamanos.join(", ") || "—"}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
          <h3>Resultados inesperados ({raros.length})</h3>
          {raros.length === 0 ? <p className="verde caja">Cada par amplifica sólo su especie, en todas sus cepas.</p> : (
            <table>
              <thead><tr><th>Par</th><th>Secuencia</th><th>Especie</th><th>Qué pasó</th><th>pb</th><th>Desajustes</th></tr></thead>
              <tbody>{raros.map((r) => (
                <tr key={r.par + r.secuencia} className="rojo">
                  <td>{r.par}</td><td>{r.secuencia}</td><td>{r.especie}</td>
                  <td>{r.esperado ? "no amplifica su propia cepa" : `amplifica otra especie (${r.combinacion})`}</td><td>{r.tamano ?? "—"}</td><td>{r.desajustes ?? "—"}</td>
                </tr>
              ))}</tbody>
            </table>
          )}
          <h3>Productos cruzados dentro de cada tubo ({cruces.size})</h3>
          {cruces.size === 0 ? <p className="verde caja">Ningún primer de un par forma producto con primers de otro par del mismo tubo.</p> : (
            <>
              <p className="nota">Bandas extra que aparecerían en el multiplex aunque cada par por separado sea específico. Rediseña una de las dos dianas (el paso 3 ya las evita si eliges sus pares en orden).</p>
              <table>
                <thead><tr><th>Tubo</th><th>Primers</th><th>Especie</th><th>Secuencias</th><th>pb</th><th>Desajustes</th></tr></thead>
                <tbody>{[...cruces.values()].map((r, i) => (
                  <tr key={i} className="rojo"><td>{r.tubo}</td><td>{r.primers}</td><td>{r.especie}</td><td>{r.n}</td><td>{r.tamano}</td><td>{r.desajustes}</td></tr>
                ))}</tbody>
              </table>
            </>
          )}
        </>
      )}
    </section>
  );
}
