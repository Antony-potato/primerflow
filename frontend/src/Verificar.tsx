import { useState } from "react";
import type { Ctx } from "./App";
import { Celda, PrimerUbicado, Verificacion, csv, descargar, especie, leerTabla, norm, post, ubicarPrimer } from "./lib";

// Tabla de posiciones (M2 / Find Motif): encabezado normalizado → valor calculado.
const POSICIONES: Record<string, (u: PrimerUbicado) => string | number> = {
  "hebra": (u) => u.hebra,
  "pegar en find motif": (u) => u.motivo,
  "inicio": (u) => u.inicio,
  "fin": (u) => u.final,
  "largo": (u) => u.motivo.length,
  "coincidencias esperadas": (u) => u.coincidencias,
  "coincidencias": (u) => u.coincidencias,
  "posiciones variables dentro del primer": (u) => u.variables,
  "codones": (u) => u.codones,
  "aminoacidos en esa region": (u) => u.aminoacidos,
};

export function Verificar({ p, setP }: Ctx) {
  const [res, setRes] = useState<Verificacion | null>(null);
  const [error, setError] = useState("");
  const [cargando, setCargando] = useState(false);
  const texto = p.verificarTexto ?? "";

  function verificarPosiciones(encabezados: string[], filas: string[][]): Verificacion {
    const k = encabezados.map(norm);
    const iArchivo = k.indexOf("archivo"), iNombre = k.findIndex((x) => x === "primer" || x === "oligo");
    const iSeq = k.findIndex((x) => x.startsWith("secuencia"));
    if (iSeq < 0) throw new Error("La tabla necesita una columna de secuencia 5′→3′.");
    const resumen = { ok: 0, aprox: 0, no: 0 };
    const salida = filas.map((f) => {
      const nombre = f[iNombre] ?? "";
      const d = p.dianas.find((x) => x.archivo === f[iArchivo]) ?? p.dianas.find((x) => especie(x.prefijo) === especie(nombre));
      const u = d?.alin ? ubicarPrimer(nombre, "", (f[iSeq] ?? "").toUpperCase().replace(/\s/g, ""), d.alin) : null;
      return encabezados.map((_, i): Celda => {
        const tuyo = f[i] ?? "", calc = POSICIONES[k[i]];
        if (!calc) return { tuyo, estado: "" };
        if (!d?.alin) return { tuyo, estado: "no", nota: d ? `Alinea ${d.prefijo} en el paso 2.` : "No hay una diana cargada para este archivo (paso 1)." };
        if (!u) return { tuyo, estado: "no", nota: "El primer no aparece en el alineamiento (ni su complemento reverso)." };
        const v = calc(u);
        const igual = String(v).toUpperCase() === tuyo.toUpperCase().trim();
        resumen[igual ? "ok" : "no"]++;
        let nota = igual ? undefined : `Calculado: ${v}.`;
        if (!igual && k[i].startsWith("aminoacidos") && String(v).startsWith(tuyo.trim()) && String(v).length === tuyo.trim().length + 1)
          nota = `Calculado: ${v}. Falta el último aminoácido: el rango ${u.codones} son ${String(v).length} codones (el último lo cubre el primer sólo en parte).`;
        return { tuyo, calculado: v, estado: igual ? "ok" : "no", nota };
      });
    });
    return { encabezados, filas: salida, resumen };
  }

  async function verificar() {
    const { encabezados, filas } = leerTabla(texto);
    if (!filas.length) return setError("Pega la tabla con su fila de encabezados (copiada directamente de Excel).");
    setCargando(true); setError(""); setRes(null);
    try {
      const k = encabezados.map(norm);
      if (k.includes("pegar en find motif") || (k.includes("inicio") && k.includes("fin")))
        setRes(verificarPosiciones(encabezados, filas));
      else
        setRes(await post<Verificacion>("verificar-oligos", { encabezados, filas, condiciones: p.condiciones }));
    } catch (e) { setError((e as Error).message); }
    setCargando(false);
  }

  const k = res?.encabezados.map(norm) ?? [];
  const iEtiqueta = Math.max(0, k.findIndex((x) => x === "oligo" || x === "primer"));
  const malas = res?.filas.flatMap((f) => f.map((c, i) => ({ f, c, h: res.encabezados[i] })).filter((x) => x.c.estado === "no")) ?? [];

  return (
    <section>
      <h2>Verificar tablas</h2>
      <p className="nota">
        Pega una tabla copiada de Excel, con su fila de encabezados. Se reconocen dos tipos:{" "}
        <b>oligos</b> (columnas como <code>secuencia 5a3</code>, <code>Tm IDT valores de fabrica</code>, <code>peso molecular g mol</code>…), que se
        recalculan con el motor; y <b>posiciones</b> (<code>Archivo</code>, <code>pegar en Find Motif</code>, <code>inicio</code>, <code>fin</code>,{" "}
        <code>codones</code>…), que se buscan en los alineamientos cargados (pasos 1–2). El orden de las columnas da igual.
      </p>
      <textarea rows={8} value={texto} onChange={(e) => setP({ ...p, verificarTexto: e.target.value })} />
      <div className="fila">
        <button className="primario" disabled={cargando} onClick={verificar}>{cargando ? "Verificando…" : "Verificar"}</button>
        {res && (
          <button onClick={() => descargar(`${p.nombre}_verificacion.csv`, "﻿" + csv([
            ["fila", "columna", "tu valor", "calculado", "diferencia", "estado", "nota"],
            ...res.filas.flatMap((f) => f.flatMap((c, i) => (c.estado ? [[f[iEtiqueta].tuyo, res.encabezados[i], c.tuyo, c.calculado, c.diferencia, c.estado, c.nota]] : []))),
          ]), "text/csv")}>Descargar verificación (CSV)</button>
        )}
      </div>
      {error && <p className="error">{error}</p>}
      {res && (
        <>
          <div className="tarjetas">
            <div className="tarjeta verde"><span>Coinciden</span><b>{res.resumen.ok}</b></div>
            <div className="tarjeta amarillo"><span>Aproximados (modelo distinto)</span><b>{res.resumen.aprox}</b></div>
            <div className="tarjeta rojo"><span>No coinciden</span><b>{res.resumen.no}</b></div>
          </div>
          {malas.length > 0 && (
            <>
              <h3>Lo que no coincide</h3>
              <ul>{malas.map(({ f, c, h }, i) => <li key={i}><b>{f[iEtiqueta].tuyo}</b> · {h}: tienes <code>{c.tuyo}</code>. {c.nota}</li>)}</ul>
            </>
          )}
          <h3>Tabla verificada</h3>
          <p className="nota">Pasa el ratón por una celda para ver el valor calculado y la explicación.</p>
          <div className="desliza">
            <table>
              <thead><tr>{res.encabezados.map((h) => <th key={h}>{h}</th>)}</tr></thead>
              <tbody>
                {res.filas.map((f, i) => (
                  <tr key={i}>
                    {f.map((c, j) => (
                      <td key={j} className={{ ok: "verde", aprox: "amarillo", no: "rojo", "": "" }[c.estado]}
                        title={c.estado ? [`calculado: ${c.calculado ?? "—"}`, c.diferencia != null ? `diferencia: ${c.diferencia}` : "", c.nota ?? ""].filter(Boolean).join("\n") : undefined}>
                        {c.tuyo}{c.estado === "no" && c.calculado != null && <> → <b>{String(c.calculado)}</b></>}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <p className="nota">
            Verde = coincide · amarillo = difiere por un modelo distinto y conocido (Tm de IDT con Mg²⁺; horquillas, que IDT
            calcula con UNAFold) · rojo = no coincide (en rojo se muestra el valor correcto).
          </p>
        </>
      )}
    </section>
  );
}
