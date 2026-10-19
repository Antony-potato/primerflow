import { useMemo, useRef, useState } from "react";
import type { Ctx } from "./App";
import { Alineamiento as Alin, Diana, PrimerUbicado, csv, descargar, leerOligos, post, rango, ubicarPrimer } from "./lib";

const COLORES = ["#84cc16", "#06b6d4", "#7e22ce", "#f97316", "#e11d48", "#2563eb"];

export function Alineamiento({ p, actualizar, d }: Ctx & { d: Diana }) {
  const [alto, setAlto] = useState(100);
  const [bajo, setBajo] = useState(50);
  const [cargando, setCargando] = useState(false);
  const [error, setError] = useState("");
  const [soloPrimers, setSoloPrimers] = useState(true);
  const [nFilas, setNFilas] = useState(99);
  const [ancho, setAncho] = useState(100);

  async function alinear(diana: Diana, archivo?: { texto: string; formato: string }) {
    setCargando(true); setError("");
    try {
      const alin = await post<Alin>("alinear", {
        secuencias: archivo ? [] : diana.secuencias, archivo: archivo?.texto, formato: archivo?.formato ?? "fasta",
        consenso_alto: alto / 100, consenso_bajo: bajo / 100,
      });
      actualizar(diana.id, { alin, candidatos: undefined, elegido: null, region: null });
    } catch (e) { setError(`${diana.prefijo}: ${(e as Error).message}`); }
    setCargando(false);
  }

  async function importar(f: File) {
    alinear(d, { texto: await f.text(), formato: f.name.toLowerCase().endsWith(".msf") ? "msf" : "fasta" });
  }

  const a = d.alin;
  const el = d.elegido;
  // Primers a comparar: lo escrito en la caja; si está vacía, el par elegido en el paso 3.
  const texto = d.comparar ?? "";
  const lista = texto.trim() ? leerOligos(texto) : el ? [
    { tubo: "Diseño", nombre: el.F.nombre, seq: el.F.seq }, { tubo: "Diseño", nombre: el.R.nombre, seq: el.R.seq },
  ] : [];
  const conjuntos = [...new Set(lista.map((x) => x.tubo))];
  const ubicados = useMemo(() => (a ? lista.map((x) => ({ x, u: ubicarPrimer(x.nombre, x.tubo, x.seq, a) })) : []),
    [a, texto, el?.id]); // la lista se deriva de texto y del par elegido
  const primers = ubicados.flatMap(({ u }) => (u ? [u] : []));
  const perdidos = ubicados.filter(({ u }) => !u).map(({ x }) => x.nombre);

  function exportarTabla() {
    descargar(`${d.prefijo}_posiciones_primers.csv`, "﻿" + csv([
      ["Archivo", "Conjunto", "Primer", "Hebra", "Secuencia 5' 3'", "pegar en Find Motif", "inicio", "fin", "largo",
        "coincidencias", "posiciones variables dentro del primer", "desajustes vs consenso", "codones", "aminoacidos en esa region"],
      ...primers.map((u) => [d.archivo, u.conjunto, u.nombre, u.hebra, u.seq, u.motivo, u.inicio, u.final, u.motivo.length,
        u.coincidencias, u.variables, u.desajustes.length, u.codones, u.aminoacidos]),
    ]), "text/csv");
  }

  return (
    <section>
      <h2>Alineamiento y conservación · {d.prefijo}</h2>
      <div className="fila">
        <button className="primario" disabled={cargando} onClick={() => alinear(d)}>{a ? "Re-alinear" : "Alinear"}</button>
        <button disabled={cargando} onClick={async () => { for (const x of p.dianas.filter((x) => !x.alin)) await alinear(x); }}>
          Alinear todas las pendientes
        </button>
        <label className="boton">Importar de MultAlin (.mul / .msf)
          <input type="file" accept=".mul,.msf,.fasta,.fa,.txt" hidden onChange={(e) => e.target.files?.[0] && importar(e.target.files[0])} />
        </label>
        <button onClick={() => descargar(`${d.prefijo}_para_multalin.fasta`, d.secuencias.map((s) => `>${s.id}\n${s.seq}`).join("\n"))}>
          Exportar FASTA para MultAlin
        </button>
        <a href="http://multalin.toulouse.inra.fr/multalin/" target="_blank" rel="noreferrer">Abrir MultAlin ↗</a>
        <label>Consenso alto <input type="number" min={50} max={100} value={alto} onChange={(e) => setAlto(+e.target.value)} />%</label>
        <label>bajo <input type="number" min={10} max={100} value={bajo} onChange={(e) => setBajo(+e.target.value)} />%</label>
        {cargando && <span>Alineando…</span>}
      </div>
      {error && <p className="error">{error}</p>}
      {!a ? <p className="nota">Si todas las secuencias miden lo mismo se comparan columna a columna; si no, se alinean con FAMSA.</p> : (
        <>
          <div className="tarjetas">
            {Object.entries(a.resumen).map(([k, v]) => (
              <div key={k} className="tarjeta"><span>{({ secuencias: "Secuencias", total: "Sitios", con_hueco: "Con hueco", analizables: "Analizables", C: "Conservados C", V: "Variables V", Pi: "Parsimonia Pi", S: "Singletons S" } as Record<string, string>)[k]}</span><b>{v}</b></div>
            ))}
          </div>
          <h3>Variabilidad a lo largo del gen</h3>
          <Pista alin={a} />
          <p className="nota">
            Barras rojas = posiciones variables (más alta = menos conservada). Verde = bloques 100 % conservados ≥ 18 nt
            ({a.bloques.length}; el más largo: {a.bloques[0]?.largo ?? 0} nt).
          </p>

          <h3>Esquema estilo MultAlin para comparar primers</h3>
          <div className="dos">
            <div>
              <textarea rows={7} value={texto} onChange={(e) => actualizar(d.id, { comparar: e.target.value })}
                placeholder={"# Artículo\nPri_artic GGATATCGTCTGGGACTTCCG\nRe_artic GCGGAGCCAGACCGAATTT\n# Convencional\nMipri_Conv GCGTGTGAAATCGATCAGTG\nMirev_Conv TTCAGGGAGGCTTTGTAGTC"} />
              <div className="fila">
                {el && <button onClick={() => actualizar(d.id, { comparar: `${texto.trim()}\n# Diseño\n${el.F.nombre} ${el.F.seq}\n${el.R.nombre} ${el.R.seq}`.trim() })}>+ par elegido</button>}
                <label><input type="checkbox" checked={soloPrimers} onChange={(e) => setSoloPrimers(e.target.checked)} /> sólo bloques con primers</label>
                <label>Cepas a mostrar <input type="number" min={1} max={a.filas.length} value={Math.min(nFilas, a.filas.length)} onChange={(e) => setNFilas(+e.target.value)} /></label>
                <label>Columnas por línea{" "}
                  <select value={ancho} onChange={(e) => setAncho(+e.target.value)}>{[60, 80, 100, 130].map((n) => <option key={n}>{n}</option>)}</select>
                </label>
              </div>
            </div>
            <p className="nota">
              Un primer por línea (<code>nombre secuencia 5′→3′</code>); <code># Conjunto</code> agrupa y da color (artículo, convencional, qPCR…).
              Los reverse se ubican solos por su complemento reverso, como en Find Motif. Se aceptan hasta 4 desajustes con el
              consenso: salen en <b style={{ color: "#dc2626" }}>rojo</b>.
              {perdidos.length > 0 && <><br /><b>No encontrados en este gen:</b> {perdidos.join(", ")}</>}
            </p>
          </div>
          <Esquema alin={a} primers={primers} colorDe={(c) => COLORES[conjuntos.indexOf(c) % COLORES.length]}
            ancho={ancho} soloPrimers={soloPrimers && primers.length > 0} nFilas={nFilas} nombre={d.prefijo} />

          {primers.length > 0 && (
            <>
              <div className="fila"><h3>Posición de cada primer</h3><button onClick={exportarTabla}>Descargar CSV</button></div>
              <div className="desliza">
                <table>
                  <thead><tr>
                    <th>Conjunto</th><th>Primer</th><th>Hebra</th><th>Secuencia 5′→3′</th><th>pegar en Find Motif</th><th>inicio</th><th>fin</th>
                    <th>largo</th><th>coincidencias</th><th>posiciones variables</th><th>desajustes vs consenso</th><th>codones</th><th>aminoácidos</th>
                  </tr></thead>
                  <tbody>
                    {primers.map((u) => (
                      <tr key={u.conjunto + u.nombre}>
                        <td><span className="punto" style={{ background: COLORES[conjuntos.indexOf(u.conjunto) % COLORES.length] }} />{u.conjunto}</td>
                        <td>{u.nombre}</td><td>{u.hebra}</td><td><code>{u.seq}</code></td><td><code>{u.motivo}</code></td>
                        <td>{u.inicio}</td><td>{u.final}</td><td>{u.motivo.length}</td>
                        <td className={u.coincidencias.split("/")[0] === u.coincidencias.split("/")[1] ? "verde" : "amarillo"}>{u.coincidencias}</td>
                        <td className={u.variables ? "amarillo" : "verde"}>{u.variables}</td>
                        <td className={u.desajustes.length ? "rojo" : "verde"}>{u.desajustes.length}</td>
                        <td>{u.codones}</td><td><code>{u.aminoacidos}</code></td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </>
          )}
        </>
      )}
    </section>
  );
}

function Pista({ alin }: { alin: Alin }) {
  const L = alin.conservacion.length;
  return (
    <svg className="pista" viewBox={`0 0 ${L} 34`} preserveAspectRatio="none" role="img" aria-label="Conservación por posición">
      <rect width={L} height={34} fill="#f8fafc" />
      {alin.bloques.map((b) => <rect key={b.inicio} x={b.inicio} width={b.largo} y={30} height={4} fill="#16a34a" />)}
      {alin.conservacion.map((c, i) => {
        const h = Math.min(28, 56 * (1 - c)); // 50 % de conservación o menos = barra completa
        return c < 1 && <rect key={i} x={i} width={Math.max(1, L / 800)} y={28 - h} height={h} fill="#dc2626" />;
      })}
    </svg>
  );
}

// Esquema como la salida de MultAlin: Courier, nombres a la derecha, regla |----+----|, y los primers contados
// como secuencias más (igual que si los metieras en MultAlin), así sólo se pintan las columnas donde coinciden.
const CW = 7.2, LH = 14, FS = 12, MONO = "'Courier New', Courier, monospace";
const ROJO = "#d40000", AZUL = "#0033cc", NEGRO = "#000";
const ALTO = 0.9, BAJO = 0.5; // umbrales por defecto de MultAlin

export function Esquema({ alin, primers, colorDe, ancho, soloPrimers, nFilas, nombre }: {
  alin: Alin; primers: PrimerUbicado[]; colorDe: (conjunto: string) => string; ancho: number; soloPrimers: boolean; nFilas: number; nombre: string;
}) {
  const ref = useRef<SVGSVGElement>(null);
  const L = alin.consenso.length;
  const filas = alin.filas.slice(0, nFilas).map((f) => ({ id: f.id.split("|")[0].slice(0, 20), seq: f.seq.toUpperCase() }));
  // Cada primer como una fila más: su motivo en sus columnas, en blanco el resto.
  const filasPrimer = primers.map((u) => {
    const fila = Array<string>(L).fill(" ");
    u.cols.forEach((c, j) => { fila[c] = u.motivo[j].toUpperCase(); });
    return { u, seq: fila.join("") };
  });
  const todas = [...filas, ...filasPrimer].map((f) => f.seq);
  // Consenso MultAlin: MAYÚSCULA si ≥ 90 % de las filas, minúscula si ≥ 50 %, «.» el resto.
  const cons = rango(L).map((c) => {
    const n: Record<string, number> = {};
    for (const s of todas) if (/[A-Z]/.test(s[c])) n[s[c]] = (n[s[c]] ?? 0) + 1;
    const [b, k] = Object.entries(n).sort((p, q) => q[1] - p[1])[0] ?? [".", 0];
    return k >= ALTO * todas.length ? b : k >= BAJO * todas.length ? b.toLowerCase() : ".";
  }).join("");
  const color = (s: string, c: number) =>
    cons[c] === "." || s[c] !== cons[c].toUpperCase() ? NEGRO : cons[c] === s[c] ? ROJO : AZUL;

  const etiquetas = [...filas.map((f) => f.id), ...primers.map((u) => u.nombre), "Consensus"];
  const izq = Math.max(...etiquetas.map((t) => t.length)) * CW + 12;
  const lineas = filas.length + primers.length + 4; // números, regla, filas, primers, consenso, renglón en blanco
  const bloques = rango(Math.ceil(L / ancho)).map((i) => [i * ancho, Math.min((i + 1) * ancho, L)])
    .filter(([a, b]) => !soloPrimers || primers.some((u) => u.cols.some((c) => c >= a && c < b)));
  const W = izq + ancho * CW + 90, H = bloques.length * lineas * LH + 6;
  const x = (c: number, desde: number) => izq + (c - desde) * CW;

  function texto(s: string, desde: number, hasta: number, y: number, colorDe: (c: number) => string) {
    // Tramos del mismo color; los espacios van aparte porque SVG los colapsaría y correría las letras.
    const clave = (c: number) => (s[c] === " " ? " " : colorDe(c));
    const tramos: [number, number][] = [];
    for (let c = desde, ini = desde; c <= hasta; c++) {
      if (c === hasta || clave(c) !== clave(ini)) { tramos.push([ini, c]); ini = c; }
    }
    return (
      <text y={y} fontFamily={MONO} fontSize={FS} fontWeight="bold">
        {tramos.filter(([i, j]) => j > i && s[i] !== " ").map(([i, j]) => (
          <tspan key={i} fill={colorDe(i)} x={rango(j - i, i).map((c) => x(c, desde).toFixed(1)).join(" ")}>{s.slice(i, j)}</tspan>
        ))}
      </text>
    );
  }
  const nombreFila = (t: string, y: number) => (
    <text x={izq - 8} y={y} fontFamily={MONO} fontSize={FS} fontWeight="bold" textAnchor="end">{t}</text>
  );

  function exportar(formato: "svg" | "png") {
    const fuente = new XMLSerializer().serializeToString(ref.current!);
    if (formato === "svg") return descargar(`${nombre}_esquema.svg`, fuente, "image/svg+xml");
    const img = new Image();
    img.onload = () => {
      const k = 300 / 96, cv = Object.assign(document.createElement("canvas"), { width: W * k, height: H * k });
      const g = cv.getContext("2d")!;
      g.scale(k, k); g.drawImage(img, 0, 0);
      cv.toBlob((bl) => bl && descargar(`${nombre}_esquema.png`, bl, "image/png"));
    };
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(fuente);
  }

  return (
    <>
      <div className="fila">
        <button onClick={() => exportar("svg")}>Descargar SVG</button>
        <button onClick={() => exportar("png")}>Descargar PNG 300 ppp</button>
        <span className="nota">
          Colores de MultAlin con cada primer contado como una secuencia más: <b style={{ color: ROJO }}>rojo</b> ≥ 90 % de las filas,{" "}
          <b style={{ color: AZUL }}>azul</b> ≥ 50 %, negro el resto (en el consenso: MAYÚSCULA, minúscula o «.»).
        </span>
      </div>
      <div className="esquema">
        <svg ref={ref} xmlns="http://www.w3.org/2000/svg" width={W} height={H} viewBox={`0 0 ${W} ${H}`}>
          <rect width={W} height={H} fill="#fff" />
          {bloques.map(([desde, hasta], bi) => {
            const y0 = bi * lineas * LH + LH;
            const yFila = (i: number) => y0 + (i + 2) * LH;
            // Números: inicio del bloque, cada 10 (si no choca con el final) y el final alineado a la derecha.
            const numeros = rango(hasta - desde, desde).filter((c) => c === desde || c === hasta - 1 || ((c + 1) % 10 === 0 && hasta - 1 - c >= 5));
            const regla = " ".repeat(desde) + rango(hasta - desde, desde).map((c) => (c === desde || c === hasta - 1 ? "|" : (c + 1) % 10 === 0 ? "+" : "-")).join("");
            const alFinal = (c: number) => c === hasta - 1 && c !== desde;
            return (
              <g key={desde}>
                {numeros.map((c) => (
                  <text key={c} x={x(alFinal(c) ? c + 1 : c, desde)} y={y0} fontFamily={MONO} fontSize={FS} fontWeight="bold"
                    textAnchor={alFinal(c) ? "end" : "start"}>{c + 1}</text>
                ))}
                {texto(regla, desde, hasta, y0 + LH, () => NEGRO)}
                {filas.map((f, fi) => (
                  <g key={f.id + fi}>{nombreFila(f.id, yFila(fi))}{texto(f.seq, desde, hasta, yFila(fi), (c) => color(f.seq, c))}</g>
                ))}
                {filasPrimer.map(({ u, seq }, pi) => {
                  const y = yFila(filas.length + pi), col = colorDe(u.conjunto), ini = u.cols[0], fin = u.cols.at(-1)!;
                  return (
                    <g key={u.conjunto + u.nombre}>
                      {nombreFila(u.nombre, y)}
                      {ini > desde && ini < hasta && ( // la línea que antes se dibujaba a mano: del inicio del bloque al primer
                        <line x1={x(desde, desde) + 2} x2={x(ini, desde) - 3} y1={y - 4} y2={y - 4} stroke={col} strokeWidth={3.5} strokeLinecap="round" />
                      )}
                      {texto(seq, desde, hasta, y, (c) => color(seq, c))}
                      {fin >= desde && fin < hasta && ( // inicio–fin en columnas del alineamiento, como en la tabla de posiciones
                        <text x={x(fin + 1, desde) + 6} y={y} fontSize={FS - 1} fontWeight="bold" fill={col} fontFamily={MONO}>{u.inicio}–{u.final}</text>
                      )}
                    </g>
                  );
                })}
                {nombreFila("Consensus", yFila(filas.length + primers.length))}
                {texto(cons, desde, hasta, yFila(filas.length + primers.length), (c) => (cons[c] === "." ? NEGRO : cons[c] === cons[c].toUpperCase() ? ROJO : AZUL))}
              </g>
            );
          })}
        </svg>
      </div>
    </>
  );
}
