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

const CW = 8.4, LH = 16, FS = 13;
const ROJO = "#d40000", AZUL = "#0033cc", NEGRO = "#1f2937";

export function Esquema({ alin, primers, colorDe, ancho, soloPrimers, nFilas, nombre }: {
  alin: Alin; primers: PrimerUbicado[]; colorDe: (conjunto: string) => string; ancho: number; soloPrimers: boolean; nFilas: number; nombre: string;
}) {
  const ref = useRef<SVGSVGElement>(null);
  const filas = alin.filas.slice(0, nFilas).map((f) => ({ ...f, id: f.id.split("|")[0].slice(0, 26) })); // etiqueta corta, como MultAlin
  const L = alin.consenso.length;
  const etiquetas = [...filas.map((f) => f.id), ...primers.map((u) => u.nombre)];
  const izq = Math.max(9, ...etiquetas.map((t) => t.length)) * 7.4 + 14;
  const lineas = filas.length + primers.length + 4; // regla, filas, primers, consenso, separación
  const bloques = rango(Math.ceil(L / ancho)).map((i) => [i * ancho, Math.min((i + 1) * ancho, L)])
    .filter(([a, b]) => !soloPrimers || primers.some((u) => u.cols.some((c) => c >= a && c < b)));
  const W = izq + ancho * CW + 50, H = bloques.length * lineas * LH + 10;
  const cons = alin.consenso_multalin;
  const x = (c: number, desde: number) => izq + (c - desde) * CW;
  // Cada primer como una "fila" del alineamiento: su motivo en sus columnas, espacios en el resto.
  const filasPrimer = primers.map((u) => {
    const fila = Array<string>(L).fill(" ");
    const malas = new Set(u.desajustes.map((j) => u.cols[j]));
    u.cols.forEach((c, j) => { fila[c] = u.motivo[j]; });
    return { u, texto: fila.join(""), malas };
  });

  function color(fila: string, c: number) {
    const m = cons[c];
    if (m === "." || fila[c].toUpperCase() !== alin.consenso[c]) return NEGRO;
    return m === m.toUpperCase() ? ROJO : AZUL;
  }

  function texto(s: string, desde: number, hasta: number, y: number, colorDe: (c: number) => string, peso?: string) {
    // Tramos del mismo color; los espacios van aparte porque SVG los colapsaría y correría las letras.
    const clave = (c: number) => (s[c] === " " ? " " : colorDe(c));
    const tramos: [number, number][] = [];
    for (let c = desde, ini = desde; c <= hasta; c++) {
      if (c === hasta || clave(c) !== clave(ini)) { tramos.push([ini, c]); ini = c; }
    }
    return (
      <text y={y} fontFamily="Consolas, 'Courier New', monospace" fontSize={FS} fontWeight={peso}>
        {tramos.filter(([i, j]) => j > i && s[i] !== " ").map(([i, j]) => (
          <tspan key={i} fill={colorDe(i)} x={rango(j - i, i).map((c) => x(c, desde).toFixed(1)).join(" ")}>{s.slice(i, j)}</tspan>
        ))}
      </text>
    );
  }

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
          <b style={{ color: ROJO }}>Rojo</b> ≥ consenso alto · <b style={{ color: AZUL }}>azul</b> ≥ consenso bajo · negro: el resto o distinto del consenso.
        </span>
      </div>
      <div className="esquema">
        <svg ref={ref} xmlns="http://www.w3.org/2000/svg" width={W} height={H} viewBox={`0 0 ${W} ${H}`}>
          <rect width={W} height={H} fill="#fff" />
          {bloques.map(([desde, hasta], bi) => {
            const y0 = bi * lineas * LH + LH;
            const yFila = (i: number) => y0 + (i + 1) * LH;
            return (
              <g key={desde}>
                <text x={x(desde, desde)} y={y0 - 2} fontSize={10} fill="#64748b" fontFamily="sans-serif">{desde + 1}</text>
                {rango(hasta - desde, desde).filter((c) => (c + 1) % 10 === 0).map((c) => (
                  <text key={c} x={x(c, desde) + CW / 2} y={y0 - 2} fontSize={10} textAnchor="middle" fill="#64748b" fontFamily="sans-serif">{c + 1}</text>
                ))}
                <line x1={x(desde, desde)} x2={x(hasta, desde)} y1={y0 + 2} y2={y0 + 2} stroke="#94a3b8" strokeDasharray="1 3" />
                {filas.map((f, fi) => (
                  <g key={f.id}>
                    <text x={4} y={yFila(fi)} fontSize={11} fontFamily="Consolas, monospace" fill="#334155">{f.id}</text>
                    {texto(f.seq, desde, hasta, yFila(fi), (c) => color(f.seq, c))}
                  </g>
                ))}
                {filasPrimer.map(({ u, texto: fila, malas }, pi) => {
                  const y = yFila(filas.length + pi), col = colorDe(u.conjunto), ini = u.cols[0];
                  return (
                    <g key={u.conjunto + u.nombre}>
                      <text x={izq - 6} y={y} fontSize={11} fontFamily="Consolas, monospace" fill={col} textAnchor="end" fontWeight="bold">{u.nombre}</text>
                      {ini >= desde && ini < hasta && ini > desde && (
                        <line x1={x(desde, desde)} x2={x(ini, desde) - 2} y1={y - 4} y2={y - 4} stroke={col} strokeWidth={4} strokeLinecap="round" />
                      )}
                      {texto(fila, desde, hasta, y, (c) => (malas.has(c) ? "#dc2626" : "#111827"), "bold")}
                    </g>
                  );
                })}
                <text x={4} y={yFila(filas.length + primers.length)} fontSize={11} fontFamily="Consolas, monospace" fontWeight="bold" fill="#334155">Consensus</text>
                {texto(cons, desde, hasta, yFila(filas.length + primers.length), (c) => (cons[c] === "." ? NEGRO : cons[c] === cons[c].toUpperCase() ? ROJO : AZUL), "bold")}
                <text x={x(hasta, desde) + 6} y={yFila(0)} fontSize={10} fill="#64748b" fontFamily="sans-serif">{hasta}</text>
              </g>
            );
          })}
        </svg>
      </div>
    </>
  );
}
