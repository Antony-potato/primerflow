import { useEffect, useState } from "react";
import { Alineamiento } from "./Alineamiento";
import { Analizador } from "./Analizador";
import { Diseno } from "./Diseno";
import { CONDICIONES, Diana, Escena, Proyecto, descargar, especieDe, leerFasta, leerOligos, post, prefijo } from "./lib";
import { PcrVirtual } from "./PcrVirtual";
import { Verificar } from "./Verificar";
import { Vista3D } from "./Vista3D";

const CLAVE = "primerflow-proyecto";
const PASOS = ["1 · Secuencias", "2 · Alineamiento", "3 · Diseño", "4 · Analizador y multiplex", "5 · PCR virtual", "6 · Vista 3D", "7 · Revisar mi Excel (opcional)"];
const NUEVO: Proyecto = { nombre: "Panel multiplex", tubos: ["Tubo 1 · Gram +", "Tubo 2 · Gram −"], dianas: [], oligos: [], condiciones: CONDICIONES };

/** Proyectos guardados con el formato viejo (oligos como texto) → lista de oligos. */
function migrar(p: Proyecto): Proyecto {
  const q = { ...NUEVO, ...p };
  if (!p.oligos && p.oligosTexto) q.oligos = leerOligos(p.oligosTexto).map((o) => ({ ...o, conjunto: "Tesis" }));
  return q;
}

function cargar(): Proyecto {
  try {
    const p = JSON.parse(localStorage.getItem(CLAVE) ?? "null");
    if (p?.dianas) return migrar(p);
  } catch { /* sin almacenamiento: proyecto nuevo */ }
  return NUEVO;
}

export type Ctx = {
  p: Proyecto; setP: React.Dispatch<React.SetStateAction<Proyecto>>;
  actualizar: (id: string, cambios: Partial<Diana>) => void; ver3d: (e: Escena) => void;
};

export default function App() {
  const [p, setP] = useState<Proyecto>(cargar);
  const [paso, setPaso] = useState(0);
  const [dianaId, setDianaId] = useState<string | null>(null);
  const [escena, setEscena] = useState<Escena | null>(null);

  useEffect(() => {
    try { localStorage.setItem(CLAVE, JSON.stringify(p)); } catch { /* lleno o bloqueado: se guarda con el botón */ }
  }, [p]);

  const diana = p.dianas.find((d) => d.id === dianaId) ?? p.dianas[0];
  const ctx: Ctx = {
    p, setP,
    actualizar: (id, cambios) => setP((q) => ({ ...q, dianas: q.dianas.map((d) => (d.id === id ? { ...d, ...cambios } : d)) })),
    ver3d: (e) => { setEscena(e); setPaso(5); },
  };
  const hechos = [
    p.dianas.length > 0,
    p.dianas.length > 0 && p.dianas.every((d) => d.alin),
    p.dianas.length > 0 && p.dianas.every((d) => d.elegido),
    false, false, false, false,
  ];

  function abrir(f: File) {
    f.text().then((t) => {
      const q = JSON.parse(t);
      if (!Array.isArray(q?.dianas)) throw new Error("no es un proyecto de PrimerFlow");
      setP(migrar(q));
    }).catch((e) => alert(`No se pudo abrir: ${e.message}`));
  }

  return (
    <div className="app">
      <header>
        <strong>PrimerFlow</strong>
        <input className="titulo" value={p.nombre} onChange={(e) => setP({ ...p, nombre: e.target.value })} aria-label="Nombre del proyecto" />
        {paso >= 1 && paso <= 2 && p.dianas.length > 0 && (
          <label>Diana{" "}
            <select value={diana?.id} onChange={(e) => setDianaId(e.target.value)}>
              {p.dianas.map((d) => <option key={d.id} value={d.id}>{d.prefijo}{d.elegido ? " ✓" : ""}</option>)}
            </select>
          </label>
        )}
        <span className="espacio" />
        <button onClick={() => descargar(`${p.nombre}.primerflow.json`, JSON.stringify(p), "application/json")}>Guardar proyecto</button>
        <label className="boton">Abrir<input type="file" accept=".json" hidden onChange={(e) => e.target.files?.[0] && abrir(e.target.files[0])} /></label>
        <button onClick={() => confirm("¿Borrar todo el proyecto?") && setP(NUEVO)}>Nuevo</button>
      </header>
      <nav>
        {PASOS.map((t, i) => (
          <button key={t} className={i === paso ? "activo" : ""} onClick={() => setPaso(i)}>
            {t} {hechos[i] ? "✓" : ""}
          </button>
        ))}
      </nav>
      <main>
        {paso === 0 && <Secuencias {...ctx} />}
        {paso === 1 && (diana ? <Alineamiento key={diana.id} {...ctx} d={diana} /> : <Vacio />)}
        {paso === 2 && (diana ? <Diseno key={diana.id} {...ctx} d={diana} /> : <Vacio />)}
        {paso === 3 && <Analizador {...ctx} />}
        {paso === 4 && <PcrVirtual {...ctx} />}
        {paso === 5 && <Vista3D key={escena?.titulo} escena={escena} />}
        {paso === 6 && <Verificar {...ctx} />}
      </main>
    </div>
  );
}

const Vacio = () => <p className="nota">Primero sube los FASTA en el paso 1.</p>;

const INICIOS = ["ATG", "TTG", "GTG"], PAROS = ["TAA", "TAG", "TGA"];

function Secuencias({ p, setP, actualizar }: Ctx) {
  const [arrastrando, setArrastrando] = useState(false);

  async function subir(archivos: FileList | File[], especie?: string) {
    const nuevas: Diana[] = [];
    for (const f of Array.from(archivos)) {
      const secuencias = leerFasta(await f.text());
      if (!secuencias.length) { alert(`${f.name}: no tiene secuencias FASTA`); continue; }
      nuevas.push({ id: `${f.name}-${Date.now()}`, archivo: f.name, prefijo: prefijo(f.name), tubo: p.tubos[0], secuencias, tamanos: [150, 600], especie });
    }
    setP((q) => ({ ...q, dianas: [...q.dianas, ...nuevas] }));
  }

  // NCBI: lo que llega se agrega igual que un FASTA arrastrado (el nombre da el prefijo y la especie).
  const [ncbi, setNcbi] = useState({ nombre: "", consulta: "", maximo: 30, trayendo: false, error: "" });
  async function traerNcbi() {
    setNcbi((n) => ({ ...n, trayendo: true, error: "" }));
    try {
      const r = await post<{ fasta: string }>("ncbi/secuencias", { consulta: ncbi.consulta, maximo: ncbi.maximo });
      // especie completa del encabezado («>NC_020547.2:… Acinetobacter baumannii D1279779»): NCBI no entiende abreviaturas
      const especie = r.fasta.match(/^>\S+ ([A-Z][a-z]+ [a-z]+)/)?.[1];
      await subir([new File([r.fasta], `${ncbi.nombre.trim().replace(/\s+/g, "_") || "NCBI"}_ncbi.fasta`)], especie);
      setNcbi((n) => ({ ...n, trayendo: false }));
    } catch (e) { setNcbi((n) => ({ ...n, trayendo: false, error: (e as Error).message })); }
  }

  return (
    <section>
      <h2>Secuencias (una diana por archivo multi-FASTA)</h2>
      <label
        className={`zona ${arrastrando ? "encima" : ""}`}
        onDragOver={(e) => { e.preventDefault(); setArrastrando(true); }}
        onDragLeave={() => setArrastrando(false)}
        onDrop={(e) => { e.preventDefault(); setArrastrando(false); subir(e.dataTransfer.files); }}
      >
        Arrastra aquí tus archivos .fasta (o haz clic)
        <input type="file" multiple accept=".fasta,.fa,.fas,.fna,.txt" hidden onChange={(e) => e.target.files && subir(e.target.files)} />
      </label>

      <details>
        <summary>…o tráelas de NCBI</summary>
        <div className="fila">
          <label>Nombre (especie_gen) <input value={ncbi.nombre} placeholder="A_baumannii_basC" onChange={(e) => setNcbi({ ...ncbi, nombre: e.target.value })} /></label>
          <label>Máximo de secuencias <input type="number" min={1} max={100} value={ncbi.maximo} onChange={(e) => setNcbi({ ...ncbi, maximo: +e.target.value })} /></label>
        </div>
        <textarea rows={4} value={ncbi.consulta} onChange={(e) => setNcbi({ ...ncbi, consulta: e.target.value })}
          placeholder={"Una búsqueda:   \"Acinetobacter baumannii\"[Organism] AND basC[Gene]\n…o accesiones, una por línea, con rango opcional (c = hebra complementaria):\nNC_020547.2:2576156-2577466\nNZ_CP012345.1:c25000-23690"} />
        <div className="fila">
          <button disabled={!ncbi.consulta.trim() || ncbi.trayendo} onClick={traerNcbi}>{ncbi.trayendo ? "Trayendo de NCBI…" : "Traer de NCBI"}</button>
          <span className="nota">Una búsqueda trae registros completos de 100–20 000 nt (p. ej. contigs que contienen el gen), no genomas.
            Para tener sólo el gen de cada cepa, como en tus FASTA, usa accesiones con rango (hasta 40).</span>
        </div>
        {ncbi.error && <p className="error">{ncbi.error}</p>}
      </details>

      <div className="fila">
        <b>Tubos:</b>
        {p.tubos.map((t, i) => (
          <input key={i} value={t} onChange={(e) => setP({ ...p, tubos: p.tubos.map((x, j) => (j === i ? e.target.value : x)) })} />
        ))}
        <button onClick={() => setP({ ...p, tubos: [...p.tubos, `Tubo ${p.tubos.length + 1}`] })}>+ tubo</button>
      </div>

      {p.dianas.length > 0 && (
        <table>
          <thead>
            <tr><th>Archivo</th><th>Prefijo de los oligos</th><th>Especie (para las tablas)</th><th>Tubo</th><th>Secuencias</th><th>Largo (nt)</th><th>Revisión del gen</th><th /></tr>
          </thead>
          <tbody>
            {p.dianas.map((d) => {
              const largos = d.secuencias.map((s) => s.seq.length);
              const sinInicio = d.secuencias.filter((s) => !INICIOS.includes(s.seq.slice(0, 3))).length;
              const sinParo = d.secuencias.filter((s) => !PAROS.includes(s.seq.slice(-3))).length;
              const no3 = d.secuencias.filter((s) => s.seq.length % 3).length;
              const avisos = [sinInicio && `${sinInicio} sin codón de inicio`, sinParo && `${sinParo} sin codón de paro`, no3 && `${no3} no múltiplo de 3`].filter(Boolean);
              return (
                <tr key={d.id}>
                  <td>{d.archivo}</td>
                  <td><input value={d.prefijo} onChange={(e) => actualizar(d.id, { prefijo: e.target.value })} /></td>
                  <td><input value={especieDe(d)} onChange={(e) => actualizar(d.id, { especie: e.target.value })} /></td>
                  <td>
                    <select value={d.tubo} onChange={(e) => actualizar(d.id, { tubo: e.target.value })}>
                      {p.tubos.map((t) => <option key={t}>{t}</option>)}
                    </select>
                  </td>
                  <td>{d.secuencias.length}</td>
                  <td>{Math.min(...largos)}{Math.min(...largos) !== Math.max(...largos) && `–${Math.max(...largos)}`}</td>
                  <td className={avisos.length ? "amarillo" : "verde"}>{avisos.length ? avisos.join(" · ") : "OK"}</td>
                  <td><button onClick={() => setP({ ...p, dianas: p.dianas.filter((x) => x.id !== d.id) })}>Quitar</button></td>
                </tr>
              );
            })}
          </tbody>
        </table>
      )}
      <p className="nota">
        Se revisa codón de inicio (ATG/TTG/GTG), de paro y que el largo sea múltiplo de 3, como en la hoja M1. El prefijo se usa para
        nombrar los oligos (<code>Saureus_nuc_cF</code>) y su primera parte identifica la especie en la PCR virtual.
      </p>
    </section>
  );
}
