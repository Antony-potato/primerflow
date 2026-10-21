// Tipos de la API (docs/ARQUITECTURA.md §4.1), llamadas y utilidades pequeñas.

export type Nivel = "verde" | "amarillo" | "rojo";
export type Bandera = { regla: string; campo: string; valor: number; nivel: Nivel };
export type Hebra = { seq: string; k: number[]; colores?: string[] };
export type Estructura = { dg: number; tm?: number | null; ascii: string[]; hebras: Hebra[]; pares: [number, number][]; lazo?: boolean };

export type Informe = {
  nombre: string; seq: string; largo: number; gc_pct: number; tm_fabrica: number; tm_pcr: number;
  pm: number; eps260: number; nmol_od: number; ug_od: number; dh: number; ds: number; dg_duplex_37: number;
  horquilla_dg: number; horquilla_tm: number | null; autodimero_dg_fabrica: number; autodimero_dg: number;
  pinza_gc: number; poly_x: number; estabilidad_3p: number;
  estructuras: { horquilla: Estructura; autodimero: Estructura };
  niveles: Record<string, Nivel>; banderas: Bandera[];
  inicio?: number; fin?: number; hebra?: number;
};

export type Par = {
  id: string; F: Informe; R: Informe; dtm: number; hetero_dg: number; anclaje_3p: number;
  tubo_dg: number | null; tubo_con: string | null; variables: number; variables_3p: number;
  cruzados: number; cruzados_detalle: string[];
  producto: { inicio: number; fin: number; tamano: number; gc_pct: number; seq: string };
  niveles: Record<string, Nivel>; banderas: Bandera[]; rojos: string[];
  desglose: { regla: string; puntos: number }[]; puntaje: number; penalizacion_primer3: number;
  estructuras: { heterodimero: Estructura };
};

export type Alineamiento = {
  filas: { id: string; seq: string }[]; consenso: string; consenso_multalin: string; conservacion_columna: number[];
  plantilla: string; conservacion: number[]; mapa: number[];
  variables: { columna: number; mayoritaria: string; composicion: Record<string, number>; cepas: Record<string, string> }[];
  bloques: { inicio: number; fin: number; largo: number }[];
  resumen: { secuencias: number; total: number; con_hueco: number; analizables: number; C: number; V: number; Pi: number; S: number };
};

export type Secuencia = { id: string; seq: string };
export type Diana = {
  id: string; archivo: string; prefijo: string; tubo: string; secuencias: Secuencia[];
  alin?: Alineamiento; region?: [number, number] | null; tamanos: [number, number];
  candidatos?: Par[]; descartados?: Par[]; elegido?: Par | null; explicacion?: Record<string, string>;
  comparar?: string; // primers a comparar en el esquema ("# Conjunto" + "nombre secuencia")
  especie?: string;  // nombre para las tablas ("Staphylococcus aureus"); por defecto sale del archivo
};
export type Condicion = { oligo_uM: number; na_mM: number; mg_mM: number; dntp_mM: number };
export type Condiciones = { fabrica: Condicion; pcr: Condicion; multiplex: Condicion };
export type Oligo = { nombre: string; seq: string; tubo: string; conjunto: string };
export type Proyecto = {
  nombre: string; tubos: string[]; dianas: Diana[]; oligos: Oligo[]; condiciones: Condiciones;
  verificarTexto?: string; oligosTexto?: string; // oligosTexto: formato viejo, se convierte al abrir
  blast?: import("./Externos").TrabajoBlast; // la búsqueda en NCBI dura minutos: se guarda para retomarla
};
export type Escena = { titulo: string; hebras: Hebra[]; pares: [number, number][]; lazo?: boolean; focos?: { nombre: string; k: number }[] };

export const CONDICIONES: Condiciones = {
  fabrica: { oligo_uM: 0.25, na_mM: 50, mg_mM: 0, dntp_mM: 0 },
  pcr: { oligo_uM: 0.4, na_mM: 50, mg_mM: 2, dntp_mM: 0.2 },
  multiplex: { oligo_uM: 0.2, na_mM: 50, mg_mM: 2.5, dntp_mM: 0.2 },
};

export async function post<T>(ruta: string, cuerpo: unknown): Promise<T> {
  const r = await fetch(`/api/${ruta}`, {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(cuerpo),
  });
  const j = await r.json().catch(() => ({}));
  if (!r.ok) throw new Error(typeof j.detail === "string" ? j.detail : JSON.stringify(j.detail ?? r.statusText));
  return j as T;
}

export function leerFasta(texto: string): Secuencia[] {
  const out: Secuencia[] = [];
  let actual: Secuencia | null = null;
  for (const linea of texto.split(/\r?\n/)) {
    if (linea.startsWith(">")) {
      actual = { id: linea.slice(1).trim().split(/\s+/)[0] || `sec${out.length + 1}`, seq: "" };
      out.push(actual);
    } else if (actual) actual.seq += linea.replace(/[\s\d]/g, "").toUpperCase();
  }
  return out.filter((s) => s.seq);
}

/** '3_S_aureus_nuc_15.fasta' → 'Saureus_nuc' (igual que backend/test_motor.py). */
export function prefijo(archivo: string) {
  let t = archivo.replace(/\.[^.]+$/, "").split("_");
  if (/^\d+[a-z]?$/.test(t[0])) t = t.slice(1);
  return t.length >= 3 ? `${t[0]}${t[1]}_${t[2]}` : archivo.replace(/\.[^.]+$/, "");
}
export const especie = (nombre: string) => nombre.split("_")[0];

const COMP: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", N: "N" };
export const complemento = (s: string) => s.split("").map((b) => COMP[b] ?? "N").join("");
export const rc = (s: string) => complemento(s).split("").reverse().join("");
export const rango = (n: number, desde = 0) => Array.from({ length: n }, (_, i) => i + desde);

export function descargar(nombre: string, contenido: BlobPart, tipo = "text/plain") {
  const url = URL.createObjectURL(new Blob([contenido], { type: tipo }));
  const a = Object.assign(document.createElement("a"), { href: url, download: nombre });
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function csv(filas: (string | number | null | undefined)[][]) {
  return filas.map((f) => f.map((v) => {
    const s = v == null ? "" : String(v);
    return /[",;\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  }).join(",")).join("\n");
}

/** Lista "# Tubo" + "nombre secuencia" → oligos con su tubo. */
export function leerOligos(texto: string) {
  let tubo = "Tubo 1";
  const out: { tubo: string; nombre: string; seq: string }[] = [];
  for (const linea of texto.split(/\r?\n/).map((l) => l.trim()).filter(Boolean)) {
    if (linea.startsWith("#")) { tubo = linea.slice(1).trim() || tubo; continue; }
    const t = linea.split(/\s+/);
    const seq = t.length > 1 ? t.slice(1).join("") : t[0];
    out.push({ tubo, nombre: t.length > 1 ? t[0] : `oligo${out.length + 1}`, seq: seq.toUpperCase() });
  }
  return out;
}

/** Empareja *_cF / *_cR (o _F / _R) por nombre. */
export function emparejar(oligos: { tubo: string; nombre: string; seq: string }[]) {
  const pares = new Map<string, { nombre: string; tubo: string; F?: string; R?: string }>();
  for (const o of oligos) {
    const m = o.nombre.match(/^(.*?)_c?([FR])$/);
    if (!m) continue;
    const p = pares.get(m[1]) ?? { nombre: m[1], tubo: o.tubo };
    p[m[2] as "F" | "R"] = o.seq;
    pares.set(m[1], p);
  }
  return [...pares.values()].filter((p) => p.F && p.R) as { nombre: string; tubo: string; F: string; R: string }[];
}

export const fmt = (v: number | null | undefined, d = 1) => (v == null ? "—" : v.toFixed(d));

// ---- Primers sobre el alineamiento (esquema comparativo y verificación de posiciones) ----

const CODIGO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; // tabla 11 (= estándar fuera del inicio)
export function traducir(s: string) {
  const i = (b: string) => "TCAG".indexOf(b);
  return rango(Math.floor(s.length / 3)).map((c) => {
    const [a, b, d] = [s[3 * c], s[3 * c + 1], s[3 * c + 2]].map(i);
    return a < 0 || b < 0 || d < 0 ? "X" : CODIGO[16 * a + 4 * b + d];
  }).join("");
}

export type PrimerUbicado = {
  nombre: string; conjunto: string; seq: string; hebra: "sentido" | "antisentido"; motivo: string;
  ini: number; fin: number;            // en la plantilla, [ini, fin)
  inicio: number; final: number;       // columnas del alineamiento, base 1 (como MEGA / Find Motif)
  cols: number[]; desajustes: number[]; // columna de cada base · índices que no coinciden con el consenso
  coincidencias: string; variables: number; codones: string; aminoacidos: string;
};

/** Busca el primer (o su complemento reverso) en la plantilla con ≤ maxMm desajustes y lo describe. */
export function ubicarPrimer(nombre: string, conjunto: string, seq: string, a: Alineamiento, maxMm = 4): PrimerUbicado | null {
  const L = a.plantilla;
  let mejor: { i: number; mm: number[]; motivo: string; hebra: PrimerUbicado["hebra"] } | null = null;
  for (const [motivo, hebra] of [[seq, "sentido"], [rc(seq), "antisentido"]] as const) {
    for (let i = 0; i + motivo.length <= L.length; i++) {
      const mm: number[] = [];
      for (let j = 0; j < motivo.length && mm.length <= maxMm; j++) if (motivo[j] !== L[i + j]) mm.push(j);
      if (mm.length <= maxMm && (!mejor || mm.length < mejor.mm.length)) mejor = { i, mm, motivo, hebra };
    }
  }
  if (!mejor) return null;
  const { i, mm, motivo, hebra } = mejor;
  const cols = rango(motivo.length).map((j) => a.mapa[i + j]);
  const exactas = a.filas.filter((f) => cols.every((c, j) => f.seq[c] === motivo[j])).length;
  const c1 = Math.ceil((i + 1) / 3), c2 = Math.ceil((i + motivo.length) / 3);
  return {
    nombre, conjunto, seq, hebra, motivo, ini: i, fin: i + motivo.length, inicio: cols[0] + 1, final: cols.at(-1)! + 1,
    cols, desajustes: mm, coincidencias: `${exactas}/${a.filas.length}`,
    variables: a.conservacion.slice(i, i + motivo.length).filter((c) => c < 1).length,
    codones: `${c1}-${c2}`, aminoacidos: traducir(L.slice((c1 - 1) * 3, c2 * 3)),
  };
}

/** Tabla pegada desde Excel (tabuladores) o desde texto (2+ espacios). Primera línea = encabezados. */
export function leerTabla(texto: string) {
  const lineas = texto.split(/\r?\n/).filter((l) => l.trim());
  const sep = lineas[0]?.includes("\t") ? /\t/ : /\s{2,}/;
  const [enc, ...filas] = lineas.map((l) => l.split(sep).map((x) => x.trim()));
  return { encabezados: enc ?? [], filas };
}
export const norm = (t: string) => t.normalize("NFKD").replace(/[\u0300-\u036f]/g, "").toLowerCase().replace(/[^a-z0-9]+/g, " ").trim();

export type Celda = { tuyo: string; calculado?: string | number | null; diferencia?: number | null; estado: "ok" | "aprox" | "no" | ""; nota?: string };
export type Verificacion = { encabezados: string[]; filas: Celda[][]; resumen: Record<"ok" | "aprox" | "no", number> };

/** '3_S_aureus_nuc_15.fasta' → 'S. aureus' (editable en el paso 1) y gen 'nuc'. */
export function especieDe(d: Diana) {
  if (d.especie) return d.especie;
  let t = d.archivo.replace(/\.[^.]+$/, "").split("_");
  if (/^\d+[a-z]?$/.test(t[0])) t = t.slice(1);
  return t.length >= 2 ? `${t[0]}. ${t[1]}` : d.prefijo;
}
export const genDe = (d: Diana) => d.prefijo.split("_")[1] ?? "";

/** Si lo pegado trae columnas «Archivo» y «Especie» (tu tabla de posiciones): archivo → especie. */
export function especiesPegadas(texto: string) {
  const { encabezados, filas } = leerTabla(texto);
  const k = encabezados.map(norm), ia = k.indexOf("archivo"), ie = k.indexOf("especie");
  return ia < 0 || ie < 0 ? new Map<string, string>() : new Map(filas.map((f) => [f[ia], f[ie]] as [string, string]));
}

/** Lo que pegues (filas de Excel de cualquiera de tus tablas, o «nombre secuencia») → lista de primers. */
export function pegarPrimers(texto: string, conjunto: string): Oligo[] {
  const out: Oligo[] = [];
  let grupo = "", iNombre = -1; // tubo vacío = no venía en lo pegado
  for (const linea of texto.split(/\r?\n/).map((l) => l.trim()).filter(Boolean)) {
    if (linea.startsWith("#")) { grupo = linea.slice(1).trim(); continue; }
    let celdas = linea.split(linea.includes("\t") ? /\t/ : /\s{2,}/).map((c) => c.trim());
    if (celdas.length === 1) celdas = linea.split(/\s+/);
    const iSeq = celdas.findIndex((c) => /^[ACGTU]{10,}$/i.test(c.replace(/\s/g, "")));
    if (iSeq < 0) { // fila de encabezados: ¿dónde está el nombre?
      iNombre = celdas.findIndex((c) => /^(oligo|primer|nombre)$/i.test(c));
      continue;
    }
    const nombre = celdas[iNombre] ?? celdas.find((c, i) => i !== iSeq && /_c?[FR]$/.test(c))
      ?? celdas.find((c, i) => i !== iSeq && !/^tubo/i.test(c) && !/\.(fasta|fa|fas)$/i.test(c)) ?? `primer${out.length + 1}`;
    const t = celdas.find((c) => /^tubo/i.test(c));
    const seq = celdas[iSeq].replace(/\s/g, "").toUpperCase();
    // ponytail: conjunto = sufijo del nombre (Pri_artic → artic) salvo que sea F/R (Saureus_nuc_cF); si no, el que llegó
    const sufijo = nombre.includes("_") ? nombre.split("_").at(-1)! : "";
    const cj = sufijo && !/^c?[FR]\d*$/i.test(sufijo) ? sufijo : conjunto;
    if (!out.some((o) => o.nombre === nombre && o.seq === seq)) out.push({ nombre, seq, tubo: t ?? grupo, conjunto: cj });
  }
  return out;
}
