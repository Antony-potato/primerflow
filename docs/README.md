# PrimerFlow — diseño de primers desde cero

> Nombre de trabajo. Estado: **MVP funcionando** · 2026-10-18

Una app web, visual e interactiva, que automatiza el flujo completo de la tesis:
**FASTA multi-cepa → región conservada → primers F/R → validación termodinámica (estilo IDT) →
multiplex por tubo → PCR virtual → amplicones, gel y tabla de control lista para Métodos.**

Hoy ese flujo pasa por 8 herramientas y se copia a mano a Excel (NCBI, MEGA, MultAlin, Primer3Plus,
IDT OligoAnalyzer, Primer-BLAST, Benchling, Excel). La app lo hace en una sola pantalla guiada
y deja cada número con las condiciones con que se calculó.

## Cómo usarlo

Una vez (instalar):

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
cd frontend && npm install && npm run build && cd ..
```

Cada vez (abrir la app en http://127.0.0.1:8765):

```bash
.venv/Scripts/python -m uvicorn app:app --app-dir backend --port 8765
```

Pruebas contra tu Excel: `.venv/Scripts/python backend/test_motor.py`.
Para tocar la interfaz con recarga en vivo: deja el servidor corriendo y en otra terminal `cd frontend && npm run dev`
(http://localhost:5173).

En internet (Vercel, plan gratis): `npx vercel deploy --prod`. `vercel.json` compila el frontend como sitio estático y
`api/index.py` sirve la API de `backend/` como función de Python en `/api/*`.

Flujo en 1 minuto: **1** arrastra los 9 FASTA y asigna tubo → **2** "Alinear todas" → **3** "Diseñar panel completo"
→ **4** "Cargar pares elegidos" + "Analizar" (tabla IDT, dímeros, CSV, Métodos) → **5** "Correr PCR virtual" → **6** 3D.

## Qué hace el MVP (comprobado de punta a punta en Edge con tus 9 FASTA)

| Paso | Hace | Resultado con tu panel |
|---|---|---|
| 1 Secuencias | Multi-FASTA por diana, revisión de codones de inicio/paro y múltiplo de 3, tubos | 134 secuencias, 9 dianas |
| 2 Alineamiento | Por columnas o FAMSA si hay indels; conteos MEGA C/V/Pi/S; bloques conservados; **esquema MultAlin comparativo** (SVG/PNG): pegas conjuntos de primers (`# Artículo`, `# Convencional`, `# qPCR`…) y cada uno sale en su fila con su color, más la tabla de posiciones (Find Motif, coincidencias, variables, codones, aminoácidos); importar/exportar MultAlin | Conteos idénticos a `M2_archivos` (incluido *E. faecalis* con su indel de 24 nt) |
| 3 Diseño | primer3 con máscara de conservación, re-evaluación estilo IDT, puntaje explicado, mapa SeqViz; **panel completo automático** por tubo (la diana más difícil primero), evitando tamaños de banda ocupados, dímeros y **productos cruzados** | 9/9 pares, 0 criterios en rojo, en segundos |
| 4 Analizador | Tabla con tus columnas (guía §6.1), presets editables, matriz de dímeros por tubo, CSV, pedido IDT, párrafo de Métodos | Reproduce `M3_IDT_oligos` y los 73 pares de `M3_IDT_dimeros` |
| 5 PCR virtual | Cada par contra las 134 secuencias, ≤ 2 desajustes, y productos cruzados dentro de cada tubo | Reproduce `M4_PCR_virtual` (15/15, 0 en otras); el panel nuevo: 0 cruzados |
| 6 Vista 3D | Hélice B del amplicón (F verde, R rojo), horquillas y dímeros en 3D | — |
| 7 Verificar tablas | Pegas una tabla de Excel (oligos estilo IDT o posiciones estilo MEGA) y recalcula **cada celda**: verde coincide, amarillo modelo distinto (Tm con Mg²⁺, horquillas UNAFold), rojo no coincide con el valor correcto | Oligos: 220 coinciden, 48 aproximados, 2 erratas. Posiciones: 159 de 162, 3 erratas |

Hallazgos en el camino:
- **2 erratas en tu Excel** (`ug por OD260` = PM / ε × 1000): `Efaecium_ddl_cR` dice 65.37 y debe ser **35.37**;
  `Efaecalis_ddl_cR` dice 30.78 y debe ser **32.78**.
- Diseñando dianas en mal orden, el multiplex puede producir una **banda extra**: en una prueba,
  `Efaecalis_ddl_cF` + `Efaecium_ddl_cR` daban 139 pb en las 15 cepas de *E. faecalis* (el R se pega con 1 desajuste).
  Ahora el diseño lo detecta y lo descarta.
- **3 filas de la tabla de posiciones con un aminoácido de menos**: `Efaecalis_ddl_cF` (71–78: falta G),
  `Abaumannii_basC_cF` (102–109: falta L) y `Abaumannii_basC_cR` (248–255: falta P). Las otras 15 filas sí
  incluyen el último codón que el primer cubre en parte.
- El diseño automático redescubrió solo primers de la tesis (`Saureus_nuc_cF`, `Paeruginosa_ecfX_cF`, `Ecoli_gadAB_cR`).

Pendiente (siguiente iteración): gel virtual y mezcla maestra (M7), exportar GenBank y XLSX, amplicones como
tabla (M6), sonda qPCR, calibrar Tm con Mg²⁺ y horquillas contra la IDT real, descarga NCBI / BLAST / InterPro (F5).

## Documentos

| Archivo | Qué contiene |
|---|---|
| [LOGICA_PRODUCTO.md](LOGICA_PRODUCTO.md) | Flujo paso a paso, pantallas, reglas, criterios, ranking, multiplex, exportaciones. |
| [ARQUITECTURA.md](ARQUITECTURA.md) | Las 3 capas, stack, contrato de la API, estructura de carpetas, proyectos de GitHub y servicios en los que nos apoyamos. |
| [MOTOR_CALCULO.md](MOTOR_CALCULO.md) | Fórmulas exactas, cómo imitamos IDT y otras calculadoras, resultados de validación. |
| [../backend/test_motor.py](../backend/test_motor.py) | Pruebas doradas: el motor reproduce las hojas M2, M3 y M4 de tu Excel. |

## Lo que ya comprobamos (con tus archivos, hoy)

Corrí el motor propuesto (Biopython + primer3-py) contra `tonycheca15_resultados_esperados_por_manual1.xlsx`:

| Métrica | Resultado | Veredicto |
|---|---|---|
| Tm valores de fábrica (18 oligos) | error máx **0.04 °C** | ✅ réplica de IDT |
| Peso molecular, ε260, nmol/OD, µg/OD, ΔH, ΔS | **idénticos** | ✅ |
| Autodímero ΔG (fábrica y PCR) | **idéntico** 18/18 | ✅ |
| Heterodímeros ΔG + anclaje 3′ (73 pares) | **idéntico** 73/73 | ✅ |
| Horquilla ΔG | error máx 0.97 kcal/mol (tolerancia del manual: ±1.5) | ✅ |
| Tm condiciones PCR | error máx 0.74 °C (medio 0.40) | ⚠️ calibrar |
| Horquilla Tm | no coincide (IDT usa UNAFold, que es de pago) | ⚠️ calibrar |
| Conteos MEGA C / V / Pi / S | idénticos en los archivos sin huecos | ✅ |

Conclusión: **no hace falta reinventar las ecuaciones**; las librerías abiertas ya reproducen a IDT
casi en todo. El trabajo real está en el flujo, la visualización y el multiplex.

## Alcance: tus manuales → módulos de la app

| Manual / fase de la guía | Herramienta que se usaba | Módulo en la app | Roadmap |
|---|---|---|---|
| M1 · NCBI | NCBI Nucleotide | **Secuencias**: subir FASTA (descarga por accesión: después) | F2 / F5 |
| M2 · MEGA, Fase 1 | MEGA / ClustalW | **Alineamiento y conservación** | F2 |
| Esquemas de alineamiento | [MultAlin](http://multalin.toulouse.inra.fr/multalin/) | **Esquema estilo MultAlin** (figura con los primers) + importar/exportar `.mul` / `.msf` | F2 |
| Fases 2–3 | a mano + Primer3Plus | **Diseño automático** con máscara de conservación | F2 |
| M3 · IDT, Fase 4 | IDT OligoAnalyzer | **Analizador** (Analyze / Hairpin / Self-Dimer / Hetero-Dimer) | F1 |
| M3 · dímeros | IDT, 73 pares a mano | **Multiplex**: matriz de dímeros por tubo | F3 |
| M4 · Primer-BLAST | NCBI | **PCR virtual** contra el panel + enlace a Primer-BLAST | F3 / F5 |
| M5 · Expasy, UniProt, InterPro | web | **Contexto proteico** (codones, dominios) | F5 |
| M6 · Benchling | Simulate PCR | **Amplicones** + exportar GenBank anotado | F3 |
| M7 · banco | Excel | **Gel virtual, mezcla maestra, programa PCR** | F3 |
| Fase 5 · tabla | Excel | **Tabla de control** con tu esquema de columnas + párrafo de Métodos | F1 |
| — | — | **Vista 3D** (amplicón, dímeros, horquillas) | F4 |

## Roadmap

Cada fase se puede usar sola; no se empieza la siguiente sin cerrar la anterior.

| Fase | Entrega | Listo cuando |
|---|---|---|
| **F0 · Motor** | `backend/` con FastAPI + motor + prueba dorada | `backend/test_motor.py` en verde; Tm PCR calibrada (ver decisiones abiertas) |
| **F1 · Analizador** | Pegas oligos → tabla IDT-like con los dos presets + exportar CSV con tu esquema | Reproduce las 18 filas de `M3_IDT_oligos` en < 2 s (hoy: ~1 h a mano) |
| **F2 · Diseño desde FASTA** | Subir multi-FASTA → alineamiento → mapa SeqViz → candidatos ordenados | Con los 9 FASTA encuentra pares iguales o mejores que los 18 de la tesis |
| **F3 · Multiplex y banco** | Matriz de dímeros por tubo, PCR virtual (134 secuencias), amplicones, gel, mezcla | Reproduce `M3_IDT_dimeros`, `M4_PCR_virtual`, `M6_amplicones`, `M7_*` |
| **F4 · 3D** | Hélice del amplicón con primers, dímeros y horquillas en 3D | Se entiende sin explicación en una presentación |
| **F5 · Extras** | Descarga NCBI, BLAST remoto, proteína/InterPro, modo offline | Según lo que pida la tesis |

## Decisiones tomadas (y por qué)

| Decisión | Por qué |
|---|---|
| **Vite + React + TypeScript**, no Next.js | El backend ya es FastAPI; no necesitamos SSR, rutas de servidor ni SEO. Menos piezas. Si luego se publica en Vercel con login, se migra. |
| **SeqViz** para el mapa | Ya soporta `primers`, `highlights`, `annotations` y `onSelection` (selección de región). MIT. |
| **React Three Fiber**, no 3Dmol.js | No tenemos coordenadas atómicas (PDB); la hélice se genera por fórmula. 3Dmol sólo si algún día hay modelos atómicos. |
| **Python (FastAPI + primer3-py + Biopython)**, no Rust/WASM | Las dos librerías ya reproducen a IDT (validado). Rust no tiene un primer3 maduro; offline queda para F5. |
| **Biopython para la Tm, primer3 para estructuras y diseño** | Biopython reproduce la Tm de IDT; primer3 reproduce exactamente tus dímeros y anclajes 3′. |
| **Sin base de datos** | Un proyecto = un archivo `.json` (guardar/abrir) + autoguardado local. Base de datos cuando haya varios usuarios. |

## Decisiones abiertas (necesito tu respuesta)

1. **¿De dónde salen las columnas `Tm IDT condiciones PCR` y `horquilla Tm` del Excel?** El
   `manual_3_idt_oligoanalyzer.md` dice 60.3 °C para `Saureus_nuc_cF` (lo reproducimos exacto) y el Excel
   dice 61.0. Si el Excel trae lo que da la web real de IDT, necesito 3–5 capturas de pantalla de IDT con el
   preset PCR para calibrar.
2. **dNTPs del preset PCR: ¿0.2 o 0.8 mM?** La guía dice 0.8 (`qPCR Default`); el Manual 3 y la mezcla
   maestra dicen 0.2. La app usará 0.2 (lo que realmente va en el tubo) salvo que digas otra cosa.
3. **Umbral de heterodímero: ¿−6 o −9 kcal/mol?** La guía dice −6; el Manual 3 dice −9 con anclaje 3′ > −5.
   Propuesta: amarillo entre −6 y −9, rojo por debajo de −9 (respeta ambos).
4. **¿Sólo PCR convencional multiplex, o también qPCR con sonda?** (M6 menciona sonda qPCR). Primer3 ya
   diseña la sonda; sólo cambia presets y criterios.
5. **¿Local en tu laptop o web pública?** Cambia el despliegue, no el diseño.
6. **Preset para evaluar los tubos multiplex.** Tu Excel (`M7_mezcla_calculadora`) dice que en multiplex los
   primers van a 0.2 µM y el MgCl₂ a 2.5 mM, pero los dímeros del Manual 3 se calcularon con 0.4 µM / 2 mM.
   Propuesta: preset `multiplex` propio para el Paso 6 (MOTOR_CALCULO.md §1).
7. **Tus pasos:** dijiste que compartirías los pasos. Mapeé tu guía y los 7 manuales; si hay pasos extra
   (o un orden distinto), los incorporo en `LOGICA_PRODUCTO.md §3`.
