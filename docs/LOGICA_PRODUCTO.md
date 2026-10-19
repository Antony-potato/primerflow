# Lógica del producto

## 1. Principios

1. **Ningún número sin sus condiciones.** Toda Tm o ΔG se muestra con su preset (`fábrica` / `PCR`)
   y la temperatura. Nunca "Tm = 60 °C" a secas (tu Manual 3 lo exige para Métodos).
2. **Explicable.** Cada rechazo dice qué regla falló, con el valor y el umbral:
   `Heterodímero −7.1 kcal/mol < −6 (amarillo)`. Nada de "score 73" sin desglose.
3. **Reproducible.** Proyecto `.json` + versión del motor + presets → mismo resultado siempre.
4. **Guiado pero saltable.** Asistente por pasos, pero se puede entrar directo al Analizador pegando oligos.
5. **Las reglas son datos.** Umbrales y presets viven en un perfil editable, no en el código.

## 2. Modelo de datos (conceptual)

```text
Proyecto
├── nombre, autor, fecha, versión del motor
├── condiciones: { fabrica: {oligo_uM 0.25, Na 50, Mg 0, dNTP 0},
│                  pcr:     {oligo_uM 0.4,  Na 50, Mg 2, dNTP 0.2} }
├── criterios:   perfil de umbrales (§4)
├── tubos:       [ {id: "T1", nombre: "Tubo 1 · Gram +"}, ... ]
└── dianas:      [ Diana ]
     Diana
     ├── especie, gen, tubo
     ├── secuencias:   [ {id, etiqueta, seq} ]          ← multi-FASTA
     ├── alineamiento: {filas, consenso, conservación[], variables[], bloques[]}
     ├── region:       {inicio, fin}                     ← elegida en el mapa
     ├── candidatos:   [ Par ]                           ← primer3 + nuestro motor
     └── elegido:      Par | null
          Par: { F: Oligo, R: Oligo, sonda?: Oligo,
                 producto: {inicio, fin, tamaño, gc, seq},
                 puntaje, banderas: [ {regla, valor, umbral, nivel} ] }
          Oligo: { nombre, seq 5'→3', hebra, inicio, fin, informe: InformeOligo }
```

Nomenclatura automática de oligos: `[Especie]_[gen]_cF` / `_cR` (como en tu tabla: `Saureus_nuc_cF`).

## 3. Flujo paso a paso

El asistente tiene un *stepper* a la izquierda. Cada paso muestra su estado (pendiente / listo / con avisos).

### Paso 0 · Proyecto
- **Entrada:** nombre, tubos del multiplex, presets de condiciones, perfil de criterios.
- **Por defecto:** los dos presets de tu Manual 3 y el perfil "Tesis" (§4).
- **Salida:** proyecto vacío, ya guardable.

### Paso 1 · Secuencias (Manual 1)
- **Entrada:** arrastrar uno o varios multi-FASTA (uno por diana). Cada archivo = una diana.
- **La app:** limpia (mayúsculas, quita espacios/números), valida alfabeto `ACGTN`, detecta duplicados,
  longitudes, codón de inicio/paro y si la longitud es múltiplo de 3 (como la hoja `M1_134_secuencias`).
- **Ve el usuario:** tabla de secuencias por diana con alertas (`N` internas, longitudes atípicas, hebra
  probablemente invertida → botón "complemento reverso").
- **Después (F5):** pegar accesión + `from`/`to` + hebra y descargar de NCBI (Entrez).

### Paso 2 · Alineamiento y conservación (Manual 2 / Fase 1)
- **La app:** si todas las secuencias miden lo mismo y no hay indels, alinea por columna (instantáneo);
  si no, corre MSA (FAMSA o MAFFT). Calcula por columna: base mayoritaria, composición, % de conservación
  y clasificación MEGA (**C** conservada, **V** variable, **Pi** parsimonia-informativa, **S** singleton).
  Ya verificado: da los mismos conteos que MEGA en tu Excel.
- **Ve el usuario:**
  - Pista de conservación (barra por columna, roja donde varía) sobre el consenso en SeqViz.
  - Posiciones variables resaltadas en naranja; al pasar el ratón: `T:10 C:5 — SAUR01=C, SAUR05=C...`.
  - Lista de **bloques 100 % conservados ≥ 18 nt** (candidatos naturales a primer), ordenados por largo.
  - Resumen tipo MEGA: sitios totales, con hueco, C / V / Pi / S.
- **Interoperabilidad con MultAlin:** botón "Exportar para MultAlin" (FASTA de entrada) y "Abrir MultAlin".
  Se puede **importar** el resultado de MultAlin (`.mul` = FASTA alineado, o `.msf`) y usarlo como alineamiento
  oficial. Así las coordenadas de la app coinciden con las de las figuras de la tesis.

### Paso 2b · Esquema estilo MultAlin
La figura del alineamiento con los primers encima, para la tesis. Se genera en la app, igual que la de MultAlin:
- **Colores MultAlin:** rojo = consenso alto (≥ 90 %), azul = consenso bajo (≥ 50 %), negro = el resto.
- **Línea de consenso:** mayúscula si ≥ alto, minúscula si ≥ bajo, `.` si no llega.
- **Encima:** caja y flecha del Forward (→) y del sitio del Reverse (←, en la hebra sentido, la misma secuencia
  que se pega en Find Motif), amplicón sombreado, numeración y nombres de cepa.
- Opciones: umbrales alto/bajo, largo de línea (60/80/100), sólo la región del amplicón ± N nt, orden de filas.
- **Cuidado con el 90 %:** con 15 cepas, 14/15 = 93 % ya sale en rojo, aunque **una cepa difiere**.
  Para que el rojo signifique "100 % conservado" (el criterio de los primers), el preset "Tesis" usa
  alto = 100 %. Con el umbral por defecto de MultAlin, la app marca aparte las columnas en rojo que no llegan al 100 %.
- Exporta SVG (editable en Inkscape/Illustrator) y PNG a 300 ppp.

### Paso 3 · Región blanco
- **Entrada:** el usuario arrastra sobre el mapa (`onSelection` de SeqViz) o elige "gen completo".
- **Parámetros del producto:** rango de tamaño (por defecto 150–600 pb convencional; 70–200 qPCR).
- **Multiplex:** si el tubo ya tiene amplicones elegidos, la app **excluye** tamaños demasiado cercanos
  para que las bandas se separen en gel (por defecto: diferencia ≥ máx(40 pb, 10 %)).
  Ejemplo con el Tubo 1 de la tesis (169, 310, 399, 584): no ofrecerá productos de 290–330 pb.

### Paso 4 · Diseño automático (Fases 2–3)
- **La app:** llama a primer3 con el consenso y una **máscara de conservación** (MOTOR_CALCULO.md §8):
  ningún primer pisa una posición variable y los 5 nt del extremo 3′ deben estar 100 % conservados.
  Pide ~50 pares, los re-analiza con el motor estilo IDT (Paso 5), aplica los filtros duros (§4),
  puntúa (§5) y devuelve los mejores 20.
- El Reverse se calcula solo como complemento reverso; el usuario nunca lo hace a mano.
- **Ve el usuario:** tabla de candidatos (rango, F, R, Tm F/R, ΔTm, tamaño, peor ΔG, semáforo) y,
  al seleccionar una fila, los primers dibujados en el mapa (F verde →, R rojo ←) y el amplicón sombreado.
- **Modo manual:** también puede seleccionar un tramo en el mapa y "usar como Forward/Reverse":
  la app calcula el oligo y lo evalúa al instante.

### Paso 5 · Evaluación termodinámica (Manual 3 / Fase 4)
Ficha por oligo con las **mismas pestañas que IDT**, con los dos presets lado a lado:

| Pestaña | Muestra |
|---|---|
| Analyze | largo, %GC, Tm fábrica, Tm PCR, PM, ε260, nmol/OD, µg/OD, ΔH, ΔS |
| Hairpin | ΔG a 25 °C, Tm, dibujo ASCII del plegamiento (como IDT), vista 3D |
| Self-Dimer | ΔG fábrica y PCR, ΔG máximo (dúplex perfecto) como referencia, dibujo, ¿3′ emparejado? |
| Hetero-Dimer | F × R del par (y vs todos los del tubo en el Paso 6) |
| Extras | pinza GC, estabilidad 3′ (ΔG últimos 5 nt), repeticiones, posiciones variables dentro del primer |

Cada valor lleva su semáforo según el perfil de criterios.

### Paso 6 · Multiplex por tubo (Manual 3, dímeros)
- **La app:** para cada tubo, calcula todos los pares de oligos (n·(n−1)/2 + autodímeros) con el preset
  `multiplex` (o `pcr`; README, decisión abierta 6).
  Tubo 1 = 28 combinaciones, Tubo 2 = 45 (como tu manual).
- **Ve el usuario:**
  - **Mapa de calor** oligo × oligo (ΔG, rojo = peor). Clic en una celda → dibujo del apareamiento y 3D.
  - Las **5 peores** parejas (como la tabla del Manual 3) y si emparejan el extremo 3′.
  - Dispersión de Tm del tubo (todos deberían caer en ± 2 °C) y tamaños de amplicón en un gel previo.
  - **Productos cruzados:** F de una diana con R de otra (PCR virtual entre pares del mismo tubo).
- Ver §6 para la estrategia de elección.

### Paso 7 · PCR virtual / especificidad (Manual 4)
- **La app:** cada par contra **todas** las secuencias del proyecto (134 en la tesis), ambas hebras,
  hasta 2 desajustes por primer y extremo 3′ estricto (MOTOR_CALCULO.md §9).
- **Ve el usuario:** matriz par × especie con `15/15` en su especie y `0` en las demás, tamaño obtenido,
  desajustes `F/R`. Debe reproducir `M4_PCR_virtual`.
- **Fuera del panel:** botón "Abrir en Primer-BLAST" con los primers precargados, y (F5) BLASTn del
  amplicón vía la API de NCBI.

### Paso 8 · Contexto proteico (Manual 5) — F5, opcional
Traducción con código genético 11, qué codones/aminoácidos cubre cada primer y en qué dominio InterPro cae.

### Paso 9 · Banco y entrega (Manuales 6–7, Fase 5)
- **Amplicones:** secuencia completa, tamaño, %GC, primeras/últimas 25 nt (como `M6_amplicones`).
- **Gel virtual:** por tubo, agarosa 2 %, marcador 100 pb, posición de cada banda.
- **Mezcla maestra:** nº de reacciones + exceso → µL por componente (como `M7_mezcla_calculadora`).
- **Programa PCR sugerido:** Ta = Tm PCR mínima del tubo − 4 °C (tu programa: Tm ~61 → 57 °C).
- **Exportar** (§8).

## 4. Criterios por defecto (perfil "Tesis", editable)

Fuente: tu guía (G) y tu Manual 3 (M3). Donde se contradicen, el amarillo cubre la zona en disputa.

| Regla | Verde | Amarillo | Rojo (rechazo) | Preset | Fuente |
|---|---|---|---|---|---|
| Largo | 18–24 nt | 16–17 / 25–28 | fuera | — | G |
| %GC | 40–60 % | 35–40 / 60–65 | fuera | — | G |
| Tm fábrica | 50–55 °C | ±2 | fuera | fábrica | G |
| Tm PCR | 58–63 °C | ±1 | fuera | PCR | G |
| ΔTm F–R | ≤ 1 °C | 1–2 | > 2 | PCR | G |
| Pinza GC (G/C en últimos 5 nt) | 1–2 | 3 | 0 o > 3 | — | G |
| Repeticiones (poly-X) | ≤ 3 | 4 | > 4 | — | práctica estándar |
| Horquilla ΔG (25 °C) | > −3 | −3 a −6 | < −6 | fábrica | G, M3 |
| Horquilla Tm | < Ta − 10 | < Ta | ≥ Ta | fábrica | G |
| Autodímero ΔG (37 °C) | > −5 | −5 a −6 | < −6 | PCR | G, M3 |
| Heterodímero ΔG (37 °C) | > −6 | −6 a −9 | < −9 | PCR | G (−6), M3 (−9) |
| Anclaje 3′ en dímeros | > −5 | — | ≤ −5 | PCR | M3 |
| Estabilidad 3′ (ΔG 5 nt) | > −9 | — | ≤ −9 | PCR | M3 |
| Posiciones variables en el primer | 0 | ≥ 1 fuera de los 5 nt 3′ | en los 5 nt 3′ | — | G |
| PCR virtual en su especie | 100 % | ≥ 90 % | < 90 % | — | M4 |
| PCR virtual en otras especies | 0 | — | ≥ 1 | — | M4 |

Perfiles extra listos para elegir: **qPCR (sonda)**, **Estricto multiplex**, **Permisivo (genes muy variables)**.

## 5. Ranking de candidatos

1. **Filtros duros:** cualquier regla en rojo descarta el par (queda visible en "descartados" con el motivo).
2. **Puntaje 0–100** = 100 − Σ penalizaciones, cada una normalizada a su rango verde:
   - distancia de la Tm PCR a 60.5 °C (F y R), ΔTm F–R,
   - distancia del %GC a 50 %,
   - horquilla, autodímero y heterodímero (cuánto se acercan a su umbral),
   - posiciones variables y conservación mínima,
   - penalización propia de primer3 (`PRIMER_PAIR_PENALTY`) como desempate.
3. **Desempate:** menos amarillos → amplicón más cercano al tamaño óptimo → mayor conservación.
4. La ficha de cada par muestra el desglose (qué restó cuántos puntos). Los pesos están en el perfil.

## 6. Multiplex: cómo se elige un panel por tubo

- **Secuencial con conciencia del tubo** *(ponytail: voraz; búsqueda exhaustiva sólo si un tubo no cierra)*:
  1. Diseña la diana más difícil primero (menos bloques conservados).
  2. Para cada diana siguiente, los candidatos se re-puntúan contra **todos los oligos ya elegidos del tubo**
     (heterodímeros, anclaje 3′, productos cruzados, separación de tamaño, ventana de Tm del tubo).
  3. Si ninguna opción pasa, retrocede una diana y prueba su siguiente candidato.
- **Ventana de Tm del tubo:** todos los primers del tubo dentro de ± 2 °C (preset PCR).
- **Separación de bandas:** diferencia ≥ máx(40 pb, 10 %) entre amplicones del mismo tubo.

## 7. Pantallas

### Distribución general

```text
┌──────────────┬────────────────────────────────────────────────┬──────────────────────┐
│ PASOS        │  MAPA (SeqViz)                                  │ INSPECTOR            │
│ 0 Proyecto ✓ │  ▁▁▁▁▂▁▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▃▁▁▁  ← conservación      │ Saureus_nuc_cF       │
│ 1 Secuen.  ✓ │  ATGACAGAATACTTATTAAGTGCTGGCATATG...            │ CTGAAGCAAGTGCATTTACGA│
│ 2 Alinear  ✓ │        ███F███▶            ◀███R███             │ [Analyze][Hairpin]   │
│ 3 Región   ● │        ░░░░░░ amplicón 169 pb ░░░░░░            │ [Self][Hetero][3D]   │
│ 4 Diseño     │                                                 │ Tm fábrica  52.1 °C ●│
│ 5 Evaluar    ├────────────────────────────────────────────────┤ Tm PCR      60.3 °C ●│
│ 6 Multiplex  │  CANDIDATOS                                     │ %GC         42.9 %  ●│
│ 7 PCR virt.  │  #  F              R              Tm   pb  ΔG  │ Horquilla  −1.14    ●│
│ 8 Proteína   │  1  CTGAAGCA...    GCAACTTT...    60.3 169 −3.9│ Autodímero −3.94    ●│
│ 9 Entrega    │  2  ...                                         │ ▸ ver desglose       │
└──────────────┴────────────────────────────────────────────────┴──────────────────────┘
```

### Multiplex (Paso 6)

```text
          Sa_F  Sa_R  Se_F  Se_R  Em_F  Em_R  El_F  El_R
Sa_F      −3.9  ...                                  ← autodímero en la diagonal
Sa_R            −0.1
...                                     −4.7  ← peor del Tubo 1 (amarillo si < −6)
[clic en celda] → dibujo del apareamiento + 3D + "¿3′ emparejado? no"
```

### Analizador directo (F1)
Caja de texto: pegar `nombre<TAB>secuencia` (o sólo secuencias), elegir tubo, botón **Analizar** →
tabla con tu esquema de columnas, semáforos y matriz de dímeros del tubo. Es el primer entregable porque
sustituye la hora de trabajo manual en IDT.

## 8. Exportaciones

| Formato | Para qué |
|---|---|
| **Tabla de control** (CSV / XLSX) | Exactamente las columnas de tu guía §6.1 (`tubo`, `oligo`, `secuencia 5a3`, `largo`, `gc pct`, `Tm IDT valores de fabrica`, ...). |
| **Pedido IDT** (texto masivo) | `nombre⇥secuencia⇥escala⇥purificación`, listo para pegar en el pedido. |
| **FASTA de primers** | Para cualquier otra herramienta. |
| **GenBank anotado** | Secuencia de referencia con primers y amplicón como *features*: se abre en Benchling/SnapGene sin re-anotar (sustituye el Manual 6). |
| **Párrafo de Métodos** | Texto con método, referencias, condiciones y versiones: *"Las Tm se calcularon por vecino más cercano (SantaLucia, 1998) con corrección de sales de Owczarzy (2008), para 0.4 µM de oligo, 50 mM de sal monovalente, 2 mM de Mg²⁺ y 0.2 mM de dNTPs..."* |
| **Proyecto** (`.json`) | Guardar / reabrir / compartir el trabajo completo. |
| **Esquema MultAlin** (SVG/PNG) | Alineamiento coloreado con los primers encima (Paso 2b). |
| **Alineamiento** (`.mul` FASTA / `.msf`) | Ida y vuelta con MultAlin. |
| **Figuras** (PNG/SVG) | Mapa, matriz de dímeros, gel virtual, captura 3D — para la tesis. |

## 9. Criterio de aceptación del producto

La app está terminada (F3) cuando, partiendo **sólo** de los 9 FASTA de esta carpeta, reproduce
automáticamente las hojas del Excel: `M2_archivos`, `M2_primers`, `M3_IDT_oligos`, `M3_IDT_dimeros`,
`M4_PCR_virtual`, `M6_amplicones` y `M7_mezcla_calculadora` — y además propone pares alternativos que
cumplan el perfil "Tesis". Estas hojas se convierten en pruebas automáticas.
