# Motor de cálculo

Objetivo: dar los mismos números que IDT OligoAnalyzer (y poder imitar otras calculadoras),
con fórmulas abiertas, citables en Métodos y verificadas contra tu Excel.

## 1. Presets de condiciones

Cada cálculo recibe un preset explícito. Nunca se usan los valores por defecto de una librería
(primer3-py trae por defecto 1.5 mM Mg²⁺ y 0.6 mM dNTPs: si no se pasan, los resultados no cuadran).

| Preset | Oligo | Na⁺ | Mg²⁺ | dNTPs | Uso | Fuente |
|---|---|---|---|---|---|---|
| `fabrica` | 0.25 µM | 50 mM | 0 | 0 | IDT "SpecSheet", horquilla a 25 °C, autodímero fábrica | IDT / Manual 3 |
| `pcr` | 0.4 µM | 50 mM | 2 mM | 0.2 mM | Tm PCR, dímeros PCR (uniplex) | Manual 3 + mezcla maestra |
| `multiplex` | 0.2 µM | 50 mM | 2.5 mM | 0.2 mM | Evaluar los tubos multiplex | nota de `M7_mezcla_calculadora` |
| `primer3` | 50 nM | 50 mM | 1.5 mM | 0.6 mM | Imitar Primer3Plus | defaults de primer3 |

> La guía describe el preset `qPCR (Default)` de IDT con 0.8 mM dNTPs; el Manual 3 usa 0.2 mM.
> Pendiente de confirmar (README, decisión abierta 2).

## 2. Perfil IDT: fórmulas

### 2.1 Tm (ANALYZE)

Vecino más cercano de dos estados con parámetros de **Allawi & SantaLucia 1997** (= SantaLucia 1998
unificado; tabla `DNA_NN3` de Biopython):

$$T_m(1\,M) = \frac{\Delta H^\circ}{\Delta S^\circ + R\,\ln(C_t/4)} - 273.15 \qquad R = 1.987\ \text{cal·K}^{-1}\text{mol}^{-1}$$

- $C_t$ = concentración del oligo. **El factor /4 es el que reproduce a IDT** (validado: con $C_t$ sin
  dividir, la Tm sale ~1.9 °C más alta).
- Corrección de sal sobre $1/T_m$, según $R = \sqrt{[\text{Mg}^{2+}]_{libre}}\,/\,[\text{Mon}^+]$:
  - $R < 0.22$ → **Owczarzy 2004** (monovalente):
    $\frac{1}{T_m} = \frac{1}{T_m(1M)} + (4.29 f_{GC} - 3.95)\times10^{-5}\ln[\text{Na}^+] + 9.40\times10^{-6}\ln^2[\text{Na}^+]$
  - $R \ge 0.22$ → **Owczarzy 2008** (Mg²⁺):
    $\frac{1}{T_m} = \frac{1}{T_m(1M)} + a + b\ln[\text{Mg}] + f_{GC}(c + d\ln[\text{Mg}]) + \frac{1}{2(N-1)}\left(e + f\ln[\text{Mg}] + g\ln^2[\text{Mg}]\right)$
    con $a=3.92\times10^{-5},\ b=-9.11\times10^{-6},\ c=6.26\times10^{-5},\ d=1.42\times10^{-5},\ e=-4.82\times10^{-4},\ f=5.25\times10^{-4},\ g=8.31\times10^{-5}$;
    si $0.22 \le R < 6$, $a$, $d$ y $g$ se ajustan con $[\text{Mon}^+]$ (Owczarzy 2008; igual que Biopython `saltcorr=7`).
- Mg²⁺ libre: los dNTPs quelan Mg²⁺ ($K_a \approx 3\times10^4\ \text{M}^{-1}$). A estas concentraciones,
  restar dNTPs y resolver el equilibrio dan la misma Tm (< 0.01 °C).
- Implementación: `Bio.SeqUtils.MeltingTemp.Tm_NN(seq, nn_table=DNA_NN3, dnac1=Ct_nM/4, dnac2=0, Na=..., Mg=..., dNTPs=..., saltcorr=7)`.

### 2.2 ΔH, ΔS y ΔG del dúplex perfecto

$\Delta H^\circ = \sum \Delta H_{NN} + \Delta H_{init} + \Delta H_{term}(5') + \Delta H_{term}(3')$ (ídem ΔS), con
iniciación terminal A·T (+2.3 kcal, +4.1 eu) o G·C (+0.1 kcal, −2.8 eu).
$\Delta G_{37} = \Delta H^\circ - 310.15\,\Delta S^\circ$ (es el "Maximum ΔG" de la pestaña Self-Dimer de IDT).

### 2.3 Peso molecular, ε260, nmol/OD, µg/OD

| Magnitud | Fórmula | Fuente |
|---|---|---|
| PM anhidro | $\sum n_i\,PM_i - 61.964$ (dA 313.209 · dC 289.184 · dG 329.208 · dT 304.196; −PO₂H +H₂, sin fosfato 5′) | IDT |
| ε260 | $\sum_{i=1}^{N-1}\varepsilon_{NN}(i,i{+}1) - \sum_{i=2}^{N-1}\varepsilon_{ind}(i)$ | Cantor 1970; Cavaluzzi & Borer 2004 |
| nmol/OD₂₆₀ | $10^6 / \varepsilon_{260}$ | IDT |
| µg/OD₂₆₀ | $\text{nmol/OD} \times PM / 1000$ | IDT |

Tabla ε (L·mol⁻¹·cm⁻¹): AA 27400 · AC 21200 · AG 25000 · AT 22800 · CA 21200 · CC 14600 · CG 18000 · CT 15200 ·
GA 25200 · GC 17600 · GG 21600 · GT 20000 · TA 23400 · TC 16200 · TG 19000 · TT 16800;
individuales A 15400 · C 7400 · G 11500 · T 8700.

### 2.4 Reglas de secuencia

| Magnitud | Definición |
|---|---|
| %GC | $(G + C)/N \times 100$ |
| Pinza GC | nº de G/C en los últimos 5 nt del 3′ (la guía: 1–2, nunca > 3). **No** es el `PRIMER_GC_CLAMP` de primer3, que exige G/C consecutivas al final (rechazaría `Saureus_nuc_cF`, que termina en A). |
| Poly-X | racha más larga de la misma base |
| Estabilidad 3′ | ΔG₃₇ de los últimos 5 nt (NN, sin iniciación). Reproduce el rango del Manual 3 (−3.64 a −6.75 kcal/mol). |
| Posiciones variables | nº de columnas no conservadas dentro del primer y dentro de sus 5 nt 3′ |

## 3. Estructuras secundarias (primer3 · thal)

| Pestaña IDT | Llamada | Preset | Temperatura |
|---|---|---|---|
| Hairpin | `calc_hairpin(seq, output_structure=True)` | `fabrica` | 25 °C |
| Self-Dimer | `calc_homodimer(seq)` | `fabrica` y `pcr` | 37 °C |
| Hetero-Dimer | `calc_heterodimer(a, b)` | `pcr` (o `multiplex`) | 37 °C |
| Anclaje 3′ en dímeros | `min(calc_end_stability(a, b), calc_end_stability(b, a))` | `pcr` | 37 °C |

Todas con `mv_conc`, `dv_conc`, `dntp_conc`, `dna_conc` (nM) pasados explícitamente.
`ascii_structure_lines` da el dibujo; el backend lo traduce a `pares_bases` para el 3D.

## 4. Perfiles para imitar otras calculadoras

Un perfil = (tabla NN, convención de concentración, corrección de sal, preset por defecto, regla de Ta).
Biopython ya trae las piezas: tablas `DNA_NN1` Breslauer 1986, `DNA_NN2` Sugimoto 1996, `DNA_NN3`
Allawi & SantaLucia 1997, `DNA_NN4` SantaLucia & Hicks 2004; correcciones de sal 1–7; `Tm_Wallace`,
`Tm_GC` (fórmulas empíricas) y `chem_correction` (DMSO, formamida).

| Perfil | Tm | Sal | Defaults | Ta sugerida | Estado |
|---|---|---|---|---|---|
| **IDT OligoAnalyzer** | NN Allawi & SantaLucia 1997, $C_t/4$ | Owczarzy 2004 / 2008 | 0.25 µM · 50 mM Na⁺ · 0 Mg · 0 dNTP | — | ✅ fábrica ±0.04 °C · ⚠️ Mg ±0.74 °C |
| **Primer3 / Primer3Plus** | `primer3.calc_tm` (SantaLucia 1998) | SantaLucia 1998 u Owczarzy | 50 nM · 50 mM · 1.5 mM · 0.6 mM | — | ✅ es el mismo código |
| NEB Tm Calculator | SantaLucia 1998 + Owczarzy 2004 *(según NEB)* | buffer de cada polimerasa | según polimerasa | Q5: Tm mín + 3 · Taq: Tm mín − 5 *(confirmar)* | ⏳ calibrar |
| Thermo (Phusion/Phire) | Breslauer 1986 modificado *(según Thermo)* | buffer HF | 0.5 µM | Tm mín + 3 si > 20 nt *(confirmar)* | ⏳ calibrar |
| Benchling / SnapGene | SantaLucia 1998 (configurable) | por verificar | por verificar | — | ⏳ calibrar |
| OligoCalc | básica $64.9 + 41(G{+}C{-}16.4)/N$; ajustada $81.5 + 16.6\log_{10}[\text{Na}^+] + 0.41\,\%GC - 600/N$; NN | 50 mM Na⁺ | 50 nM | — | ⏳ calibrar |

La UI permite ver un oligo con varios perfiles a la vez ("¿por qué IDT y NEB me dan distinto?"),
lo que además sirve para la discusión de la tesis.

## 5. Validación contra tu Excel (2026-10-18)

`.venv/Scripts/python backend/test_motor.py`

| Comprobación | n | Error máx | Estado |
|---|---|---|---|
| Tm fábrica | 18 | 0.04 °C | ✅ |
| Tm PCR | 18 | 0.74 °C (medio 0.40) | ⚠️ calibrar |
| Peso molecular | 18 | 0.08 (redondeo) | ✅ |
| ε260 · ΔH · ΔS | 18 | 0 | ✅ |
| Horquilla ΔG (25 °C) | 18 | 0.97 kcal/mol (tolerancia ±1.5) | ✅ |
| Horquilla Tm | 14 | muy distinta | ⚠️ calibrar |
| Autodímero ΔG fábrica / PCR | 18 | 0 | ✅ |
| ΔG dúplex perfecto | 18 | 0.05 (redondeo) | ✅ |
| Heterodímero ΔG PCR | 73 | 0 | ✅ |
| Peor anclaje 3′ | 73 | 0 | ✅ |
| Conteos C / V / Pi / S (MEGA) | 3 archivos | 0 | ✅ |

Dato útil: la tabla de `manual_3_idt_oligoanalyzer.md` (Tm PCR 60.3, 60.2, 60.0…, horquillas −1.14, −0.42…)
se reproduce **exacta**. Lo que no cuadra son algunas columnas del Excel, que parecen venir de otra fuente
(¿la web de IDT?).

## 6. Pendientes de calibración

1. **Tm con Mg²⁺.** Si el Excel trae valores reales de IDT, nos separamos +0.40 °C de media, y no es un
   desfase constante (máx 0.74). Plan: obtener de IDT (API SciTools Plus o capturas) las Tm de 5 oligos con
   Mg = 0, 1, 2, 3 mM y ver qué término de Owczarzy 2008 implementa distinto.
2. **Horquilla (Tm y algunas ΔG).** IDT usa UNAFold (comercial). Probar en este orden: thal (actual),
   seqfold, ViennaRNA con parámetros de ADN; quedarse con el que mejor cuadre contra la API de IDT.
   Mientras tanto, la Tm de horquilla se muestra como "estimada".
3. **NEB, Thermo, Benchling, OligoCalc:** sin calibrar (§7).

### 6.1 Intento de calibración (2026-10-18, con tu tabla de 18 oligos)

- **Tm con Mg²⁺:** busqué en una rejilla de oligo 0.2–1 µM × Mg 1.5–3 mM × dNTPs 0/0.2/0.8 mM × Ct/1, /2, /4.
  Ninguna combinación reproduce tu columna (la mejor queda a ±0.43 °C). No es un problema de condiciones, sino de
  implementación de IDT. Con el preset PCR quedamos a −0.40 °C de media (máx 0.74); 6 de 18 coinciden exacto.
- **Horquillas:** comparé thal (Primer3), seqfold (Zuker, SantaLucia 2004) y ViennaRNA 2.7 con parámetros de ADN de
  Mathews 2004. Error medio en ΔG: thal 0.64, seqfold 0.71, ViennaRNA 0.82 kcal/mol. Ninguno reproduce a UNAFold,
  que además reporta estructuras con ΔG > 0 (p. ej. +1.59 kcal/mol, Tm −15.8 °C) que los motores de MFE descartan.
  Se queda thal. El verificador trata "IDT ΔG > 0 y motor sin horquilla" como equivalente.

## 7. Cómo calibrar una calculadora nueva

1. Pasar por la calculadora los 18 oligos de la tesis + 5 extremos (muy AT, muy GC, 15 y 30 nt) con sus
   valores por defecto. Guardar en `backend/calibracion/<calculadora>.csv` (`seq, condiciones, tm, ...`).
2. La prueba compara nuestro perfil contra ese CSV. Se acepta con error máx ≤ 0.5 °C en Tm.
3. Si no cuadra: ajustar **sólo** parámetros con nombre (tabla NN, convención de $C_t$, método de sal).
   Si queda un residuo sistemático, se permite un desfase explícito por perfil, documentado con su valor y n.

## 8. Diseño con primer3: mapeo de parámetros

primer3 propone; nuestro motor decide. Primer3 filtra con una ventana algo más ancha y luego cada candidato
se re-evalúa con el perfil IDT y los criterios del perfil (LOGICA_PRODUCTO.md §4–5).

| Necesidad | Parámetro primer3 | Valor (perfil Tesis) |
|---|---|---|
| Plantilla | `SEQUENCE_TEMPLATE` | consenso sin huecos |
| Región blanco | `SEQUENCE_INCLUDED_REGION` | `[inicio, largo]` elegido en el mapa |
| **Conservación** | `SEQUENCE_QUALITY` | `round(100 × conservación)` por base |
| Sin variables en el primer | `PRIMER_MIN_QUALITY` | 100 (perfil permisivo: 90) |
| 3′ 100 % conservado | `PRIMER_MIN_END_QUALITY` | 100 (aplica a los 5 nt del 3′) |
| Escala de calidad | `PRIMER_QUALITY_RANGE_MIN/MAX` | 0 / 100 |
| Largo | `PRIMER_MIN/OPT/MAX_SIZE` | 18 / 20 / 24 |
| Tm | `PRIMER_MIN/OPT/MAX_TM` | 57 / 60.5 / 64 (1 °C más ancho que el criterio) |
| Condiciones | `PRIMER_SALT_MONOVALENT`, `_DIVALENT`, `PRIMER_DNTP_CONC`, `PRIMER_DNA_CONC` | preset `pcr` (DNA en nM) |
| Fórmulas | `PRIMER_TM_FORMULA`, `PRIMER_SALT_CORRECTIONS` | 1 (SantaLucia), 2 (Owczarzy) |
| %GC | `PRIMER_MIN/MAX_GC` | 40 / 60 |
| Pinza GC | `PRIMER_GC_CLAMP` / `PRIMER_MAX_END_GC` | 0 / 3 (el mínimo de 1 lo aplica nuestro motor) |
| Repeticiones | `PRIMER_MAX_POLY_X` | 4 |
| ΔTm del par | `PRIMER_PAIR_MAX_DIFF_TM` | 2 |
| Estabilidad 3′ | `PRIMER_MAX_END_STABILITY` | 9.0 kcal/mol |
| Estructuras (prefiltro) | `PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT` + `*_TH` | 1; Tm de estructuras ≤ 47 °C (Ta − 10) |
| Tamaño de producto | `PRIMER_PRODUCT_SIZE_RANGE` | lista de rangos, sin las ventanas ocupadas del tubo |
| Cuántos | `PRIMER_NUM_RETURN` | 50 |
| Sonda qPCR | `PRIMER_PICK_INTERNAL_OLIGO` + `PRIMER_INTERNAL_*` | 1; Tm sonda = Tm primers + 8–10 °C |

Coordenadas: se guardan en columnas del alineamiento **y** en la secuencia de referencia (fila 1), que es
lo que piden `M2_primers` (MEGA) y `M6_Benchling`.

## 9. PCR virtual (algoritmo)

Una sola regla cubre pares correctos, pares invertidos, productos de un solo primer y productos cruzados
del multiplex:

1. Juntar el *pool* de primers (un par, o todo el tubo si `cruzados`).
2. Para cada secuencia `s`: sitios directos de cada primer `x` en `s` y sitios de `rc(y)` en `s`, con
   ≤ `max_desajustes` (por defecto 2, como `M4_PCR_virtual`; opción "3′ estricto": 0 desajustes en los últimos 3 nt).
3. Hay producto si un sitio de `x` en la posición *i* precede a un sitio de `rc(y)` que termina en *j*,
   con $j - i \le$ `tamano_max`. Tamaño = $j - i$.
4. Implementación: `numpy.lib.stride_tricks.sliding_window_view` + comparación vectorizada
   (18 primers × 134 secuencias × 2 kb ≈ fracciones de segundo).

Resultado esperado con la tesis: cada par 15/15 (o 14/14) en su especie, 0 en las otras 119–120 secuencias.

## 10. Conservación (algoritmo)

Por columna del alineamiento (validado contra MEGA):

- **Con hueco:** alguna fila tiene `-` → fuera de C/V (MEGA las excluye) y conservación 0 para el diseño.
- **C** (conservada): un solo estado. **V** (variable): ≥ 2 estados.
- **Pi** (parsimonia-informativa): ≥ 2 estados que aparecen cada uno en ≥ 2 secuencias. **S** = V − Pi.
- Conservación (0–1) = frecuencia de la base mayoritaria.
- **Bloques** = rachas máximas de columnas con conservación 1, sin huecos, de ≥ 18 nt.
- **Consenso estilo MultAlin** (para el esquema): por columna, frecuencia *f* de la base mayoritaria.
  *f* ≥ alto → mayúscula y rojo; *f* ≥ bajo → minúscula y azul; si no → `.` y negro.
  MultAlin usa por defecto alto = 90 % y bajo = 50 %; el preset "Tesis" usa alto = 100 %, porque con 15 cepas
  el 90 % deja pasar una cepa distinta (14/15 = 93 %).
- **Alineamiento importado de MultAlin:** si el usuario sube un `.mul` / `.msf`, se usa tal cual (no se
  re-alinea) y todas las coordenadas se refieren a él.
- Atajo: si todas las secuencias miden lo mismo, se compara por columna sin MSA (8 de tus 9 archivos;
  sólo `1_E_faecalis_ddl` tiene un indel de 24 nt).

## 11. Banco: gel, mezcla, programa

- **Gel virtual:** migración $d = a - b\log_{10}(\text{pb})$ ajustada al marcador de 100 pb (100–1000 + 1500).
- **Mezcla maestra:** $V = C_{final} \cdot V_{rx} / C_{stock}$ por reacción; total $= V \cdot n \cdot (1 + \text{exceso})$.
  Ej.: MgCl₂ 25 mM → 2 mM en 25 µL = 2 µL; × 12 × 1.1 = 26.4 µL (igual que tu Excel).
- **Ta sugerida:** Tm PCR mínima del tubo − 4 °C (configurable). Con la tesis: ~57 °C, como tu programa.
  Siempre con la recomendación de confirmar con gradiente.
