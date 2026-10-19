// Hélice B idealizada generada por fórmula (docs/ARQUITECTURA.md §3). Ilustrativa, no es una predicción estructural.
import { Html, OrbitControls, PerspectiveCamera } from "@react-three/drei";
import { Canvas } from "@react-three/fiber";
import { useLayoutEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { Escena, Estructura, Par, complemento, rango } from "./lib";

const RISE = 3.38, TWIST = (34.3 * Math.PI) / 180, R = 10, FASE = (140 * Math.PI) / 180;
const BASE: Record<string, string> = { A: "#16a34a", T: "#dc2626", G: "#ca8a04", C: "#2563eb" };
const VERDE = "#22c55e", ROJO = "#ef4444", GRIS = "#94a3b8", AZUL = "#3b82f6", NARANJA = "#f97316";

const punto = (k: number, h: number, r = R) => {
  const a = k * TWIST + h * FASE;
  return new THREE.Vector3(k * RISE, r * Math.cos(a), r * Math.sin(a));
};

export function escenaAmplicon(par: Par): Escena {
  const s = par.producto.seq, n = s.length;
  return {
    titulo: `Amplicón ${par.F.nombre} → ${par.R.nombre} · ${n} pb`,
    hebras: [
      { seq: s, k: rango(n), colores: rango(n).map((i) => (i < par.F.largo ? VERDE : GRIS)) },
      { seq: complemento(s), k: rango(n), colores: rango(n).map((i) => (i >= n - par.R.largo ? ROJO : "#cbd5e1")) },
    ],
    pares: rango(n).map((i) => [i, i]),
    focos: [{ nombre: "Forward", k: par.F.largo / 2 }, { nombre: "Reverse", k: n - par.R.largo / 2 }],
  };
}

export function escenaEstructura(e: Estructura, titulo: string): Escena {
  return {
    titulo, lazo: e.lazo, pares: e.pares,
    hebras: e.hebras.map((h, i) => ({ ...h, colores: h.seq.split("").map(() => (i ? NARANJA : AZUL)) })),
  };
}

type Item = { m: THREE.Matrix4; color: string };

function Instancias({ items, children }: { items: Item[]; children: React.ReactNode }) {
  const ref = useRef<THREE.InstancedMesh>(null!);
  useLayoutEffect(() => {
    const c = new THREE.Color();
    items.forEach((it, i) => { ref.current.setMatrixAt(i, it.m); ref.current.setColorAt(i, c.set(it.color)); });
    ref.current.instanceMatrix.needsUpdate = true;
    if (ref.current.instanceColor) ref.current.instanceColor.needsUpdate = true;
  }, [items]);
  return (
    <instancedMesh key={items.length} ref={ref} args={[undefined, undefined, items.length]}>
      {children}
      <meshStandardMaterial />
    </instancedMesh>
  );
}

function Helice({ escena }: { escena: Escena }) {
  const { esferas, palos, curvas, etiquetas } = useMemo(() => {
    const esferas: Item[] = [], palos: Item[] = [], curvas: { pts: THREE.Vector3[]; color: string }[] = [];
    const etiquetas: { p: THREE.Vector3; t: string }[] = [];
    const apareada = [new Set(escena.pares.map((x) => x[0])), new Set(escena.pares.map((x) => x[1]))];
    const arriba = new THREE.Vector3(0, 1, 0);
    escena.hebras.forEach((h, hi) => {
      const pts = h.k.map((k) => punto(k, hi));
      curvas.push({ pts, color: hi ? "#64748b" : "#475569" });
      h.seq.split("").forEach((b, i) => {
        const p = pts[i], k = h.k[i];
        esferas.push({ m: new THREE.Matrix4().setPosition(p), color: h.colores?.[i] ?? GRIS });
        const hasta = apareada[hi].has(i) ? p.clone().lerp(punto(k, 1 - hi), 0.5) : p.clone().lerp(new THREE.Vector3(k * RISE, 0, 0), 0.5);
        const dir = hasta.clone().sub(p);
        const m = new THREE.Matrix4().compose(
          p.clone().add(hasta).multiplyScalar(0.5),
          new THREE.Quaternion().setFromUnitVectors(arriba, dir.clone().normalize()),
          new THREE.Vector3(1, dir.length(), 1),
        );
        palos.push({ m, color: BASE[b] ?? GRIS });
      });
      // Arriba: 5'→3'. Abajo: escrita 3'→5'.
      etiquetas.push({ p: punto(h.k[0] - 1.5, hi, R + 3), t: hi ? "3′" : "5′" });
      etiquetas.push({ p: punto(h.k[h.k.length - 1] + 1.5, hi, R + 3), t: hi ? "5′" : "3′" });
    });
    if (escena.lazo && escena.hebras.length === 2) {
      const a = curvas[0].pts.at(-1)!, b = curvas[1].pts.at(-1)!;
      const medio = a.clone().add(b).multiplyScalar(0.5).add(new THREE.Vector3(RISE * 4, 0, 0));
      curvas.push({ pts: new THREE.QuadraticBezierCurve3(a, medio, b).getPoints(16), color: "#475569" });
    }
    return { esferas, palos, curvas, etiquetas };
  }, [escena]);

  return (
    <>
      <Instancias items={esferas}><sphereGeometry args={[1.7, 14, 14]} /></Instancias>
      <Instancias items={palos}><cylinderGeometry args={[0.7, 0.7, 1, 8]} /></Instancias>
      {curvas.filter((c) => c.pts.length > 1).map((c, i) => (
        <mesh key={i}>
          <tubeGeometry args={[new THREE.CatmullRomCurve3(c.pts), c.pts.length * 4, 0.6, 6, false]} />
          <meshStandardMaterial color={c.color} />
        </mesh>
      ))}
      {etiquetas.map((e, i) => (
        <Html key={i} position={e.p} center><span className="etiqueta3d">{e.t}</span></Html>
      ))}
    </>
  );
}

export function Vista3D({ escena }: { escena: Escena | null }) {
  const [foco, setFoco] = useState<number | null>(null);
  if (!escena) return <p className="nota">Abre una vista 3D desde el diseño (amplicón, heterodímero) o desde el analizador (horquilla, dímeros).</p>;
  const ks = escena.hebras.flatMap((h) => h.k);
  const kmin = Math.min(...ks), kmax = Math.max(...ks);
  const cx = (foco ?? (kmin + kmax) / 2) * RISE;
  const dist = foco != null ? 110 : Math.max(90, (kmax - kmin) * RISE * 0.75);
  return (
    <section>
      <h2>{escena.titulo}</h2>
      <div className="fila">
        {escena.focos && <button onClick={() => setFoco(null)}>Ver todo</button>}
        {escena.focos?.map((f) => <button key={f.nombre} onClick={() => setFoco(f.k)}>Ir al {f.nombre}</button>)}
        <span className="nota">
          Arrastra para girar, rueda para acercar. Bases: <b style={{ color: BASE.A }}>A</b> <b style={{ color: BASE.T }}>T</b>{" "}
          <b style={{ color: BASE.G }}>G</b> <b style={{ color: BASE.C }}>C</b>.{" "}
          {escena.focos ? "Verde = Forward, rojo = Reverse." : "Azul = oligo de arriba (5′→3′), naranja = el de abajo (3′→5′); peldaño completo = par de bases."}
          {" "}Geometría B idealizada: ilustrativa.
        </span>
      </div>
      <div className="lienzo">
        <Canvas>
          <PerspectiveCamera key={`c${foco}${escena.titulo}`} makeDefault position={[cx, 0, dist]} fov={45} near={1} far={20000} />
          <OrbitControls key={`o${foco}${escena.titulo}`} target={[cx, 0, 0]} />
          <ambientLight intensity={0.7} />
          <directionalLight position={[50, 80, 120]} intensity={1.2} />
          <Helice escena={escena} />
        </Canvas>
      </div>
    </section>
  );
}
