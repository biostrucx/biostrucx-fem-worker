// src/components/FEMViewer.jsx
import React, { useEffect, useRef } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";

// Colormap sencillo: azul → magenta → rojo
function colorMap(t) {
  const clamp = (x) => Math.max(0, Math.min(1, x));
  t = clamp(t);
  const r = clamp(1.6 * t);
  const g = clamp(0.2 * (1 - Math.abs(t - 0.5) * 2));
  const b = clamp(1.0 * (1 - t) + 0.2 * (1 - Math.abs(t - 0.5) * 2));
  return new THREE.Color(r, g, b);
}

/**
 * Props:
 *  - viz: { vertices:number[], indices:number[], u_mag?:number[], marker?:[x,y,z] }
 *  - marker?: [x,y,z]
 *  - height?: number
 *  - showMesh?: boolean  (nuevo: true por defecto)
 *  - meshOpacity?: number (0..1, default 0.18)
 */
export default function FEMViewer({
  viz,
  marker,
  height = 220,
  showMesh = true,
  meshOpacity = 0.18,
}) {
  const wrapRef = useRef(null);
  const rendererRef = useRef(null);
  const sceneRef = useRef(new THREE.Scene());
  const cameraRef = useRef(null);
  const controlsRef = useRef(null);

  const meshRef = useRef(null);
  const wireRef = useRef(null);   // NUEVO: capa wireframe
  const markerRef = useRef(null);

  const animRef = useRef(0);

  // init
  useEffect(() => {
    const wrap = wrapRef.current;
    const scene = sceneRef.current;

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.setClearAlpha(0);
    wrap.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    const cam = new THREE.PerspectiveCamera(35, 1, 0.01, 1000);
    cameraRef.current = cam;
    scene.add(cam);

    scene.add(new THREE.AmbientLight(0xffffff, 0.35));
    const hemi = new THREE.HemisphereLight(0xffffff, 0x111111, 0.9);
    scene.add(hemi);
    const dir = new THREE.DirectionalLight(0xffffff, 0.85);
    dir.position.set(2, 3, 2);
    cam.add(dir);

    const controls = new OrbitControls(cam, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controlsRef.current = controls;

    const onResize = () => {
      const w = wrap.clientWidth;
      const h = wrap.clientHeight;
      renderer.setSize(w, h, false);
      cam.aspect = Math.max(1e-6, w / h);
      cam.updateProjectionMatrix();
      if (meshRef.current) frameObject(meshRef.current);
    };
    onResize();
    window.addEventListener("resize", onResize);

    const tick = () => {
      controls.update();
      renderer.render(scene, cam);
      animRef.current = requestAnimationFrame(tick);
    };
    tick();

    return () => {
      cancelAnimationFrame(animRef.current);
      window.removeEventListener("resize", onResize);
      controls.dispose();
      renderer.dispose();
      if (renderer.domElement.parentNode) {
        renderer.domElement.parentNode.removeChild(renderer.domElement);
      }
      // limpiar objetos
      [meshRef.current, wireRef.current, markerRef.current].forEach((o) => {
        if (o) scene.remove(o);
      });
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Construir/actualizar la geometría
  useEffect(() => {
    const scene = sceneRef.current;
    if (!viz || !viz.vertices || !viz.indices) return;

    // limpiar anterior
    if (meshRef.current) scene.remove(meshRef.current);
    if (wireRef.current) scene.remove(wireRef.current);

    const verts = Array.isArray(viz.vertices[0]) ? viz.vertices.flat() : viz.vertices;
    const idx = viz.indices;

    const geom = new THREE.BufferGeometry();
    geom.setAttribute("position", new THREE.Float32BufferAttribute(verts, 3));
    geom.setIndex(idx);

    // Colores por vértice si u_mag disponible
    let material;
    if (Array.isArray(viz.u_mag) && viz.u_mag.length * 3 === verts.length) {
      const vals = viz.u_mag.map(Number).filter(Number.isFinite);
      const vmin = Math.min(...vals);
      const vmax = Math.max(...vals);
      const colors = new Float32Array(viz.u_mag.length * 3);
      for (let i = 0; i < viz.u_mag.length; i++) {
        const t = vmax - vmin > 1e-8 ? (viz.u_mag[i] - vmin) / (vmax - vmin) : 0.5;
        const c = colorMap(t);
        colors[i * 3] = c.r;
        colors[i * 3 + 1] = c.g;
        colors[i * 3 + 2] = c.b;
      }
      geom.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
      material = new THREE.MeshStandardMaterial({
        side: THREE.DoubleSide,
        flatShading: true,
        vertexColors: true,
        metalness: 0.0,
        roughness: 0.9,
      });
    } else {
      material = new THREE.MeshStandardMaterial({
        color: 0x7aa2ff,
        side: THREE.DoubleSide,
        flatShading: true,
        metalness: 0.0,
        roughness: 0.9,
      });
    }

    const mesh = new THREE.Mesh(geom, material);
    scene.add(mesh);
    meshRef.current = mesh;

    // === NUEVO: Wireframe (malla) con todas las aristas ===
    if (showMesh) {
      const wfGeom = new THREE.WireframeGeometry(geom); // incluye TODAS las aristas
      const wfMat = new THREE.LineBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: meshOpacity, // sutil
      });
      const wf = new THREE.LineSegments(wfGeom, wfMat);
      scene.add(wf);
      wireRef.current = wf;
    }

    frameObject(mesh);
    updateMarker();

    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [viz, showMesh, meshOpacity]);

  // marcador
  useEffect(() => {
    updateMarker();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [marker]);

  function frameObject(object3D) {
    const camera = cameraRef.current;
    const controls = controlsRef.current;

    const box = new THREE.Box3().setFromObject(object3D);
    const sphere = box.getBoundingSphere(new THREE.Sphere());
    const center = sphere.center;
    const r = Math.max(sphere.radius, 1e-3);

    const fov = THREE.MathUtils.degToRad(camera.fov);
    const dist = r / Math.sin(fov / 2);

    camera.near = r / 100;
    camera.far = r * 100;
    camera.updateProjectionMatrix();

    camera.position.set(center.x + dist * 0.9, center.y + dist * 0.6, center.z + dist * 0.9);
    controls.target.copy(center);
    controls.minDistance = r * 0.2;
    controls.maxDistance = r * 20;
    controls.update();
  }

  function updateMarker() {
    const scene = sceneRef.current;
    const m = marker || viz?.marker;
    if (!m || m.length < 3 || !meshRef.current) {
      if (markerRef.current) {
        scene.remove(markerRef.current);
        markerRef.current = null;
      }
      return;
    }

    const box = new THREE.Box3().setFromObject(meshRef.current);
    const r = Math.max(box.getBoundingSphere(new THREE.Sphere()).radius, 1e-3);

    const geo = new THREE.SphereGeometry(r * 0.03, 24, 24);
    const mat = new THREE.MeshStandardMaterial({
      color: 0xff4444,
      emissive: 0xaa0000,
      emissiveIntensity: 0.8,
    });
    const sph = new THREE.Mesh(geo, mat);
    sph.position.set(m[0], m[1], m[2]);

    if (markerRef.current) scene.remove(markerRef.current);
    scene.add(sph);
    markerRef.current = sph;
  }

  return (
    <div
      ref={wrapRef}
      style={{ width: "100%", height, position: "relative", borderRadius: 12, overflow: "hidden" }}
    >
      {/* Botón FIT */}
      <button
        onClick={() => meshRef.current && frameObject(meshRef.current)}
        style={{
          position: "absolute",
          right: 8,
          bottom: 8,
          fontSize: 12,
          padding: "4px 8px",
          background: "rgba(255,255,255,0.08)",
          border: "1px solid rgba(255,255,255,0.2)",
          color: "#fff",
          borderRadius: 6,
          backdropFilter: "blur(4px)",
          cursor: "pointer",
        }}
      >
        Fit
      </button>
    </div>
  );
}





