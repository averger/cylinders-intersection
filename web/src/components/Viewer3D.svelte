<script lang="ts">
  import * as THREE from "three";
  import { onMount } from "svelte";
  import { store } from "../lib/store.svelte";
  import type { IntersectionPayload } from "../lib/api";

  let container = $state<HTMLDivElement | undefined>(undefined);

  let renderer: THREE.WebGLRenderer | null = null;
  let scene: THREE.Scene | null = null;
  let camera: THREE.PerspectiveCamera | null = null;
  let raf = 0;

  // Scene-managed objects we replace on every recompute:
  let mainGroup: THREE.Group | null = null;
  let branchGroup: THREE.Group | null = null;
  let curveLine: THREE.Line | null = null;
  let planeMesh: THREE.Mesh | null = null;

  // Camera target / orbit state (lightweight orbit controls — no extra dep).
  let target = new THREE.Vector3(0, 0, 0);
  let phi = Math.PI / 3;       // polar angle from +Y down
  let theta = Math.PI / 4;     // azimuth around +Y
  let radius = 360;
  let dragging = false;
  let lastX = 0,
    lastY = 0;

  function setCamera() {
    if (!camera) return;
    const x = target.x + radius * Math.sin(phi) * Math.cos(theta);
    const z = target.z + radius * Math.sin(phi) * Math.sin(theta);
    const y = target.y + radius * Math.cos(phi);
    camera.position.set(x, y, z);
    camera.lookAt(target);
  }

  function fitToScene(payload: IntersectionPayload) {
    if (payload.curve3d.length === 0) return;
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;
    for (const p of payload.curve3d) {
      minX = Math.min(minX, p.x); maxX = Math.max(maxX, p.x);
      minY = Math.min(minY, p.y); maxY = Math.max(maxY, p.y);
      minZ = Math.min(minZ, p.z); maxZ = Math.max(maxZ, p.z);
    }
    target.set(
      (minX + maxX) / 2,
      (minZ + maxZ) / 2,   // we map z (math) -> y (three) so this is "vertical" centre
      -(minY + maxY) / 2,
    );
    const span = Math.max(maxX - minX, maxY - minY, maxZ - minZ, payload.r1 * 2);
    radius = Math.max(280, span * 2.6);
    setCamera();
  }

  function disposeGroup(group: THREE.Group | null) {
    if (!group) return;
    group.traverse((obj) => {
      const o = obj as THREE.Mesh;
      if (o.geometry) o.geometry.dispose();
      if (o.material) {
        if (Array.isArray(o.material)) o.material.forEach((m) => m.dispose());
        else o.material.dispose();
      }
    });
    group.parent?.remove(group);
  }

  function buildCylinderMesh(radius: number, height: number, color: number, opacity: number): THREE.Mesh {
    const geom = new THREE.CylinderGeometry(radius, radius, height, 96, 1, true);
    const mat = new THREE.MeshStandardMaterial({
      color,
      transparent: true,
      opacity,
      metalness: 0.25,
      roughness: 0.32,
      side: THREE.DoubleSide,
    });
    return new THREE.Mesh(geom, mat);
  }

  function rebuildScene(payload: IntersectionPayload) {
    if (!scene) return;
    disposeGroup(mainGroup); mainGroup = null;
    disposeGroup(branchGroup); branchGroup = null;
    if (curveLine) {
      curveLine.geometry.dispose();
      (curveLine.material as THREE.Material).dispose();
      scene.remove(curveLine);
      curveLine = null;
    }
    if (planeMesh) {
      planeMesh.geometry.dispose();
      (planeMesh.material as THREE.Material).dispose();
      scene.remove(planeMesh);
      planeMesh = null;
    }

    // ----- main cylinder (along world +Y in three coordinates) -----
    const r1 = payload.r1;
    const heightMain = Math.max(payload.r1 * 4, 200);
    mainGroup = new THREE.Group();
    const main = buildCylinderMesh(r1, heightMain, 0x88dffa, 0.18);
    mainGroup.add(main);

    // wireframe outline rings to keep a "drafting" feeling
    const ringGeom = new THREE.TorusGeometry(r1, 0.4, 6, 96);
    const ringMat = new THREE.MeshBasicMaterial({ color: 0x88dffa, transparent: true, opacity: 0.4 });
    for (let i = -1; i <= 1; i++) {
      const ring = new THREE.Mesh(ringGeom, ringMat);
      ring.rotation.x = Math.PI / 2;
      ring.position.y = i * (heightMain / 2 - 0.2);
      mainGroup.add(ring);
    }

    scene.add(mainGroup);

    // ----- branch cylinder OR plane -----
    if (payload.mode === "cyl_cyl" && payload.r2 != null) {
      const r2 = payload.r2;
      const heightBranch = Math.max(payload.r2 * 6, 240);
      branchGroup = new THREE.Group();
      const branch = buildCylinderMesh(r2, heightBranch, 0xff7a3a, 0.22);
      branchGroup.add(branch);

      // Tilt around X by phi: in math, +z (axis) becomes (0, sinφ, cosφ).
      // Three uses Y-up; our math z maps to three's Y. So rotation is around
      // three's X axis as well.
      branchGroup.rotation.x = -payload.phi;
      scene.add(branchGroup);
    } else if (payload.mode === "cyl_plane") {
      const size = Math.max(payload.r1 * 6, 360);
      const planeGeom = new THREE.PlaneGeometry(size, size);
      const planeMat = new THREE.MeshStandardMaterial({
        color: 0xff7a3a,
        transparent: true,
        opacity: 0.25,
        side: THREE.DoubleSide,
        metalness: 0.1,
        roughness: 0.6,
      });
      planeMesh = new THREE.Mesh(planeGeom, planeMat);
      // Plane equation in math coords: z = z0 − y·tan(phi).  In math the
      // normal is (0, sin φ, cos φ); after the (mx, my, mz) → (mx, mz, −my)
      // mapping, the three-space normal is (0, cos φ, −sin φ).  A
      // PlaneGeometry has its default normal at +Z; the rotation around
      // three.X that maps (0, 0, 1) → (0, cos φ, −sin φ) is α = −π/2 − φ.
      planeMesh.rotation.x = -Math.PI / 2 - payload.phi;
      planeMesh.position.y = store.params.z0;
      scene.add(planeMesh);
    }

    // ----- intersection curve -----
    const positions = new Float32Array(payload.curve3d.length * 3);
    for (let i = 0; i < payload.curve3d.length; i++) {
      const p = payload.curve3d[i];
      // math (x, y, z) -> three (x, z, -y)
      positions[i * 3 + 0] = p.x;
      positions[i * 3 + 1] = p.z;
      positions[i * 3 + 2] = -p.y;
    }
    const curveGeom = new THREE.BufferGeometry();
    curveGeom.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    const curveMat = new THREE.LineBasicMaterial({ color: 0xff5b1a, linewidth: 2 });
    curveLine = payload.curve3d.length > 1 ? new THREE.LineLoop(curveGeom, curveMat) : new THREE.Line(curveGeom, curveMat);
    scene.add(curveLine);

    fitToScene(payload);
  }

  function tick() {
    if (renderer && scene && camera) renderer.render(scene, camera);
    raf = requestAnimationFrame(tick);
  }

  onMount(() => {
    if (!container) return;

    scene = new THREE.Scene();
    scene.background = null;

    camera = new THREE.PerspectiveCamera(38, 1, 0.5, 5000);
    setCamera();

    renderer = new THREE.WebGLRenderer({ alpha: true, antialias: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    container.appendChild(renderer.domElement);

    // lights
    const hemi = new THREE.HemisphereLight(0xb8d8ff, 0x0a0a0b, 0.45);
    scene.add(hemi);
    const key = new THREE.DirectionalLight(0xffffff, 1.1);
    key.position.set(2, 4, 3);
    scene.add(key);
    const rim = new THREE.DirectionalLight(0xff8a50, 0.6);
    rim.position.set(-3, 2, -4);
    scene.add(rim);

    // ground grid for context
    const grid = new THREE.GridHelper(800, 40, 0x223040, 0x111820);
    (grid.material as THREE.Material).opacity = 0.25;
    (grid.material as THREE.Material).transparent = true;
    grid.position.y = -0.001;
    scene.add(grid);

    function resize() {
      if (!renderer || !camera || !container) return;
      const r = container.getBoundingClientRect();
      const w = Math.max(8, r.width);
      const h = Math.max(8, r.height);
      renderer.setSize(w, h, false);
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
    }
    const ro = new ResizeObserver(resize);
    ro.observe(container);
    resize();

    // Pointer-based orbit
    const el = renderer.domElement;
    el.style.touchAction = "none";
    el.style.cursor = "grab";

    el.addEventListener("pointerdown", (e: PointerEvent) => {
      dragging = true;
      el.setPointerCapture(e.pointerId);
      el.style.cursor = "grabbing";
      lastX = e.clientX;
      lastY = e.clientY;
    });
    el.addEventListener("pointerup", (e: PointerEvent) => {
      dragging = false;
      el.releasePointerCapture(e.pointerId);
      el.style.cursor = "grab";
    });
    el.addEventListener("pointermove", (e: PointerEvent) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      theta -= dx * 0.008;
      phi -= dy * 0.008;
      phi = Math.max(0.05, Math.min(Math.PI - 0.05, phi));
      lastX = e.clientX;
      lastY = e.clientY;
      setCamera();
    });
    el.addEventListener(
      "wheel",
      (e: WheelEvent) => {
        e.preventDefault();
        const factor = Math.exp(e.deltaY * 0.0015);
        radius = Math.max(40, Math.min(4000, radius * factor));
        setCamera();
      },
      { passive: false },
    );

    raf = requestAnimationFrame(tick);

    // Re-render whenever the result changes
    return () => {
      cancelAnimationFrame(raf);
      ro.disconnect();
      if (renderer) renderer.dispose();
      disposeGroup(mainGroup);
      disposeGroup(branchGroup);
    };
  });

  $effect(() => {
    if (!scene) return;
    if (!store.result) return;
    rebuildScene(store.result);
  });
</script>

<div class="relative w-full h-full">
  <div bind:this={container} class="w-full h-full"></div>

  <!-- HUD overlay -->
  <div class="absolute inset-x-0 top-0 p-4 flex items-start justify-between pointer-events-none">
    <div class="pill pointer-events-auto">vue 3D</div>
    <div class="text-right text-[11px] text-ash/80 num leading-tight">
      {#if store.result}
        Ø₁ {store.result.r1 * 2} mm{#if store.result.r2}
          · Ø₂ {store.result.r2 * 2} mm{/if}<br />
        φ = {((store.result.phi * 180) / Math.PI).toFixed(2)}°
      {/if}
    </div>
  </div>

  <!-- Loading shimmer -->
  {#if store.loading}
    <div class="absolute inset-0 grid place-items-center pointer-events-none">
      <div class="pill bg-black/50 backdrop-blur">calcul…</div>
    </div>
  {/if}
</div>
