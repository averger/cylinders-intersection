<script lang="ts">
  import * as THREE from "three";
  import { onMount } from "svelte";
  import { store } from "../lib/store.svelte";
  import type { IntersectionPayload } from "../lib/api";

  let container = $state<HTMLDivElement | undefined>(undefined);
  let labelsEl = $state<HTMLDivElement | undefined>(undefined);

  let renderer: THREE.WebGLRenderer | null = null;
  let scene: THREE.Scene | null = null;
  let camera: THREE.PerspectiveCamera | null = null;
  let raf = 0;

  // Scene-managed objects we replace on every recompute:
  let mainGroup: THREE.Group | null = null;
  let branchGroup: THREE.Group | null = null;
  let curveMesh: THREE.Mesh | null = null;
  let planeMesh: THREE.Mesh | null = null;
  let annotGroup: THREE.Group | null = null;

  // Projected HTML labels (Ø₁, Ø₂, φ…) — DOM managed imperatively for speed.
  interface Anchor {
    text: string;
    pos: THREE.Vector3;
    color: string;
    el?: HTMLDivElement;
  }
  let anchors: Anchor[] = [];

  // Camera target / orbit state (lightweight orbit controls — no extra dep).
  let target = new THREE.Vector3(0, 0, 0);
  let phi = Math.PI / 3;       // polar angle from +Y down
  let theta = Math.PI / 4;     // azimuth around +Y
  let radius = 360;
  let dragging = false;
  let lastX = 0,
    lastY = 0;
  let autoRotate = $state(true);
  let capturing = $state(false);

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
      (minZ + maxZ) / 2,   // math z -> three y
      -(minY + maxY) / 2,
    );
    const span = Math.max(maxX - minX, maxY - minY, maxZ - minZ, payload.r1 * 2);
    radius = Math.max(280, span * 2.6);
    setCamera();
  }

  function disposeObject(obj: THREE.Object3D | null) {
    if (!obj) return;
    obj.traverse((o) => {
      const m = o as THREE.Mesh;
      if (m.geometry) m.geometry.dispose();
      if (m.material) {
        if (Array.isArray(m.material)) m.material.forEach((mm) => mm.dispose());
        else m.material.dispose();
      }
    });
    obj.parent?.remove(obj);
  }

  function buildCylinderMesh(r: number, height: number, color: number, opacity: number): THREE.Mesh {
    const geom = new THREE.CylinderGeometry(r, r, height, 128, 1, true);
    const mat = new THREE.MeshStandardMaterial({
      color,
      transparent: true,
      opacity,
      metalness: 0.35,
      roughness: 0.28,
      side: THREE.DoubleSide,
      depthWrite: false,
    });
    return new THREE.Mesh(geom, mat);
  }

  function axisLine(from: THREE.Vector3, to: THREE.Vector3, color: number): THREE.Line {
    const geom = new THREE.BufferGeometry().setFromPoints([from, to]);
    const mat = new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.45 });
    return new THREE.Line(geom, mat);
  }

  function rebuildLabels() {
    if (!labelsEl) return;
    labelsEl.innerHTML = "";
    for (const a of anchors) {
      const el = document.createElement("div");
      el.textContent = a.text;
      el.style.cssText =
        `position:absolute;transform:translate(-50%,-130%);white-space:nowrap;` +
        `font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:0.04em;` +
        `color:${a.color};background:rgba(6,6,8,0.55);border:1px solid ${a.color}44;` +
        `padding:2px 7px;border-radius:999px;pointer-events:none;backdrop-filter:blur(6px);`;
      labelsEl.appendChild(el);
      a.el = el;
    }
  }

  function updateLabels(cam: THREE.Camera) {
    if (!renderer || !labelsEl) return;
    const w = renderer.domElement.clientWidth;
    const h = renderer.domElement.clientHeight;
    const v = new THREE.Vector3();
    for (const a of anchors) {
      if (!a.el) continue;
      v.copy(a.pos).project(cam);
      const visible = v.z < 1 && Math.abs(v.x) < 1.2 && Math.abs(v.y) < 1.2;
      a.el.style.display = visible ? "block" : "none";
      if (visible) {
        a.el.style.left = `${((v.x + 1) / 2) * w}px`;
        a.el.style.top = `${((1 - v.y) / 2) * h}px`;
      }
    }
  }

  function rebuildScene(payload: IntersectionPayload) {
    if (!scene) return;
    disposeObject(mainGroup); mainGroup = null;
    disposeObject(branchGroup); branchGroup = null;
    disposeObject(curveMesh); curveMesh = null;
    disposeObject(planeMesh); planeMesh = null;
    disposeObject(annotGroup); annotGroup = null;
    anchors = [];

    const r1 = payload.r1;
    const heightMain = Math.max(r1 * 4, 200);
    const phiAngle = payload.phi;

    // ----- main cylinder (three +Y axis) -----
    mainGroup = new THREE.Group();
    const main = buildCylinderMesh(r1, heightMain, 0x88dffa, 0.52);
    mainGroup.add(main);
    const ringGeom = new THREE.TorusGeometry(r1, Math.max(0.5, r1 * 0.012), 8, 128);
    const ringMat = new THREE.MeshBasicMaterial({ color: 0x88dffa, transparent: true, opacity: 0.6 });
    for (let i = -1; i <= 1; i += 2) {
      const ring = new THREE.Mesh(ringGeom, ringMat);
      ring.rotation.x = Math.PI / 2;
      ring.position.y = i * (heightMain / 2);
      mainGroup.add(ring);
    }
    scene.add(mainGroup);

    // ----- annotations group (axes, angle arc) -----
    annotGroup = new THREE.Group();
    const axisLen = heightMain * 0.62;
    annotGroup.add(
      axisLine(new THREE.Vector3(0, -axisLen, 0), new THREE.Vector3(0, axisLen, 0), 0x88dffa),
    );

    anchors.push({
      text: `Ø₁ ${(r1 * 2).toFixed(1)} mm`,
      pos: new THREE.Vector3(r1 * 0.75, heightMain * 0.42, r1 * 0.6),
      color: "#88dffa",
    });

    // Tilted axis direction in three coords: rotate +Y by −φ around X.
    const branchDir = new THREE.Vector3(0, 1, 0).applyEuler(new THREE.Euler(-phiAngle, 0, 0));

    if (payload.mode === "cyl_cyl" && payload.r2 != null) {
      const r2 = payload.r2;
      const heightBranch = Math.max(r2 * 6, 240);
      branchGroup = new THREE.Group();
      branchGroup.add(buildCylinderMesh(r2, heightBranch, 0xff7a3a, 0.55));
      branchGroup.rotation.x = -phiAngle;
      scene.add(branchGroup);

      annotGroup.add(
        axisLine(
          branchDir.clone().multiplyScalar(-heightBranch * 0.62),
          branchDir.clone().multiplyScalar(heightBranch * 0.62),
          0xff7a3a,
        ),
      );
      anchors.push({
        text: `Ø₂ ${(r2 * 2).toFixed(1)} mm`,
        pos: new THREE.Vector3(0, heightBranch * 0.4, r2 * 1.05).applyEuler(
          new THREE.Euler(-phiAngle, 0, 0),
        ),
        color: "#ffb28a",
      });
    } else if (payload.mode === "cyl_plane") {
      const size = Math.max(r1 * 6, 360);
      const planeGeom = new THREE.PlaneGeometry(size, size);
      const planeMat = new THREE.MeshStandardMaterial({
        color: 0xff7a3a,
        transparent: true,
        opacity: 0.4,
        side: THREE.DoubleSide,
        metalness: 0.15,
        roughness: 0.5,
        depthWrite: false,
      });
      planeMesh = new THREE.Mesh(planeGeom, planeMat);
      // three-space plane normal: (0, cosφ, −sinφ) — rotation −π/2 − φ on X.
      planeMesh.rotation.x = -Math.PI / 2 - phiAngle;
      planeMesh.position.y = store.params.z0;
      scene.add(planeMesh);
      planeMesh.updateMatrixWorld();
      anchors.push({
        text: `plan · φ ${((phiAngle * 180) / Math.PI).toFixed(1)}°`,
        pos: planeMesh.localToWorld(new THREE.Vector3(size * 0.28, size * 0.22, 0)),
        color: "#ffb28a",
      });
    }

    // ----- angle arc between the two axes (plane x = 0) -----
    if (Math.abs(phiAngle) > 1e-3) {
      const rArc = Math.max(r1 * 1.9, (payload.r2 ?? 0) * 1.9, 55);
      const arcPts: THREE.Vector3[] = [];
      const steps = 48;
      for (let i = 0; i <= steps; i++) {
        const s = (phiAngle * i) / steps;
        arcPts.push(new THREE.Vector3(0, Math.cos(s) * rArc, -Math.sin(s) * rArc));
      }
      const arcGeom = new THREE.BufferGeometry().setFromPoints(arcPts);
      const arc = new THREE.Line(
        arcGeom,
        new THREE.LineBasicMaterial({ color: 0xf5f5f7, transparent: true, opacity: 0.8 }),
      );
      annotGroup.add(arc);
      const mid = phiAngle / 2;
      anchors.push({
        text: `φ = ${((phiAngle * 180) / Math.PI).toFixed(1)}°`,
        pos: new THREE.Vector3(0, Math.cos(mid) * rArc * 1.18, -Math.sin(mid) * rArc * 1.18),
        color: "#f5f5f7",
      });
    }
    scene.add(annotGroup);

    // ----- intersection curve as a fat tube -----
    if (payload.curve3d.length > 2) {
      const pts = payload.curve3d.map((p) => new THREE.Vector3(p.x, p.z, -p.y));
      const isLoop = payload.mode === "cyl_plane" || (payload.dev_main_closed ?? false) ||
        pts[0].distanceTo(pts[pts.length - 1]) < r1 * 0.05;
      const curve = new THREE.CatmullRomCurve3(pts, isLoop);
      const tubeR = Math.max(0.6, Math.max(r1, payload.r2 ?? 0) * 0.014);
      const tubeGeom = new THREE.TubeGeometry(curve, Math.min(720, pts.length), tubeR, 10, isLoop);
      const tubeMat = new THREE.MeshStandardMaterial({
        color: 0xff5b1a,
        emissive: 0xff3c00,
        emissiveIntensity: 0.55,
        metalness: 0.1,
        roughness: 0.35,
      });
      curveMesh = new THREE.Mesh(tubeGeom, tubeMat);
      scene.add(curveMesh);
    }

    rebuildLabels();
    fitToScene(payload);
  }

  function tick() {
    if (autoRotate && !dragging) {
      theta += 0.0022;
      setCamera();
    }
    if (renderer && scene && camera) {
      renderer.render(scene, camera);
      updateLabels(camera);
    }
    raf = requestAnimationFrame(tick);
  }

  /// Render one Full HD frame off-screen (labels composited on a 2D canvas)
  /// and download it as PNG — 1920×1080, ready for social media.
  function capturePNG() {
    if (!scene || !camera || capturing) return;
    capturing = true;
    try {
      const W = 1920, H = 1080;
      const shotRenderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
      shotRenderer.setPixelRatio(1);
      shotRenderer.setSize(W, H, false);
      shotRenderer.outputColorSpace = THREE.SRGBColorSpace;
      shotRenderer.setClearColor(new THREE.Color("#0a0a0b"), 1);
      const shotCamera = camera.clone();
      shotCamera.aspect = W / H;
      shotCamera.updateProjectionMatrix();
      shotRenderer.render(scene, shotCamera);

      const out = document.createElement("canvas");
      out.width = W;
      out.height = H;
      const ctx = out.getContext("2d")!;
      ctx.drawImage(shotRenderer.domElement, 0, 0);
      shotRenderer.dispose();

      // Composite the annotation labels at their projected positions.
      const v = new THREE.Vector3();
      ctx.font = "500 24px 'JetBrains Mono', monospace";
      ctx.textBaseline = "middle";
      for (const a of anchors) {
        v.copy(a.pos).project(shotCamera);
        if (v.z >= 1 || Math.abs(v.x) > 1.1 || Math.abs(v.y) > 1.1) continue;
        const x = ((v.x + 1) / 2) * W;
        const y = ((1 - v.y) / 2) * H;
        const tw = ctx.measureText(a.text).width;
        ctx.fillStyle = "rgba(6,6,8,0.6)";
        ctx.beginPath();
        ctx.roundRect(x - tw / 2 - 12, y - 38, tw + 24, 34, 17);
        ctx.fill();
        ctx.fillStyle = a.color;
        ctx.fillText(a.text, x - tw / 2, y - 21);
      }

      const a = document.createElement("a");
      a.href = out.toDataURL("image/png");
      a.download = "cylix-3d-1920x1080.png";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } finally {
      capturing = false;
    }
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
    const hemi = new THREE.HemisphereLight(0xb8d8ff, 0x0a0a0b, 0.5);
    scene.add(hemi);
    const key = new THREE.DirectionalLight(0xffffff, 1.2);
    key.position.set(2, 4, 3);
    scene.add(key);
    const rim = new THREE.DirectionalLight(0xff8a50, 0.7);
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
      autoRotate = false;
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
      theta -= (e.clientX - lastX) * 0.008;
      phi -= (e.clientY - lastY) * 0.008;
      phi = Math.max(0.05, Math.min(Math.PI - 0.05, phi));
      lastX = e.clientX;
      lastY = e.clientY;
      setCamera();
    });
    el.addEventListener(
      "wheel",
      (e: WheelEvent) => {
        e.preventDefault();
        radius = Math.max(40, Math.min(4000, radius * Math.exp(e.deltaY * 0.0015)));
        setCamera();
      },
      { passive: false },
    );

    raf = requestAnimationFrame(tick);

    return () => {
      cancelAnimationFrame(raf);
      ro.disconnect();
      if (renderer) renderer.dispose();
      disposeObject(mainGroup);
      disposeObject(branchGroup);
      disposeObject(curveMesh);
      disposeObject(planeMesh);
      disposeObject(annotGroup);
    };
  });

  $effect(() => {
    if (!scene) return;
    if (!store.result) return;
    rebuildScene(store.result);
  });
</script>

<div class="absolute inset-0">
  <div bind:this={container} class="w-full h-full"></div>
  <div bind:this={labelsEl} class="absolute inset-0 overflow-hidden pointer-events-none"></div>

  <!-- Bottom toolbar -->
  <div class="absolute left-4 bottom-4 flex items-center gap-2">
    <button
      class="pill hover:border-white/25 transition-colors {autoRotate ? 'text-pearl' : ''}"
      style={autoRotate ? "border-color: rgba(255,91,26,0.45)" : ""}
      onclick={() => (autoRotate = !autoRotate)}
      aria-pressed={autoRotate}
    >
      ⟳ rotation {autoRotate ? "on" : "off"}
    </button>
    <button
      class="pill hover:border-white/25 transition-colors"
      onclick={capturePNG}
      disabled={capturing}
      aria-label="Capturer la scène en PNG 1920×1080"
    >
      ◉ png · 1920×1080
    </button>
  </div>

  <!-- Loading shimmer -->
  {#if store.loading}
    <div class="absolute inset-0 grid place-items-center pointer-events-none">
      <div class="pill bg-black/50 backdrop-blur">calcul…</div>
    </div>
  {/if}
</div>
