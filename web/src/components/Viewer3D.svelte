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
  let solidsGroup: THREE.Group | null = null;
  let annotGroup: THREE.Group | null = null;
  let grid: THREE.GridHelper | null = null;

  // Projected HTML labels (Ø₁, Ø₂, φ…) — DOM managed imperatively for speed.
  interface Anchor {
    text: string;
    pos: THREE.Vector3;
    /** CSS color expression (var(--…) allowed — labels are DOM). */
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

  /** Theme-dependent scene colors. */
  function palette() {
    const light = store.theme === "light";
    return {
      main: light ? 0x2f9fce : 0x88dffa,
      mainRing: light ? 0x1f7ea8 : 0x88dffa,
      branch: light ? 0xf4691f : 0xff7a3a,
      curve: light ? 0xe04c0d : 0xff5b1a,
      curveEmissive: light ? 0xb33a05 : 0xff3c00,
      axisMain: light ? 0x1f7ea8 : 0x88dffa,
      axisBranch: light ? 0xd15c17 : 0xff7a3a,
      arc: light ? 0x1b1b1d : 0xf5f5f7,
      grid1: light ? 0xc4c4be : 0x223040,
      grid2: light ? 0xdedad4 : 0x111820,
      hemiSky: light ? 0xffffff : 0xb8d8ff,
      hemiGround: light ? 0xd8d8d2 : 0x0a0a0b,
      bg: light ? "#f2f2f0" : "#0a0a0b",
    };
  }

  // ---------------------------------------------------------------------
  // Geometry helpers — everything is trimmed AT the intersection.
  // math (x, y, z) → three (x, z, −y)
  // ---------------------------------------------------------------------

  /** Point on the tilted cylinder (formula (1) of the theory). */
  function branchPoint(r2: number, phiA: number, th: number, t: number): THREE.Vector3 {
    const sp = Math.sin(phiA), cp = Math.cos(phiA);
    const x = r2 * Math.cos(th);
    const y = r2 * Math.sin(th) * cp - t * sp;
    const z = r2 * Math.sin(th) * sp + t * cp;
    return new THREE.Vector3(x, z, -y);
  }

  /** Outward normal of the tilted cylinder at angle θ (unit, three coords). */
  function branchNormal(phiA: number, th: number): THREE.Vector3 {
    const sp = Math.sin(phiA), cp = Math.cos(phiA);
    const nx = Math.cos(th);
    const ny = Math.sin(th) * cp;
    const nz = Math.sin(th) * sp;
    return new THREE.Vector3(nx, nz, -ny);
  }

  /**
   * Branch tube trimmed at its FIRST contact with the main cylinder:
   * each generator runs from the selected root t_sel(θ) to the far end,
   * on the side given by the chosen branch (outer → +t, inner → −t).
   */
  function buildTrimmedBranch(payload: IntersectionPayload): THREE.BufferGeometry | null {
    const samples = payload.dev_branch;
    if (samples.length < 3 || payload.r2 == null) return null;
    const r2 = payload.r2;
    const phiA = payload.phi;
    const outer = payload.branch !== "inner";

    let tLo = Infinity, tHi = -Infinity;
    for (const s of samples) {
      tLo = Math.min(tLo, s.v);
      tHi = Math.max(tHi, s.v);
    }
    const ext = Math.max(3 * r2, 90);
    const tEnd = outer ? tHi + ext : tLo - ext;

    const stride = Math.max(1, Math.floor(samples.length / 420));
    const cols: { th: number; t: number }[] = [];
    for (let i = 0; i < samples.length; i += stride) {
      cols.push({ th: samples[i].theta, t: samples[i].v });
    }

    const positions: number[] = [];
    const normals: number[] = [];
    const indices: number[] = [];
    for (const c of cols) {
      const a = branchPoint(r2, phiA, c.th, c.t);
      const b = branchPoint(r2, phiA, c.th, tEnd);
      const n = branchNormal(phiA, c.th);
      positions.push(a.x, a.y, a.z, b.x, b.y, b.z);
      normals.push(n.x, n.y, n.z, n.x, n.y, n.z);
    }
    const meanStep = (2 * Math.PI) / cols.length;
    const quad = (i: number, j: number) => {
      const a0 = 2 * i, a1 = 2 * i + 1, b0 = 2 * j, b1 = 2 * j + 1;
      indices.push(a0, b0, a1, b0, b1, a1);
    };
    for (let i = 0; i + 1 < cols.length; i++) {
      if (Math.abs(cols[i + 1].th - cols[i].th) < 4 * meanStep) quad(i, i + 1);
    }
    const wrapGap = cols[0].th + 2 * Math.PI - cols[cols.length - 1].th;
    if (wrapGap < 4 * meanStep) quad(cols.length - 1, 0);

    const geom = new THREE.BufferGeometry();
    geom.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    geom.setAttribute("normal", new THREE.Float32BufferAttribute(normals, 3));
    geom.setIndex(indices);
    return geom;
  }

  /**
   * Main cylinder with the gueule de loup opening carved out.
   * Per α-column, the wall is drawn outside the [zLow, zHigh] hole interval
   * computed from the developed contour (crossings of u = R1·α).
   */
  function buildHoledMain(
    payload: IntersectionPayload,
    heightMain: number,
  ): THREE.BufferGeometry {
    const r1 = payload.r1;
    const H = heightMain / 2;
    const nA = 481;

    // Hole interval per column, from the developed contour (u, v).
    let holeAt: (u0: number) => [number, number] | null = () => null;
    if (payload.dev_main && payload.dev_main_closed && payload.dev_main.length > 8) {
      const poly = payload.dev_main.map((p) => ({ u: p.u, v: p.v }));
      const circ = payload.circumference_main;
      holeAt = (u0: number) => {
        for (const shift of [0, -circ, circ]) {
          const u = u0 + shift;
          const vs: number[] = [];
          for (let i = 0; i < poly.length; i++) {
            const a = poly[i];
            const b = poly[(i + 1) % poly.length];
            if ((a.u - u) * (b.u - u) < 0) {
              const s = (u - a.u) / (b.u - a.u);
              vs.push(a.v + s * (b.v - a.v));
            }
          }
          if (vs.length >= 2) {
            return [Math.min(...vs), Math.max(...vs)];
          }
        }
        return null;
      };
    }

    interface Col {
      alpha: number;
      hole: [number, number] | null;
    }
    const cols: Col[] = [];
    for (let i = 0; i < nA; i++) {
      const alpha = (2 * Math.PI * i) / (nA - 1);
      cols.push({ alpha, hole: holeAt(r1 * alpha) ?? holeAt(r1 * (alpha - 2 * Math.PI)) });
    }

    const positions: number[] = [];
    const normals: number[] = [];
    const indices: number[] = [];
    let vi = 0;

    const push = (alpha: number, z: number) => {
      positions.push(r1 * Math.cos(alpha), z, -r1 * Math.sin(alpha));
      normals.push(Math.cos(alpha), 0, -Math.sin(alpha));
      return vi++;
    };
    const quadZ = (a: Col, b: Col, zaLo: number, zaHi: number, zbLo: number, zbHi: number) => {
      const p0 = push(a.alpha, zaLo);
      const p1 = push(a.alpha, zaHi);
      const p2 = push(b.alpha, zbLo);
      const p3 = push(b.alpha, zbHi);
      indices.push(p0, p2, p1, p2, p3, p1);
    };

    for (let i = 0; i + 1 < cols.length; i++) {
      const a = cols[i];
      const b = cols[i + 1];
      if (!a.hole && !b.hole) {
        quadZ(a, b, -H, H, -H, H);
        continue;
      }
      // Degenerate closure at the hole extremities: a hole-less column that
      // borders a holed one collapses its interval to the neighbour's middle.
      const ah: [number, number] = a.hole ?? [mid(b.hole!), mid(b.hole!)];
      const bh: [number, number] = b.hole ?? [mid(a.hole!), mid(a.hole!)];
      quadZ(a, b, -H, ah[0], -H, bh[0]); // below the opening
      quadZ(a, b, ah[1], H, bh[1], H);   // above the opening
    }

    const geom = new THREE.BufferGeometry();
    geom.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    geom.setAttribute("normal", new THREE.Float32BufferAttribute(normals, 3));
    geom.setIndex(indices);
    return geom;
  }

  function mid(h: [number, number]): number {
    return (h[0] + h[1]) / 2;
  }

  /** Tube cut by the inclined plane: kept below the cut, per α-column. */
  function buildPlaneCutMain(payload: IntersectionPayload): {
    wall: THREE.BufferGeometry;
    bottomZ: number;
  } {
    const r1 = payload.r1;
    const tan = Math.tan(payload.phi);
    const z0 = store.params.z0;
    const drop = Math.max(2.4 * r1, 140);
    const bottomZ = z0 - Math.abs(tan) * r1 - drop;
    const nA = 361;

    const positions: number[] = [];
    const normals: number[] = [];
    const indices: number[] = [];
    let vi = 0;
    const push = (alpha: number, z: number) => {
      positions.push(r1 * Math.cos(alpha), z, -r1 * Math.sin(alpha));
      normals.push(Math.cos(alpha), 0, -Math.sin(alpha));
      return vi++;
    };
    let prevLo = -1, prevHi = -1;
    for (let i = 0; i < nA; i++) {
      const alpha = (2 * Math.PI * i) / (nA - 1);
      const zTop = z0 - r1 * Math.sin(alpha) * tan;
      const lo = push(alpha, bottomZ);
      const hi = push(alpha, zTop);
      if (i > 0) indices.push(prevLo, lo, prevHi, lo, hi, prevHi);
      prevLo = lo;
      prevHi = hi;
    }
    const wall = new THREE.BufferGeometry();
    wall.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    wall.setAttribute("normal", new THREE.Float32BufferAttribute(normals, 3));
    wall.setIndex(indices);
    return { wall, bottomZ };
  }

  /** Elliptical cap filling the plane cut (triangle fan over the curve). */
  function buildPlaneCap(payload: IntersectionPayload): THREE.BufferGeometry | null {
    const pts = payload.curve3d;
    if (pts.length < 3) return null;
    const z0 = store.params.z0;
    const positions: number[] = [0, z0, 0];
    for (const p of pts) positions.push(p.x, p.z, -p.y);
    const indices: number[] = [];
    for (let i = 1; i <= pts.length; i++) {
      const j = i === pts.length ? 1 : i + 1;
      indices.push(0, i, j);
    }
    const geom = new THREE.BufferGeometry();
    geom.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    geom.setIndex(indices);
    geom.computeVertexNormals();
    return geom;
  }

  function solidMaterial(color: number): THREE.MeshStandardMaterial {
    return new THREE.MeshStandardMaterial({
      color,
      metalness: 0.3,
      roughness: 0.38,
      side: THREE.DoubleSide,
    });
  }

  function ring(r: number, tube: number, color: number): THREE.Mesh {
    const geom = new THREE.TorusGeometry(r, tube, 8, 128);
    const mat = new THREE.MeshBasicMaterial({ color, transparent: true, opacity: 0.75 });
    return new THREE.Mesh(geom, mat);
  }

  function axisLine(from: THREE.Vector3, to: THREE.Vector3, color: number): THREE.Line {
    const geom = new THREE.BufferGeometry().setFromPoints([from, to]);
    const mat = new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.5 });
    return new THREE.Line(geom, mat);
  }

  // ---------------------------------------------------------------------

  function setCamera() {
    if (!camera) return;
    const x = target.x + radius * Math.sin(phi) * Math.cos(theta);
    const z = target.z + radius * Math.sin(phi) * Math.sin(theta);
    const y = target.y + radius * Math.cos(phi);
    camera.position.set(x, y, z);
    camera.lookAt(target);
  }

  function fitToScene(payload: IntersectionPayload) {
    if (payload.curve3d.length === 0 || !solidsGroup) return;
    let cx = 0, cy = 0, cz = 0;
    for (const p of payload.curve3d) {
      cx += p.x; cy += p.z; cz += -p.y;
    }
    const n = payload.curve3d.length;
    // Aim at the joint, frame the whole assembly.
    const box = new THREE.Box3().setFromObject(solidsGroup);
    const size = box.getSize(new THREE.Vector3());
    const diag = Math.max(size.length(), payload.r1 * 4);
    target.set(cx / n, cy / n, cz / n).lerp(box.getCenter(new THREE.Vector3()), 0.45);
    radius = Math.max(300, diag * 1.15);
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

  function rebuildLabels() {
    if (!labelsEl) return;
    labelsEl.innerHTML = "";
    for (const a of anchors) {
      const el = document.createElement("div");
      el.textContent = a.text;
      el.style.cssText =
        `position:absolute;transform:translate(-50%,-130%);white-space:nowrap;` +
        `font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:0.04em;` +
        `color:${a.color};background:var(--panel-2);border:1px solid var(--line-2);` +
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
    const pal = palette();
    disposeObject(solidsGroup); solidsGroup = null;
    disposeObject(annotGroup); annotGroup = null;
    anchors = [];

    solidsGroup = new THREE.Group();
    annotGroup = new THREE.Group();

    const r1 = payload.r1;
    const phiAngle = payload.phi;
    let maxAbsZ = 0;
    for (const p of payload.curve3d) maxAbsZ = Math.max(maxAbsZ, Math.abs(p.z));

    if (payload.mode === "cyl_cyl" && payload.r2 != null) {
      const r2 = payload.r2;
      const heightMain = Math.max(4 * r1, 2.4 * maxAbsZ + 2 * r2, 200);

      // Main tube, pierced with the exact gueule de loup opening.
      const mainMesh = new THREE.Mesh(buildHoledMain(payload, heightMain), solidMaterial(pal.main));
      solidsGroup.add(mainMesh);
      const rt = Math.max(0.5, r1 * 0.012);
      for (const s of [-1, 1]) {
        const rg = ring(r1, rt, pal.mainRing);
        rg.rotation.x = Math.PI / 2;
        rg.position.y = s * (heightMain / 2);
        solidsGroup.add(rg);
      }

      // Branch tube, stopped at its first contact with the main tube.
      const branchGeom = buildTrimmedBranch(payload);
      if (branchGeom) {
        solidsGroup.add(new THREE.Mesh(branchGeom, solidMaterial(pal.branch)));
        // far-end rim of the branch
        let tLo = Infinity, tHi = -Infinity;
        for (const s of payload.dev_branch) {
          tLo = Math.min(tLo, s.v);
          tHi = Math.max(tHi, s.v);
        }
        const outer = payload.branch !== "inner";
        const ext = Math.max(3 * r2, 90);
        const tEnd = outer ? tHi + ext : tLo - ext;
        const rim = ring(r2, Math.max(0.4, r2 * 0.012), pal.branch);
        rim.rotation.x = Math.PI / 2;
        rim.position.y = tEnd;
        const holder = new THREE.Group();
        holder.add(rim);
        // Branch axis in three coords is (0, cosφ, +sinφ) — rotation +φ on X.
        holder.rotation.x = phiAngle;
        solidsGroup.add(holder);

        annotGroup.add(
          axisLine(
            new THREE.Vector3(0, 1, 0)
              .applyEuler(new THREE.Euler(phiAngle, 0, 0))
              .multiplyScalar(tEnd * 1.12),
            new THREE.Vector3(0, 0, 0),
            pal.axisBranch,
          ),
        );
        anchors.push({
          text: `Ø₂ ${(r2 * 2).toFixed(1)} mm`,
          pos: new THREE.Vector3(0, tEnd * 0.72, r2 * 1.05).applyEuler(
            new THREE.Euler(phiAngle, 0, 0),
          ),
          color: "var(--ember)",
        });
      }

      const axisLen = heightMain * 0.62;
      annotGroup.add(
        axisLine(new THREE.Vector3(0, -axisLen, 0), new THREE.Vector3(0, axisLen, 0), pal.axisMain),
      );
      anchors.push({
        text: `Ø₁ ${(r1 * 2).toFixed(1)} mm`,
        pos: new THREE.Vector3(r1 * 0.55, -heightMain * 0.4, r1 * 0.7),
        color: "var(--cyan)",
      });
    } else {
      // Cylinder × plane: the tube stops exactly at the cut.
      const { wall, bottomZ } = buildPlaneCutMain(payload);
      solidsGroup.add(new THREE.Mesh(wall, solidMaterial(pal.main)));
      const rg = ring(r1, Math.max(0.5, r1 * 0.012), pal.mainRing);
      rg.rotation.x = Math.PI / 2;
      rg.position.y = bottomZ;
      solidsGroup.add(rg);

      const cap = buildPlaneCap(payload);
      if (cap) {
        const capMat = new THREE.MeshStandardMaterial({
          color: pal.branch,
          metalness: 0.25,
          roughness: 0.4,
          side: THREE.DoubleSide,
        });
        solidsGroup.add(new THREE.Mesh(cap, capMat));
        anchors.push({
          text: `plan · φ ${((phiAngle * 180) / Math.PI).toFixed(1)}°`,
          pos: new THREE.Vector3(0, store.params.z0 + r1 * Math.abs(Math.tan(phiAngle)) * 0.5 + 14, 0),
          color: "var(--ember)",
        });
      }
      const axisLen = Math.abs(bottomZ) + r1 * 2;
      annotGroup.add(
        axisLine(new THREE.Vector3(0, bottomZ, 0), new THREE.Vector3(0, axisLen * 0.4, 0), pal.axisMain),
      );
      anchors.push({
        text: `Ø₁ ${(r1 * 2).toFixed(1)} mm`,
        pos: new THREE.Vector3(r1 * 0.55, bottomZ * 0.55, r1 * 0.7),
        color: "var(--cyan)",
      });
    }

    // Angle arc between the two axes (in the x = 0 plane).
    if (Math.abs(phiAngle) > 1e-3 && payload.mode === "cyl_cyl") {
      const rArc = Math.max(r1 * 1.9, (payload.r2 ?? 0) * 1.9, 55);
      const arcPts: THREE.Vector3[] = [];
      for (let i = 0; i <= 48; i++) {
        const s = (phiAngle * i) / 48;
        arcPts.push(new THREE.Vector3(0, Math.cos(s) * rArc, Math.sin(s) * rArc));
      }
      const arc = new THREE.Line(
        new THREE.BufferGeometry().setFromPoints(arcPts),
        new THREE.LineBasicMaterial({ color: pal.arc, transparent: true, opacity: 0.8 }),
      );
      annotGroup.add(arc);
      const midA = phiAngle / 2;
      anchors.push({
        text: `φ = ${((phiAngle * 180) / Math.PI).toFixed(1)}°`,
        pos: new THREE.Vector3(0, Math.cos(midA) * rArc * 1.18, Math.sin(midA) * rArc * 1.18),
        color: "var(--text)",
      });
    } else if (payload.mode === "cyl_plane") {
      anchors.push({
        text: `φ = ${((phiAngle * 180) / Math.PI).toFixed(1)}°`,
        pos: new THREE.Vector3(-r1 * 1.4, store.params.z0 + r1, 0),
        color: "var(--text)",
      });
    }

    // Intersection curve as an emissive tube — the star of the scene.
    if (payload.curve3d.length > 2) {
      const pts = payload.curve3d.map((p) => new THREE.Vector3(p.x, p.z, -p.y));
      const isLoop =
        payload.mode === "cyl_plane" ||
        payload.dev_main_closed ||
        pts[0].distanceTo(pts[pts.length - 1]) < r1 * 0.05;
      const curve = new THREE.CatmullRomCurve3(pts, isLoop);
      const tubeR = Math.max(0.6, Math.max(r1, payload.r2 ?? 0) * 0.012);
      const tubeGeom = new THREE.TubeGeometry(curve, Math.min(720, pts.length), tubeR, 10, isLoop);
      const pal2 = palette();
      const tubeMat = new THREE.MeshStandardMaterial({
        color: pal2.curve,
        emissive: pal2.curveEmissive,
        emissiveIntensity: 0.5,
        metalness: 0.1,
        roughness: 0.35,
      });
      solidsGroup.add(new THREE.Mesh(tubeGeom, tubeMat));
    }

    scene.add(solidsGroup);
    scene.add(annotGroup);
    rebuildLabels();
    fitToScene(payload);
  }

  function applySceneTheme() {
    if (!scene) return;
    const pal = palette();
    if (grid) {
      scene.remove(grid);
      grid.geometry.dispose();
      (grid.material as THREE.Material).dispose();
    }
    grid = new THREE.GridHelper(800, 40, pal.grid1, pal.grid2);
    (grid.material as THREE.Material).opacity = 0.3;
    (grid.material as THREE.Material).transparent = true;
    grid.position.y = -0.001;
    scene.add(grid);
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

  /** Resolve a CSS color expression (e.g. `var(--ember)`) to a concrete value. */
  function resolveCssColor(expr: string): string {
    const m = expr.match(/var\((--[a-z0-9-]+)\)/i);
    if (!m) return expr;
    return getComputedStyle(document.documentElement).getPropertyValue(m[1]).trim() || "#888";
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
      shotRenderer.setClearColor(new THREE.Color(palette().bg), 1);
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

      const light = store.theme === "light";
      const v = new THREE.Vector3();
      ctx.font = "500 24px 'JetBrains Mono', monospace";
      ctx.textBaseline = "middle";
      for (const a of anchors) {
        v.copy(a.pos).project(shotCamera);
        if (v.z >= 1 || Math.abs(v.x) > 1.1 || Math.abs(v.y) > 1.1) continue;
        const x = ((v.x + 1) / 2) * W;
        const y = ((1 - v.y) / 2) * H;
        const tw = ctx.measureText(a.text).width;
        ctx.fillStyle = light ? "rgba(255,255,255,0.85)" : "rgba(6,6,8,0.6)";
        ctx.beginPath();
        ctx.roundRect(x - tw / 2 - 12, y - 38, tw + 24, 34, 17);
        ctx.fill();
        ctx.fillStyle = resolveCssColor(a.color);
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

    const hemi = new THREE.HemisphereLight(0xffffff, 0x888888, 0.55);
    scene.add(hemi);
    const key = new THREE.DirectionalLight(0xffffff, 1.15);
    key.position.set(2, 4, 3);
    scene.add(key);
    const rim = new THREE.DirectionalLight(0xffb28a, 0.5);
    rim.position.set(-3, 2, -4);
    scene.add(rim);

    applySceneTheme();

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
      disposeObject(solidsGroup);
      disposeObject(annotGroup);
    };
  });

  $effect(() => {
    if (!scene || !store.result) return;
    void store.theme; // rebuild materials & grid when the theme flips
    applySceneTheme();
    rebuildScene(store.result);
  });
</script>

<div class="absolute inset-0">
  <div bind:this={container} class="w-full h-full"></div>
  <div bind:this={labelsEl} class="absolute inset-0 overflow-hidden pointer-events-none"></div>

  <!-- Bottom toolbar -->
  <div class="absolute left-4 bottom-4 flex items-center gap-2">
    <button
      class="pill hover:border-mist transition-colors {autoRotate ? 'text-pearl' : ''}"
      style={autoRotate ? "border-color: color-mix(in srgb, var(--ember) 45%, transparent)" : ""}
      onclick={() => (autoRotate = !autoRotate)}
      aria-pressed={autoRotate}
    >
      ⟳ rotation {autoRotate ? "on" : "off"}
    </button>
    <button
      class="pill hover:border-mist transition-colors"
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
      <div class="pill glass-strong">calcul…</div>
    </div>
  {/if}
</div>
