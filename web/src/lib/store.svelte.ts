/**
 * Cylix application state — multi-study (études en onglets, à la pilegroupx),
 * persisted in localStorage.  Each study owns its geometry parameters and its
 * export-editor state; the compute result is held for the active study only
 * and recomputed on switch (the backend answers in milliseconds).
 */

import { api, type Branch, type IntersectionPayload, type MultiPayload } from "./api";
import type {
  Annotation,
  Layers,
  PageSpec,
  ScaleMode,
  TitleBlock,
} from "./export";

export type Mode = "cyl_cyl" | "cyl_plane" | "multi";
export type View = "3d" | "2d";
export type Theme = "light" | "dark";

/** One branch of a multi-branch node, in UI units (mm / degrees).
 * `z = 0` (the default) puts the axis through the CENTRE of the main tube:
 * concurrent axes, like in the two-cylinder mode.  A non-zero `z` slides the
 * axis along the chord and creates the ECCENTRICITY of a real K / N joint. */
export interface BranchParam {
  d: number;          // branch diameter, mm
  z: number;          // axis crossing point along the chord axis, mm
  angleDeg: number;   // inclination from the main axis, (0°, 180°)
  azimutDeg: number;  // azimuth around the main tube, degrees
}

export interface Params {
  mode: Mode;
  d1: number;       // diameter of main cyl, mm
  d2: number;       // diameter of branch (cyl-cyl only), mm
  angleDeg: number; // angle between axes (cyl-cyl) or plane X-tilt (cyl-plane)
  angleYDeg: number; // plane Y-tilt (oriented plane, cyl-plane only)
  branch: Branch;
  z0: number;       // plane offset (cyl-plane), mm
  samples: number;
  /** Branches of the multi mode ("châssis"). */
  branches: BranchParam[];
}

export interface EditorState {
  page: PageSpec;
  scale: ScaleMode;
  layers: Layers;
  titleBlock: TitleBlock;
  annotations: Annotation[];
  cutWidth: number;
}

export interface Study {
  id: string;
  name: string;
  params: Params;
  editor: EditorState;
}

/** Default node: three tubes fanning onto the main tube, all axes through
 * its centre, no azimuth — a classic truss node.  The list order gives the
 * priority: 2 dies on 1, 3 dies on 1 and 2. */
export function defaultBranches(): BranchParam[] {
  return [
    { d: 60, z: 0, angleDeg: 45, azimutDeg: 0 },
    { d: 50, z: 0, angleDeg: 90, azimutDeg: 0 },
    { d: 45, z: 0, angleDeg: 135, azimutDeg: 0 },
  ];
}

const DEFAULT_PARAMS: Params = {
  mode: "cyl_cyl",
  d1: 100,
  d2: 70,
  angleDeg: 45,
  angleYDeg: 0,
  branch: "outer",
  z0: 0,
  samples: 1440,
  branches: defaultBranches(),
};

function today(): string {
  const d = new Date();
  const mm = String(d.getMonth() + 1).padStart(2, "0");
  const dd = String(d.getDate()).padStart(2, "0");
  return `${d.getFullYear()}-${mm}-${dd}`;
}

function defaultEditor(): EditorState {
  return {
    page: { format: "a4", orientation: "landscape", margin_mm: 10 },
    scale: "one_to_one",
    layers: { grid: true, frame: true, axis: true, labels: true, scale_bar: true },
    titleBlock: { title: "", project: "", author: "", date: today(), notes: "" },
    annotations: [],
    cutWidth: 0.35,
  };
}

let seq = 0;
function newId(): string {
  seq += 1;
  return `s${Date.now().toString(36)}${seq}`;
}

function makeStudy(name: string): Study {
  return { id: newId(), name, params: { ...DEFAULT_PARAMS }, editor: defaultEditor() };
}

const STORAGE_KEY = "cylix.studies.v1";

/** Visibility of the 3D scene layers — what the eye needs right now. */
export interface Show3D {
  main: boolean;     // main tube (wall, end rings)
  cutters: boolean;  // inclined tube(s) / plane
  curves: boolean;   // cut rims and opening contours
  angles: boolean;   // angle protractors (φ, ψ) with their chips
  labels: boolean;   // Ø chips with leader lines
  axes: boolean;     // axis lines and diameter lines
  grid: boolean;     // ground grid
}

export function defaultShow3D(): Show3D {
  return {
    main: true,
    cutters: true,
    curves: true,
    angles: true,
    labels: true,
    axes: true,
    grid: true,
  };
}

interface PersistShape {
  studies: Study[];
  activeId: string;
  view: View;
  theme?: Theme;
  show3d?: Show3D;
}

function load(): PersistShape | null {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return null;
    const data = JSON.parse(raw) as PersistShape;
    if (!Array.isArray(data.studies) || data.studies.length === 0) return null;
    // Merge with defaults so older payloads gain new fields gracefully.
    data.studies = data.studies.map((s) => {
      const params = { ...DEFAULT_PARAMS, ...s.params };
      // Le domaine gueule de loup impose Ø₂ ≤ Ø₁ — pour chaque piquage aussi.
      params.d2 = Math.min(params.d2, params.d1);
      if (!Array.isArray(params.branches) || params.branches.length === 0) {
        params.branches = defaultBranches();
      }
      params.branches = params.branches.map((b) => ({
        ...b,
        d: Math.min(b.d, params.d1),
        z: b.z ?? 0, // studies saved before the eccentricity slider
      }));
      return {
        id: s.id ?? newId(),
        name: s.name ?? "Étude",
        params,
        editor: { ...defaultEditor(), ...s.editor },
      };
    });
    return data;
  } catch {
    return null;
  }
}

class AppStore {
  studies = $state<Study[]>([makeStudy("Étude 1")]);
  activeId = $state<string>("");
  view = $state<View>("3d");
  /** Light is the default — dark is the opt-in "mission control" mode. */
  theme = $state<Theme>("light");
  /** 3D layer visibility (persisted, shared by all studies). */
  show3d = $state<Show3D>(defaultShow3D());

  result = $state<IntersectionPayload | null>(null);
  multiResult = $state<MultiPayload | null>(null);
  loading = $state(false);
  error = $state<string | null>(null);
  lastComputedAt = $state<number>(0);

  /** Monotonic token — drops stale responses while scrubbing sliders. */
  private token = 0;
  private persistTimer: ReturnType<typeof setTimeout> | null = null;

  constructor() {
    const saved = typeof localStorage !== "undefined" ? load() : null;
    if (saved) {
      this.studies = saved.studies;
      this.activeId = saved.studies.some((s) => s.id === saved.activeId)
        ? saved.activeId
        : saved.studies[0].id;
      this.view = saved.view === "2d" ? "2d" : "3d";
      this.theme = saved.theme === "dark" ? "dark" : "light";
      if (saved.show3d) this.show3d = { ...defaultShow3D(), ...saved.show3d };
    } else {
      this.activeId = this.studies[0].id;
    }
  }

  get active(): Study {
    return this.studies.find((s) => s.id === this.activeId) ?? this.studies[0];
  }

  get params(): Params {
    return this.active.params;
  }

  set params(p: Params) {
    this.active.params = p;
    this.persist();
  }

  get editor(): EditorState {
    return this.active.editor;
  }

  // ----- studies ----------------------------------------------------------
  addStudy() {
    const n = this.studies.length + 1;
    const study = makeStudy(`Étude ${n}`);
    // The new study starts from the current one — the usual workflow is
    // "same piquage, variant angle".
    study.params = { ...this.active.params };
    this.studies.push(study);
    this.switchStudy(study.id);
  }

  switchStudy(id: string) {
    if (id === this.activeId) return;
    if (!this.studies.some((s) => s.id === id)) return;
    this.activeId = id;
    this.persist();
    this.compute();
  }

  closeStudy(id: string) {
    if (this.studies.length <= 1) return;
    const idx = this.studies.findIndex((s) => s.id === id);
    if (idx < 0) return;
    this.studies.splice(idx, 1);
    if (this.activeId === id) {
      this.activeId = this.studies[Math.max(0, idx - 1)].id;
      this.compute();
    }
    this.persist();
  }

  renameStudy(id: string, name: string) {
    const s = this.studies.find((x) => x.id === id);
    if (s) {
      s.name = name.trim() || s.name;
      this.persist();
    }
  }

  setView(view: View) {
    this.view = view;
    this.persist();
  }

  toggleTheme() {
    this.theme = this.theme === "light" ? "dark" : "light";
    this.persist();
  }

  toggle3d(key: keyof Show3D) {
    this.show3d = { ...this.show3d, [key]: !this.show3d[key] };
    this.persist();
  }

  persist() {
    if (typeof localStorage === "undefined") return;
    if (this.persistTimer) clearTimeout(this.persistTimer);
    this.persistTimer = setTimeout(() => {
      const payload: PersistShape = {
        studies: JSON.parse(JSON.stringify(this.studies)),
        activeId: this.activeId,
        view: this.view,
        theme: this.theme,
        show3d: { ...this.show3d },
      };
      try {
        localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
      } catch {
        /* quota — non-fatal */
      }
    }, 250);
  }

  /** Whether the active mode has a usable computation result. */
  get hasResult(): boolean {
    return this.params.mode === "multi" ? this.multiResult !== null : this.result !== null;
  }

  // ----- compute ----------------------------------------------------------
  async compute() {
    const myToken = ++this.token;
    this.loading = true;
    this.error = null;
    try {
      const p = this.params;
      const phi = (p.angleDeg * Math.PI) / 180;
      if (p.mode === "multi") {
        const res = await api.multi({
          r1: p.d1 / 2,
          branches: p.branches.map((b) => ({
            r: b.d / 2,
            z: b.z ?? 0,
            phi: (b.angleDeg * Math.PI) / 180,
            psi: (b.azimutDeg * Math.PI) / 180,
          })),
          n_samples: Math.min(p.samples, 720),
        });
        if (myToken !== this.token) return;
        this.multiResult = res;
        this.result = null;
        this.lastComputedAt = Date.now();
        return;
      }
      let res: IntersectionPayload;
      if (p.mode === "cyl_cyl") {
        res = await api.cylCyl({
          r1: p.d1 / 2,
          r2: p.d2 / 2,
          phi,
          n_samples: p.samples,
          branch: p.branch,
        });
      } else {
        res = await api.cylPlane({
          r1: p.d1 / 2,
          phi,
          phi_y: (p.angleYDeg * Math.PI) / 180,
          z0: p.z0,
          n_samples: p.samples,
        });
      }
      if (myToken !== this.token) return;
      this.result = res;
      this.multiResult = null;
      this.lastComputedAt = Date.now();
    } catch (e: unknown) {
      if (myToken !== this.token) return;
      this.error = e instanceof Error ? e.message : String(e);
      this.result = null;
      this.multiResult = null;
    } finally {
      if (myToken === this.token) this.loading = false;
    }
  }
}

export const store = new AppStore();
