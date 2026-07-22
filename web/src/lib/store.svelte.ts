/**
 * Cylix application state — multi-study (études en onglets, à la pilegroupx),
 * persisted in localStorage.  Each study owns its geometry parameters and its
 * export-editor state; the compute result is held for the active study only
 * and recomputed on switch (the backend answers in milliseconds).
 */

import { api, type Branch, type IntersectionPayload } from "./api";
import type {
  Annotation,
  Layers,
  PageSpec,
  ScaleMode,
  TitleBlock,
} from "./export";

export type Mode = "cyl_cyl" | "cyl_plane";
export type View = "3d" | "2d";

export interface Params {
  mode: Mode;
  d1: number;       // diameter of main cyl, mm
  d2: number;       // diameter of branch (cyl-cyl only), mm
  angleDeg: number; // angle between axes (cyl-cyl) or plane tilt (cyl-plane)
  branch: Branch;
  z0: number;       // plane offset (cyl-plane), mm
  samples: number;
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

const DEFAULT_PARAMS: Params = {
  mode: "cyl_cyl",
  d1: 100,
  d2: 70,
  angleDeg: 45,
  branch: "outer",
  z0: 0,
  samples: 1440,
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

interface PersistShape {
  studies: Study[];
  activeId: string;
  view: View;
}

function load(): PersistShape | null {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return null;
    const data = JSON.parse(raw) as PersistShape;
    if (!Array.isArray(data.studies) || data.studies.length === 0) return null;
    // Merge with defaults so older payloads gain new fields gracefully.
    data.studies = data.studies.map((s) => ({
      id: s.id ?? newId(),
      name: s.name ?? "Étude",
      params: { ...DEFAULT_PARAMS, ...s.params },
      editor: { ...defaultEditor(), ...s.editor },
    }));
    return data;
  } catch {
    return null;
  }
}

class AppStore {
  studies = $state<Study[]>([makeStudy("Étude 1")]);
  activeId = $state<string>("");
  view = $state<View>("3d");

  result = $state<IntersectionPayload | null>(null);
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

  persist() {
    if (typeof localStorage === "undefined") return;
    if (this.persistTimer) clearTimeout(this.persistTimer);
    this.persistTimer = setTimeout(() => {
      const payload: PersistShape = {
        studies: JSON.parse(JSON.stringify(this.studies)),
        activeId: this.activeId,
        view: this.view,
      };
      try {
        localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
      } catch {
        /* quota — non-fatal */
      }
    }, 250);
  }

  // ----- compute ----------------------------------------------------------
  async compute() {
    const myToken = ++this.token;
    this.loading = true;
    this.error = null;
    try {
      const p = this.params;
      const phi = (p.angleDeg * Math.PI) / 180;
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
        res = await api.cylPlane({ r1: p.d1 / 2, phi, z0: p.z0, n_samples: p.samples });
      }
      if (myToken !== this.token) return;
      this.result = res;
      this.lastComputedAt = Date.now();
    } catch (e: unknown) {
      if (myToken !== this.token) return;
      this.error = e instanceof Error ? e.message : String(e);
      this.result = null;
    } finally {
      if (myToken === this.token) this.loading = false;
    }
  }
}

export const store = new AppStore();
