/**
 * Export-studio facade — reads and writes the editor state of the ACTIVE
 * study (each study keeps its own page setup, layers and annotations), and
 * drives the backend PDF / DXF exports.
 */

import { store } from "./store.svelte";
import {
  exportApi,
  planTiles,
  tileLabel,
  type ExportDocument,
  type Layers,
  type PageSpec,
  type PatternKind,
  type ScaleMode,
  type SourceSpec,
  type TitleBlock,
} from "./export";

class EditorStore {
  /** UI-only state (not persisted per study). */
  selected = $state<number | null>(null);
  busy = $state<null | "pdf" | "dxf">(null);
  error = $state<string | null>(null);

  // ----- persisted, per active study --------------------------------------
  get page(): PageSpec {
    return store.editor.page;
  }
  set page(v: PageSpec) {
    store.editor.page = v;
    store.persist();
  }

  get scale(): ScaleMode {
    return store.editor.scale;
  }
  set scale(v: ScaleMode) {
    store.editor.scale = v;
    store.persist();
  }

  get layers(): Layers {
    return store.editor.layers;
  }
  set layers(v: Layers) {
    store.editor.layers = v;
    store.persist();
  }

  get titleBlock(): TitleBlock {
    return store.editor.titleBlock;
  }
  set titleBlock(v: TitleBlock) {
    store.editor.titleBlock = v;
    store.persist();
  }

  get annotations() {
    return store.editor.annotations;
  }

  get cutWidth(): number {
    return store.editor.cutWidth;
  }
  set cutWidth(v: number) {
    store.editor.cutWidth = v;
    store.persist();
  }

  /** Patterns available for the current computation result. */
  get available(): PatternKind[] {
    return store.result?.dev_main ? ["branch", "main"] : ["branch"];
  }

  /** Rebuild the backend source spec from the compute parameters. */
  source(): SourceSpec {
    const p = store.params;
    const phi = (p.angleDeg * Math.PI) / 180;
    return p.mode === "cyl_cyl"
      ? {
          mode: "cyl_cyl",
          r1: p.d1 / 2,
          r2: p.d2 / 2,
          phi,
          n_samples: p.samples,
          branch: p.branch,
        }
      : { mode: "cyl_plane", r1: p.d1 / 2, phi, z0: p.z0, n_samples: p.samples };
  }

  document(): ExportDocument {
    const e = store.editor;
    return {
      source: this.source(),
      page: { ...e.page },
      scale: e.scale,
      layers: { ...e.layers },
      title_block: { ...e.titleBlock },
      annotations: e.annotations.map((a) => ({ ...a })),
      patterns: this.available,
      cut_width_mm: e.cutWidth,
    };
  }

  /** Total exported page count, matching the backend layout logic. */
  pageCount(): number {
    if (!store.result) return 0;
    if (this.scale === "fit") return this.available.length;
    let pages = 0;
    for (const kind of this.available) pages += this.tiles(kind).length;
    return pages;
  }

  /** Points of a pattern's cut curve. */
  private cutPoints(kind: PatternKind): { u: number; v: number }[] {
    const r = store.result;
    if (!r) return [];
    const pts = kind === "branch" ? r.dev_branch : r.dev_main;
    return pts ?? [];
  }

  /**
   * Bounding box of what must reach the paper: cut curve + annotations.
   * Twin of the backend `Sheet::cut_bbox`.
   */
  cutBox(kind: PatternKind): { uMin: number; vMin: number; uMax: number; vMax: number } | null {
    const pts = this.cutPoints(kind);
    if (pts.length === 0) return null;
    let uMin = pts[0].u,
      uMax = uMin,
      vMin = pts[0].v,
      vMax = vMin;
    for (const p of pts) {
      if (p.u < uMin) uMin = p.u;
      if (p.u > uMax) uMax = p.u;
      if (p.v < vMin) vMin = p.v;
      if (p.v > vMax) vMax = p.v;
    }
    for (const a of this.annotations) {
      if (a.pattern !== kind) continue;
      uMin = Math.min(uMin, a.u);
      uMax = Math.max(uMax, a.u);
      vMin = Math.min(vMin, a.v);
      vMax = Math.max(vMax, a.v);
    }
    return { uMin, vMin, uMax, vMax };
  }

  /**
   * Non-empty 1:1 tiles for a pattern, in pattern coordinates — the exact
   * pages the backend will emit (empty tiles are skipped on both sides).
   */
  tiles(kind: PatternKind): { u0: number; vTop: number; w: number; h: number; label: string }[] {
    const box = this.cutBox(kind);
    if (!box || this.scale !== "one_to_one") return [];
    const plan = planTiles(this.page, box.uMax - box.uMin, box.vMax - box.vMin);
    const pts = this.cutPoints(kind);
    const annots = this.annotations.filter((a) => a.pattern === kind);
    const out: { u0: number; vTop: number; w: number; h: number; label: string }[] = [];
    for (let row = 0; row < plan.rows; row++) {
      for (let col = 0; col < plan.cols; col++) {
        const u0 = box.uMin + col * plan.stepX;
        const vTop = box.vMax - row * plan.stepY;
        const m = plan.overlap;
        const inside = (u: number, v: number) =>
          u >= u0 - m && u <= u0 + plan.viewW + m && v >= vTop - plan.viewH - m && v <= vTop + m;
        const has =
          pts.some((p) => inside(p.u, p.v)) || annots.some((a) => inside(a.u, a.v));
        if (has) out.push({ u0, vTop, w: plan.viewW, h: plan.viewH, label: tileLabel(col, row) });
      }
    }
    return out;
  }

  /** Drawing extents of a pattern (mm), circumference included. */
  patternBox(kind: PatternKind): { uMin: number; vMin: number; w: number; h: number } | null {
    const r = store.result;
    if (!r) return null;
    const pts = kind === "branch" ? r.dev_branch : r.dev_main;
    if (!pts || pts.length === 0) return null;
    const circumference =
      kind === "branch" ? (r.circumference_branch ?? r.circumference_main) : r.circumference_main;
    let uMin = pts[0].u,
      uMax = uMin,
      vMin = pts[0].v,
      vMax = vMin;
    for (const p of pts) {
      if (p.u < uMin) uMin = p.u;
      if (p.u > uMax) uMax = p.u;
      if (p.v < vMin) vMin = p.v;
      if (p.v > vMax) vMax = p.v;
    }
    const w = Math.max(uMax - uMin, circumference);
    return { uMin, vMin, w, h: Math.max(vMax - vMin, 0.1) };
  }

  addAnnotation(pattern: PatternKind, u: number, v: number) {
    store.editor.annotations.push({ pattern, u, v, text: "annotation", size_mm: 4 });
    this.selected = store.editor.annotations.length - 1;
    store.persist();
  }

  removeAnnotation(index: number) {
    store.editor.annotations.splice(index, 1);
    if (this.selected === index) this.selected = null;
    else if (this.selected !== null && this.selected > index) this.selected -= 1;
    store.persist();
  }

  async exportPdf() {
    await this.run("pdf", () => exportApi.pdf(this.document(), `cylix-${slug()}.pdf`));
  }

  async exportDxf() {
    await this.run("dxf", () => exportApi.dxf(this.document(), `cylix-${slug()}.dxf`));
  }

  private async run(kind: "pdf" | "dxf", fn: () => Promise<void>) {
    this.busy = kind;
    this.error = null;
    try {
      await fn();
    } catch (e: unknown) {
      this.error = e instanceof Error ? e.message : String(e);
    } finally {
      this.busy = null;
    }
  }
}

function slug(): string {
  return store.active.name
    .toLowerCase()
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-|-$/g, "") || "etude";
}

export const editor = new EditorStore();
