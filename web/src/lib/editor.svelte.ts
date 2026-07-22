/**
 * State of the export studio (the pre-export SVG editor).  Owns everything
 * the user can adjust before producing a PDF / DXF: page setup, layers,
 * title block, annotations — and mirrors it into the backend document model.
 */

import { store } from "./store.svelte";
import {
  exportApi,
  planTiles,
  type Annotation,
  type ExportDocument,
  type Layers,
  type PageSpec,
  type PatternKind,
  type ScaleMode,
  type SourceSpec,
  type TitleBlock,
} from "./export";

function today(): string {
  const d = new Date();
  const mm = String(d.getMonth() + 1).padStart(2, "0");
  const dd = String(d.getDate()).padStart(2, "0");
  return `${d.getFullYear()}-${mm}-${dd}`;
}

class EditorStore {
  page = $state<PageSpec>({ format: "a4", orientation: "landscape", margin_mm: 10 });
  scale = $state<ScaleMode>("one_to_one");
  layers = $state<Layers>({ grid: true, frame: true, axis: true, labels: true, scale_bar: true });
  titleBlock = $state<TitleBlock>({
    title: "",
    project: "",
    author: "",
    date: today(),
    notes: "",
  });
  annotations = $state<Annotation[]>([]);
  cutWidth = $state(0.35);
  /** Pattern currently shown in the preview. */
  active = $state<PatternKind>("branch");
  selected = $state<number | null>(null);
  busy = $state<null | "pdf" | "dxf">(null);
  error = $state<string | null>(null);

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
    return {
      source: this.source(),
      page: { ...this.page },
      scale: this.scale,
      layers: { ...this.layers },
      title_block: { ...this.titleBlock },
      annotations: this.annotations.map((a) => ({ ...a })),
      patterns: this.available,
      cut_width_mm: this.cutWidth,
    };
  }

  /** Total exported page count, matching the backend layout logic. */
  pageCount(): number {
    if (!store.result) return 0;
    if (this.scale === "fit") return this.available.length;
    let pages = 0;
    for (const kind of this.available) {
      const box = this.patternBox(kind);
      if (!box) continue;
      const plan = planTiles(this.page, box.w, box.h);
      pages += plan.cols * plan.rows;
    }
    return pages;
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

  addAnnotation(u: number, v: number) {
    this.annotations.push({
      pattern: this.active,
      u,
      v,
      text: "annotation",
      size_mm: 4,
    });
    this.selected = this.annotations.length - 1;
  }

  removeAnnotation(index: number) {
    this.annotations.splice(index, 1);
    if (this.selected === index) this.selected = null;
    else if (this.selected !== null && this.selected > index) this.selected -= 1;
  }

  async exportPdf() {
    await this.run("pdf", () => exportApi.pdf(this.document()));
  }

  async exportDxf() {
    await this.run("dxf", () => exportApi.dxf(this.document()));
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

export const editor = new EditorStore();
