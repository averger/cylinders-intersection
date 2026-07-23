/**
 * Export document model — the exact mirror of the backend `ExportDocument`
 * (src/export/mod.rs). The SVG editor edits this structure; the backend is
 * the single source of truth that recomputes the geometry and renders the
 * final PDF / DXF from it.
 */

import type { Branch, MultiBranchSpec } from "./api";

export type SourceSpec =
  | { mode: "cyl_cyl"; r1: number; r2: number; phi: number; n_samples: number; branch: Branch }
  | { mode: "cyl_plane"; r1: number; phi: number; phi_y: number; z0: number; n_samples: number }
  | { mode: "multi"; r1: number; branches: MultiBranchSpec[]; n_samples: number };

export type PageFormat = "a4" | "a3" | "a2";
export type PageOrientation = "portrait" | "landscape";
export type ScaleMode = "fit" | "one_to_one";
export type PatternKind = "branch" | "main";

export interface PageSpec {
  format: PageFormat;
  orientation: PageOrientation;
  margin_mm: number;
}

export interface Layers {
  grid: boolean;
  frame: boolean;
  axis: boolean;
  labels: boolean;
  scale_bar: boolean;
}

export interface TitleBlock {
  title: string;
  project: string;
  author: string;
  date: string;
  notes: string;
}

export interface Annotation {
  pattern: PatternKind;
  u: number;
  v: number;
  text: string;
  size_mm: number;
}

export interface ExportDocument {
  source: SourceSpec;
  page: PageSpec;
  scale: ScaleMode;
  layers: Layers;
  title_block: TitleBlock;
  annotations: Annotation[];
  patterns: PatternKind[];
  cut_width_mm: number;
}

/** Portrait dimensions of the ISO formats, mm. */
export const PAGE_DIMS: Record<PageFormat, [number, number]> = {
  a4: [210, 297],
  a3: [297, 420],
  a2: [420, 594],
};

/** Page size in mm with orientation applied. */
export function pageSize(page: PageSpec): [number, number] {
  const [w, h] = PAGE_DIMS[page.format];
  return page.orientation === "portrait" ? [w, h] : [h, w];
}

/** Must stay in sync with TITLE_BLOCK_H_MM / TILE_OVERLAP_MM (mod.rs). */
export const TITLE_BLOCK_H_MM = 24;
export const TILE_OVERLAP_MM = 12;

export interface TilePlan {
  cols: number;
  rows: number;
  stepX: number;
  stepY: number;
  viewW: number;
  viewH: number;
  overlap: number;
}

/** TS twin of `plan_tiles` — used to preview the 1:1 page tiling. */
export function planTiles(page: PageSpec, w: number, h: number): TilePlan {
  const [pw, ph] = pageSize(page);
  const viewW = Math.max(pw - 2 * page.margin_mm, 20);
  const viewH = Math.max(ph - 2 * page.margin_mm - TITLE_BLOCK_H_MM, 20);
  const stepX = Math.max(viewW - TILE_OVERLAP_MM, 10);
  const stepY = Math.max(viewH - TILE_OVERLAP_MM, 10);
  const cols = w <= viewW ? 1 : 1 + Math.ceil((w - viewW) / stepX);
  const rows = h <= viewH ? 1 : 1 + Math.ceil((h - viewH) / stepY);
  return { cols, rows, stepX, stepY, viewW, viewH, overlap: TILE_OVERLAP_MM };
}

/** Spreadsheet-style tile label, matching the backend ("A1", "B3", …). */
export function tileLabel(col: number, row: number): string {
  let letters = "";
  let c = col;
  for (;;) {
    letters = String.fromCharCode(65 + (c % 26)) + letters;
    if (c < 26) break;
    c = Math.floor(c / 26) - 1;
  }
  return `${letters}${row + 1}`;
}

async function postForBlob(path: string, body: unknown): Promise<Blob> {
  const res = await fetch(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    let msg = `HTTP ${res.status}`;
    try {
      const j = await res.json();
      if (j && typeof j === "object" && "error" in j) msg = String(j.error);
    } catch {
      /* ignore */
    }
    throw new Error(msg);
  }
  return await res.blob();
}

function saveBlob(filename: string, blob: Blob) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 0);
}

export const exportApi = {
  async pdf(doc: ExportDocument, filename = "gabarits.pdf") {
    saveBlob(filename, await postForBlob("/api/export/pdf", doc));
  },
  async dxf(doc: ExportDocument, filename = "gabarits.dxf") {
    saveBlob(filename, await postForBlob("/api/export/dxf", doc));
  },
  async stl(doc: ExportDocument, filename = "maquette.stl") {
    saveBlob(filename, await postForBlob("/api/export/stl", doc));
  },
};
