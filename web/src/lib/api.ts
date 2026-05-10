import { invoke } from "@tauri-apps/api/core";

export type Branch = "outer" | "inner";

export interface DevPoint {
  theta: number;
  u: number;
  v: number;
}

export interface Point3 {
  x: number;
  y: number;
  z: number;
}

export interface BBox2 {
  u_min: number;
  u_max: number;
  v_min: number;
  v_max: number;
}

export interface IntersectionPayload {
  mode: "cyl_cyl" | "cyl_plane";
  r1: number;
  r2: number | null;
  phi: number;
  branch: Branch | null;
  curve3d: Point3[];
  dev_branch: DevPoint[];
  dev_main: DevPoint[] | null;
  bbox_branch: BBox2 | null;
  bbox_main: BBox2 | null;
  circumference_branch: number | null;
  circumference_main: number;
}

export interface CylCylInput {
  r1: number;
  r2: number;
  phi: number;
  n_samples?: number;
  branch?: Branch;
}

export interface CylPlaneInput {
  r1: number;
  phi: number;
  z0?: number;
  n_samples?: number;
}

async function call<T>(cmd: string, input: unknown): Promise<T> {
  try {
    return await invoke<T>(cmd, { input });
  } catch (err: unknown) {
    const msg = typeof err === "string" ? err : err instanceof Error ? err.message : String(err);
    throw new Error(msg);
  }
}

export const api = {
  cylCyl: (input: CylCylInput) =>
    call<IntersectionPayload>("intersect_cyl_cyl", input),
  cylPlane: (input: CylPlaneInput) =>
    call<IntersectionPayload>("intersect_cyl_plane", input),
};
