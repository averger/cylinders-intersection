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
  dev_main_closed: boolean;
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

async function postJSON<T>(path: string, body: unknown): Promise<T> {
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
  return (await res.json()) as T;
}

export const api = {
  cylCyl: (input: CylCylInput) =>
    postJSON<IntersectionPayload>("/api/intersect/cyl-cyl", input),
  cylPlane: (input: CylPlaneInput) =>
    postJSON<IntersectionPayload>("/api/intersect/cyl-plane", input),
};
