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
  phi_y: number;
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
  phi_y?: number;
  z0?: number;
  n_samples?: number;
}

/** One branch of a multi-branch node ("châssis"). Radians / mm. */
export interface MultiBranchSpec {
  r: number;
  z: number;
  phi: number;
  psi: number;
}

export interface MultiInput {
  r1: number;
  branches: MultiBranchSpec[];
  n_samples?: number;
}

export interface MultiBranchResult {
  r: number;
  z: number;
  phi: number;
  psi: number;
  /** Landing curve on the main tube (open, one period). */
  dev: DevPoint[];
  /** Crossing contours carved by neighbours (closed), same dev plane. */
  holes: HoleResult[];
  curve3d: Point3[];
  bbox: BBox2 | null;
  circumference: number;
  cut_by_neighbor: boolean;
}

export interface HoleResult {
  branch: number;
  pts: DevPoint[];
  closed: boolean;
  bbox: BBox2 | null;
}

/** Joint metrics of a coplanar pair of branches (EN 1993-1-8 quantities). */
export interface NodePair {
  i: number;
  j: number;
  same_side: boolean;
  /** Distance from the chord axis to the brace-axes crossing point, mm. */
  eccentricity: number | null;
  /** Gap between the footprints in the joint plane, mm (null if overlapping). */
  gap: number | null;
  /** Overlap ratio λov, % (null when there is a gap). */
  overlap_pct: number | null;
}

export interface MultiPayload {
  mode: "multi";
  r1: number;
  circumference_main: number;
  branches: MultiBranchResult[];
  holes: HoleResult[];
  pairs: NodePair[];
  warnings: string[];
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
  multi: (input: MultiInput) => postJSON<MultiPayload>("/api/intersect/multi", input),
};
