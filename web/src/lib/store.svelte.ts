import { api, type Branch, type IntersectionPayload } from "./api";

export type Mode = "cyl_cyl" | "cyl_plane";

export interface Params {
  mode: Mode;
  d1: number;       // diameter of main cyl, mm
  d2: number;       // diameter of branch (cyl-cyl only), mm
  angleDeg: number; // angle between axes (cyl-cyl) or plane tilt (cyl-plane)
  branch: Branch;
  z0: number;       // plane offset (cyl-plane), mm
  samples: number;
}

const DEFAULTS: Params = {
  mode: "cyl_cyl",
  d1: 100,
  d2: 70,
  angleDeg: 45,
  branch: "outer",
  z0: 0,
  samples: 1440,
};

class Store {
  params = $state<Params>({ ...DEFAULTS });
  result = $state<IntersectionPayload | null>(null);
  loading = $state(false);
  error = $state<string | null>(null);
  lastComputedAt = $state<number>(0);

  // Monotonic token used to drop stale responses when the user is actively
  // scrubbing a slider.
  private token = 0;

  reset() {
    this.params = { ...DEFAULTS };
  }

  async compute() {
    const myToken = ++this.token;
    this.loading = true;
    this.error = null;
    try {
      const phi = (this.params.angleDeg * Math.PI) / 180;
      let res: IntersectionPayload;
      if (this.params.mode === "cyl_cyl") {
        res = await api.cylCyl({
          r1: this.params.d1 / 2,
          r2: this.params.d2 / 2,
          phi,
          n_samples: this.params.samples,
          branch: this.params.branch,
        });
      } else {
        res = await api.cylPlane({
          r1: this.params.d1 / 2,
          phi,
          z0: this.params.z0,
          n_samples: this.params.samples,
        });
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

export const store = new Store();
