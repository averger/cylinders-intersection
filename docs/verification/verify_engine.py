# Cross-validation: Rust engine (HTTP API) vs the original NumPy reference.
import json
import sys
import urllib.request

import numpy as np

sys.path.insert(0, "/home/user/cylinders-intersection")
from tubes_intersection_and_unwrap import intersect_cylinders, unwrap_on_cylinder1

BASE = "http://127.0.0.1:8787"

def api(path, payload):
    req = urllib.request.Request(
        BASE + path,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as r:
        return json.load(r)

worst = 0.0
cases = [
    (50.0, 35.0, np.deg2rad(45.0), "outer"),
    (50.0, 35.0, np.deg2rad(45.0), "inner"),
    (80.0, 80.0, np.deg2rad(90.0), "outer"),
    (60.0, 25.0, np.deg2rad(20.0), "outer"),
    (30.0, 45.0, np.deg2rad(60.0), "outer"),
    (100.0, 99.5, np.deg2rad(75.0), "inner"),
]
n = 720
for r1, r2, ph, branch in cases:
    rust = api("/api/intersect/cyl-cyl", {"r1": r1, "r2": r2, "phi": ph, "n_samples": n, "branch": branch})
    ref = intersect_cylinders(R1=r1, R2=r2, phi=ph, n_samples=n, branch=branch)

    # 3D curve: same θ sampling (0..2π, endpoint=False) — align valid samples.
    ref_xyz = ref.xyz[ref.valid_mask]
    rust_xyz = np.array([[p["x"], p["y"], p["z"]] for p in rust["curve3d"]])
    assert len(ref_xyz) == len(rust_xyz), (len(ref_xyz), len(rust_xyz))
    d3 = np.max(np.linalg.norm(ref_xyz - rust_xyz, axis=1))

    # branch development
    ref_uv = ref.uv[ref.valid_mask]
    rust_uv = np.array([[p["u"], p["v"]] for p in rust["dev_branch"]])
    d2 = np.max(np.linalg.norm(ref_uv - rust_uv, axis=1))

    # gueule de loup: compare as point sets (orders differ: Rust keeps θ order,
    # the reference sorts by u) — use nearest-neighbour distance both ways.
    u1, v1 = unwrap_on_cylinder1(ref, r1)
    ref_m = np.stack([u1, v1], axis=1)
    rust_m = np.array([[p["u"], p["v"]] for p in rust["dev_main"]])
    def set_dist(A, B):
        return max(
            np.max([np.min(np.linalg.norm(B - a, axis=1)) for a in A[::7]]),
            np.max([np.min(np.linalg.norm(A - b, axis=1)) for b in B[::7]]),
        )
    dm = set_dist(ref_m, rust_m)

    worst = max(worst, d3, d2, dm)
    print(f"R1={r1:6.1f} R2={r2:6.1f} phi={np.rad2deg(ph):5.1f}° {branch:5s} | "
          f"3D {d3:.2e}  dev2 {d2:.2e}  gueule {dm:.2e} mm")

print(f"\nécart maximal Rust vs NumPy : {worst:.3e} mm")
assert worst < 1e-9, "divergence moteur/référence"
print("Moteur Rust identique à la référence NumPy (à l'epsilon machine).")
