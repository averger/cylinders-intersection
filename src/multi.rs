//! Multi-branch nodes ("châssis tubulaire") — N inclined tubes landing on one
//! main tube, with mutual first-contact cuts between neighbouring branches.
//!
//! Geometry
//! --------
//! * Main cylinder: axis `Oz`, radius `r1`, equation `x² + y² = r1²`.
//! * Branch `i`: radius `r_i ≤ r1`, axis through `C_i = (0, 0, z_i)`, with
//!   direction `d_i = Rz(ψ_i)·Rx(φ_i)·e_z` — `φ` is the inclination from the
//!   main axis (rotation around `Ox`), `ψ` the azimuth around the main tube.
//! * A surface generatrix of branch `i` at angle `θ` is the straight line
//!   `P(t) = C_i + r_i(cosθ·u_i + sinθ·w_i) + t·d_i` with `(u_i, w_i)` the
//!   rotated `(e_x, e_y)` frame.
//!
//! The physical cut of branch `i` follows the *first contact* rule: coming
//! from `t = +∞`, the tube stops at the first surface met — the main
//! cylinder or any neighbouring branch.  Both intersections are quadratics
//! in `t` (closed form, no numerical solving); the cut is the maximum of the
//! admissible entry roots.  A neighbour contact is only admissible if it
//! happens outside the main cylinder — inside there is no branch material.
//!
//! The development of the main tube carries one opening per branch: the
//! classic gueule de loup of `(r1, r_i, φ_i)` translated by `(r1·ψ_i, z_i)`
//! in the unrolled plane (rotation around the axis unrolls to a horizontal
//! shift, translation along the axis to a vertical one — both exact).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::geometry::{rot_x, rot_z, BBox2, Branch, Point3};
use crate::intersection::{cyl_cyl, CylCylInput, DevPoint};

/// One branch of the node.
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct MultiBranchSpec {
    /// Branch radius, mm (`r ≤ r1`).
    pub r: f64,
    /// Height of the branch axis on the main axis, mm.
    pub z: f64,
    /// Inclination from the main axis, radians, in `(0, π)` — values above
    /// `π/2` lean the branch downwards.
    pub phi: f64,
    /// Azimuth around the main tube, radians.
    #[serde(default)]
    pub psi: f64,
}

/// Input of the multi-branch computation.
#[derive(Debug, Clone, Deserialize)]
pub struct MultiInput {
    /// Main tube radius, mm.
    pub r1: f64,
    pub branches: Vec<MultiBranchSpec>,
    #[serde(default = "default_samples")]
    pub n_samples: usize,
}

fn default_samples() -> usize {
    1440
}

/// Developed template of one branch plus its 3D rim.
#[derive(Debug, Clone, Serialize)]
pub struct MultiBranchResult {
    pub r: f64,
    pub z: f64,
    pub phi: f64,
    pub psi: f64,
    /// Template: `u = r·θ`, `v = t_cut(θ)` — open curve over one period.
    pub dev: Vec<DevPoint>,
    /// 3D rim of the cut (world coordinates, mm).
    pub curve3d: Vec<Point3>,
    pub bbox: Option<BBox2>,
    pub circumference: f64,
    /// True where a neighbouring branch (not the main tube) shapes part of
    /// the cut — the template carries the mutual seam.
    pub cut_by_neighbor: bool,
}

/// One opening on the developed main tube.
#[derive(Debug, Clone, Serialize)]
pub struct HoleResult {
    /// Which branch the opening belongs to (index in `branches`).
    pub branch: usize,
    pub pts: Vec<DevPoint>,
    pub closed: bool,
    pub bbox: Option<BBox2>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiPayload {
    pub mode: &'static str,
    pub r1: f64,
    pub circumference_main: f64,
    pub branches: Vec<MultiBranchResult>,
    /// All openings of the main tube, positioned in its unrolled plane.
    pub holes: Vec<HoleResult>,
    pub warnings: Vec<String>,
}

struct BranchFrame {
    c: Vector3<f64>,
    d: Vector3<f64>,
    u: Vector3<f64>,
    w: Vector3<f64>,
    r: f64,
}

impl BranchFrame {
    fn new(spec: &MultiBranchSpec) -> Self {
        let m = rot_z(spec.psi) * rot_x(spec.phi);
        BranchFrame {
            c: Vector3::new(0.0, 0.0, spec.z),
            d: m * Vector3::z(),
            u: m * Vector3::x(),
            w: m * Vector3::y(),
            r: spec.r,
        }
    }

    /// Base point of the generatrix at angle `θ` (on the branch surface).
    fn base(&self, theta: f64) -> Vector3<f64> {
        let (st, ct) = theta.sin_cos();
        self.c + self.r * (ct * self.u + st * self.w)
    }
}

/// Larger root of `a·t² + b·t + c = 0` — the entry point coming from
/// `t = +∞`.  `None` when the line misses the cylinder.
fn entry_root(a: f64, b: f64, c: f64) -> Option<f64> {
    if a.abs() < 1e-14 {
        return None;
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    Some((-b + disc.sqrt()) / (2.0 * a))
}

/// Entry of the generatrix line `B + t·d` into the main cylinder.
fn entry_into_main(b: &Vector3<f64>, d: &Vector3<f64>, r1: f64) -> Option<f64> {
    let a = d.x * d.x + d.y * d.y;
    let bb = 2.0 * (b.x * d.x + b.y * d.y);
    let c = b.x * b.x + b.y * b.y - r1 * r1;
    entry_root(a, bb, c)
}

/// Entry of the generatrix line `B + t·d` into the (infinite) cylinder of
/// branch `other`: distance to its axis equals its radius.
fn entry_into_branch(b: &Vector3<f64>, d: &Vector3<f64>, other: &BranchFrame) -> Option<f64> {
    let m0 = b - other.c;
    let e = d - other.d * d.dot(&other.d);
    let m = m0 - other.d * m0.dot(&other.d);
    entry_root(e.dot(&e), 2.0 * m.dot(&e), m.dot(&m) - other.r * other.r)
}

/// Compute the whole node: one template per branch (first-contact cut) and
/// the developed main tube with every opening in place.
pub fn multi(input: &MultiInput) -> MultiPayload {
    let n = input.n_samples.max(64);
    let r1 = input.r1;
    let frames: Vec<BranchFrame> = input.branches.iter().map(BranchFrame::new).collect();

    let mut branches = Vec::with_capacity(frames.len());
    let mut neighbor_pairs: Vec<(usize, usize)> = Vec::new();

    for (i, f) in frames.iter().enumerate() {
        let mut dev = Vec::with_capacity(n);
        let mut curve3d = Vec::with_capacity(n);
        let mut cut_by_neighbor = false;

        for k in 0..n {
            let theta = std::f64::consts::TAU * (k as f64) / (n as f64);
            let base = f.base(theta);
            let Some(mut t_cut) = entry_into_main(&base, &f.d, r1) else {
                continue;
            };

            for (j, other) in frames.iter().enumerate() {
                if j == i {
                    continue;
                }
                if let Some(t_j) = entry_into_branch(&base, &f.d, other) {
                    if t_j > t_cut {
                        // Only material outside the main tube can stop us.
                        let p = base + t_j * f.d;
                        if p.x * p.x + p.y * p.y >= r1 * r1 * (1.0 - 1e-9) {
                            t_cut = t_j;
                            cut_by_neighbor = true;
                            if !neighbor_pairs.contains(&(i.min(j), i.max(j))) {
                                neighbor_pairs.push((i.min(j), i.max(j)));
                            }
                        }
                    }
                }
            }

            let p = base + t_cut * f.d;
            curve3d.push(Point3::from(p));
            dev.push(DevPoint { theta, u: f.r * theta, v: t_cut });
        }

        let bbox = BBox2::from_points(dev.iter().map(|p| (p.u, p.v)));
        let spec = &input.branches[i];
        branches.push(MultiBranchResult {
            r: spec.r,
            z: spec.z,
            phi: spec.phi,
            psi: spec.psi,
            dev,
            curve3d,
            bbox,
            circumference: std::f64::consts::TAU * spec.r,
            cut_by_neighbor,
        });
    }

    // Openings of the main tube: the isolated gueule de loup of each branch,
    // shifted by (r1·ψ, z) in the unrolled plane.
    let mut holes = Vec::with_capacity(input.branches.len());
    for (i, spec) in input.branches.iter().enumerate() {
        let iso = cyl_cyl(CylCylInput {
            r1,
            r2: spec.r,
            phi: spec.phi,
            n_samples: n,
            branch: Branch::Outer,
        });
        if let Some(pts) = iso.dev_main {
            let du = r1 * spec.psi;
            let shifted: Vec<DevPoint> = pts
                .into_iter()
                .map(|p| DevPoint { theta: p.theta + spec.psi, u: p.u + du, v: p.v + spec.z })
                .collect();
            let bbox = BBox2::from_points(shifted.iter().map(|p| (p.u, p.v)));
            holes.push(HoleResult {
                branch: i,
                pts: shifted,
                closed: iso.dev_main_closed,
                bbox,
            });
        }
    }

    // Warnings: mutual seams and overlapping openings.
    let mut warnings = Vec::new();
    for &(i, j) in &neighbor_pairs {
        warnings.push(format!(
            "Les piquages {} et {} se rencontrent avant le tube principal — leurs gabarits portent la couture mutuelle.",
            i + 1,
            j + 1
        ));
    }
    let circ = std::f64::consts::TAU * r1;
    for a in 0..holes.len() {
        for b in (a + 1)..holes.len() {
            if let (Some(ba), Some(bb)) = (holes[a].bbox, holes[b].bbox) {
                if bboxes_overlap_on_tube(&ba, &bb, circ) {
                    warnings.push(format!(
                        "Les lumières des piquages {} et {} se chevauchent sur le tube principal — la découpe résultante est l'union des deux contours.",
                        holes[a].branch + 1,
                        holes[b].branch + 1
                    ));
                }
            }
        }
    }

    MultiPayload {
        mode: "multi",
        r1,
        circumference_main: circ,
        branches,
        holes,
        warnings,
    }
}

/// Overlap test between two developed bboxes, `u` being periodic (2πR).
fn bboxes_overlap_on_tube(a: &BBox2, b: &BBox2, circ: f64) -> bool {
    if a.v_min > b.v_max || b.v_min > a.v_max {
        return false;
    }
    // Compare the u intervals modulo the circumference: shift b by the
    // multiple of `circ` that brings its centre closest to a's centre.
    let ca = 0.5 * (a.u_min + a.u_max);
    let cb = 0.5 * (b.u_min + b.u_max);
    let shift = ((ca - cb) / circ).round() * circ;
    let (b0, b1) = (b.u_min + shift, b.u_max + shift);
    a.u_min <= b1 && b0 <= a.u_max
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v_at(dev: &[DevPoint], theta: f64) -> f64 {
        dev.iter()
            .min_by(|a, b| {
                (a.theta - theta)
                    .abs()
                    .partial_cmp(&(b.theta - theta).abs())
                    .unwrap()
            })
            .unwrap()
            .v
    }

    #[test]
    fn single_branch_matches_cyl_cyl_outer() {
        // One branch through the origin with ψ = 0 must reproduce the classic
        // two-cylinder outer development exactly.
        let phi = 0.9;
        let multi_res = multi(&MultiInput {
            r1: 50.0,
            branches: vec![MultiBranchSpec { r: 30.0, z: 0.0, phi, psi: 0.0 }],
            n_samples: 720,
        });
        let iso = cyl_cyl(CylCylInput {
            r1: 50.0,
            r2: 30.0,
            phi,
            n_samples: 720,
            branch: Branch::Outer,
        });
        let dev = &multi_res.branches[0].dev;
        assert_eq!(dev.len(), iso.dev_branch.len());
        for (a, b) in dev.iter().zip(iso.dev_branch.iter()) {
            assert!((a.v - b.v).abs() < 1e-9, "θ={}: {} vs {}", a.theta, a.v, b.v);
        }
        assert!(!multi_res.branches[0].cut_by_neighbor);
        assert!(multi_res.warnings.is_empty());
    }

    #[test]
    fn distant_branches_do_not_interact() {
        let spec = |z: f64| MultiBranchSpec { r: 20.0, z, phi: std::f64::consts::FRAC_PI_2, psi: 0.0 };
        let res = multi(&MultiInput {
            r1: 40.0,
            branches: vec![spec(-200.0), spec(200.0)],
            n_samples: 360,
        });
        assert!(!res.branches[0].cut_by_neighbor);
        assert!(!res.branches[1].cut_by_neighbor);
        let iso = multi(&MultiInput {
            r1: 40.0,
            branches: vec![spec(-200.0)],
            n_samples: 360,
        });
        for (a, b) in res.branches[0].dev.iter().zip(iso.branches[0].dev.iter()) {
            assert!((a.v - b.v).abs() < 1e-12);
        }
    }

    #[test]
    fn v_node_branches_cut_each_other() {
        // Two branches leaning towards each other (φ and π−φ) close enough to
        // meet above the main tube: the mutual seam must lengthen the cut
        // (larger t) somewhere, and never shorten it.
        let phi = std::f64::consts::FRAC_PI_4;
        let mk = |z: f64, phi: f64| MultiBranchSpec { r: 25.0, z, phi, psi: 0.0 };
        let node = multi(&MultiInput {
            r1: 40.0,
            branches: vec![mk(-45.0, phi), mk(45.0, std::f64::consts::PI - phi)],
            n_samples: 720,
        });
        assert!(node.branches[0].cut_by_neighbor, "la branche 1 doit être coupée par la 2");
        assert!(node.branches[1].cut_by_neighbor);
        assert!(!node.warnings.is_empty());

        let alone = multi(&MultiInput {
            r1: 40.0,
            branches: vec![mk(-45.0, phi)],
            n_samples: 720,
        });
        let mut raised = 0usize;
        for (a, b) in node.branches[0].dev.iter().zip(alone.branches[0].dev.iter()) {
            assert!(a.v >= b.v - 1e-9, "la couture ne peut que raccourcir le tube");
            if a.v > b.v + 1e-6 {
                raised += 1;
            }
        }
        assert!(raised > 10, "la couture mutuelle doit modifier une plage de θ ({raised})");
    }

    #[test]
    fn mutual_seam_is_symmetric_for_a_symmetric_v() {
        // Same V-node: by symmetry (z ↔ −z), the two templates must be
        // identical up to the θ ↦ −θ reflection of their seam.
        let phi = std::f64::consts::FRAC_PI_4;
        let mk = |z: f64, phi: f64| MultiBranchSpec { r: 25.0, z, phi, psi: 0.0 };
        let node = multi(&MultiInput {
            r1: 40.0,
            branches: vec![mk(-45.0, phi), mk(45.0, std::f64::consts::PI - phi)],
            n_samples: 720,
        });
        let d0 = &node.branches[0].dev;
        let d1 = &node.branches[1].dev;
        for p in d0.iter().step_by(16) {
            let mirrored = (std::f64::consts::TAU - p.theta).rem_euclid(std::f64::consts::TAU);
            let v1 = v_at(d1, mirrored);
            assert!((p.v - v1).abs() < 1e-6, "θ={} : {} vs {}", p.theta, p.v, v1);
        }
    }

    #[test]
    fn holes_follow_azimuth_and_height() {
        let base = MultiBranchSpec { r: 20.0, z: 0.0, phi: 1.0, psi: 0.0 };
        let moved = MultiBranchSpec { r: 20.0, z: 35.0, phi: 1.0, psi: 0.8 };
        let a = multi(&MultiInput { r1: 45.0, branches: vec![base], n_samples: 720 });
        let b = multi(&MultiInput { r1: 45.0, branches: vec![moved], n_samples: 720 });
        let (ba, bb) = (a.holes[0].bbox.unwrap(), b.holes[0].bbox.unwrap());
        let du = 45.0 * 0.8;
        assert!((bb.u_min - ba.u_min - du).abs() < 1e-9);
        assert!((bb.u_max - ba.u_max - du).abs() < 1e-9);
        assert!((bb.v_min - ba.v_min - 35.0).abs() < 1e-9);
        assert!(a.holes[0].closed && b.holes[0].closed);
    }

    #[test]
    fn overlapping_holes_raise_a_warning() {
        // Two branches at the same height, azimuths 0 and a small angle:
        // their openings overlap on the main tube.
        let mk = |psi: f64| MultiBranchSpec { r: 25.0, z: 0.0, phi: 1.2, psi };
        let res = multi(&MultiInput {
            r1: 40.0,
            branches: vec![mk(0.0), mk(0.5)],
            n_samples: 720,
        });
        assert!(
            res.warnings.iter().any(|w| w.contains("chevauchent")),
            "warnings: {:?}",
            res.warnings
        );
    }
}
