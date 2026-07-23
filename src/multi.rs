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
    let tau = std::f64::consts::TAU;
    let frames: Vec<BranchFrame> = input.branches.iter().map(BranchFrame::new).collect();

    // Axial floor of each branch: the lowest wall-entry of its generatrices.
    // Below its floor a tube certainly has no material (it would be beyond
    // the far side of the main tube) — contacts there are phantoms of the
    // infinite cylinder.  Above it, the workshop convention applies: each
    // tube is cut on the FULL cylinder surface of its neighbours, so the
    // walls kiss along the seam and never cross inside the joint.
    let floors: Vec<f64> = frames
        .iter()
        .map(|f| {
            (0..n)
                .filter_map(|k| {
                    entry_into_main(&f.base(tau * (k as f64) / (n as f64)), &f.d, r1)
                })
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    let above_floor = |j: usize, p: &Vector3<f64>| -> bool {
        let fr = &frames[j];
        (p - fr.c).dot(&fr.d) >= floors[j] - 1e-6
    };

    // --- First-contact cuts. ----------------------------------------------
    let mut branches = Vec::with_capacity(frames.len());
    let mut neighbor_pairs: Vec<(usize, usize)> = Vec::new();

    for (i, f) in frames.iter().enumerate() {
        let mut dev = Vec::with_capacity(n);
        let mut curve3d = Vec::with_capacity(n);
        let mut cut_by_neighbor = false;

        for k in 0..n {
            let theta = tau * (k as f64) / (n as f64);
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
                        // A neighbour contact only counts outside the main
                        // tube and above the neighbour's axial floor.
                        let p = base + t_j * f.d;
                        if p.x * p.x + p.y * p.y >= r1 * r1 * (1.0 - 1e-9)
                            && above_floor(j, &p)
                        {
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

    // Overlapping openings are merged into their envelope: the template of
    // the main tube must only ever show the actual cut contour.
    let circ = std::f64::consts::TAU * r1;
    let (holes, merged) = merge_overlapping_holes(holes, circ);

    // Warnings: mutual seams and merged openings.
    let mut warnings = Vec::new();
    for &(i, j) in &neighbor_pairs {
        warnings.push(format!(
            "Les piquages {} et {} se rencontrent avant le tube principal — leurs gabarits portent la couture mutuelle.",
            i + 1,
            j + 1
        ));
    }
    for (i, j) in merged {
        warnings.push(format!(
            "Les lumières des piquages {} et {} se chevauchent — la feuille du tube principal porte leur contour d'enveloppe.",
            i + 1,
            j + 1
        ));
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

// ---------------------------------------------------------------------
// Envelope of overlapping openings (2D boolean union of closed loops).
// ---------------------------------------------------------------------

/// Even-odd point-in-polygon test on a closed developed loop.
fn point_in_loop(u: f64, v: f64, pts: &[DevPoint]) -> bool {
    let mut inside = false;
    let n = pts.len();
    for i in 0..n {
        let a = &pts[i];
        let b = &pts[(i + 1) % n];
        if (a.v > v) != (b.v > v) {
            let x = a.u + (v - a.v) / (b.v - a.v) * (b.u - a.u);
            if x > u {
                inside = !inside;
            }
        }
    }
    inside
}

fn mean_step(pts: &[DevPoint]) -> f64 {
    if pts.len() < 2 {
        return 1.0;
    }
    let total: f64 = pts
        .windows(2)
        .map(|w| (w[1].u - w[0].u).hypot(w[1].v - w[0].v))
        .sum();
    total / (pts.len() - 1) as f64
}

/// Boolean union of two dense closed loops: keep the runs of each loop whose
/// segments are not strictly inside the other, then stitch them back into a
/// single closed contour by endpoint proximity.
fn union_two_loops(a: &[DevPoint], b: &[DevPoint]) -> Option<Vec<DevPoint>> {
    let runs_outside = |own: &[DevPoint], other: &[DevPoint]| -> Vec<Vec<DevPoint>> {
        let n = own.len();
        let keep: Vec<bool> = (0..n)
            .map(|i| {
                let p = &own[i];
                let q = &own[(i + 1) % n];
                !point_in_loop(0.5 * (p.u + q.u), 0.5 * (p.v + q.v), other)
            })
            .collect();
        if keep.iter().all(|&k| k) {
            return vec![own.to_vec()];
        }
        // Start each run just after a dropped segment.
        let mut runs = Vec::new();
        let start = (0..n).find(|&i| !keep[i]).unwrap();
        let mut current: Vec<DevPoint> = Vec::new();
        for off in 1..=n {
            let i = (start + off) % n;
            if keep[i] {
                if current.is_empty() {
                    current.push(own[i]);
                }
                current.push(own[(i + 1) % n]);
            } else if !current.is_empty() {
                runs.push(std::mem::take(&mut current));
            }
        }
        if !current.is_empty() {
            runs.push(current);
        }
        runs
    };

    let mut runs = runs_outside(a, b);
    runs.extend(runs_outside(b, a));
    runs.retain(|r| r.len() >= 2);
    if runs.is_empty() {
        return None;
    }
    let tol = 6.0 * mean_step(a).max(mean_step(b));

    // Stitch: grow a contour by appending the run whose endpoint matches.
    let mut contour = runs.swap_remove(0);
    while !runs.is_empty() {
        let end = *contour.last().unwrap();
        let mut best: Option<(usize, bool, f64)> = None;
        for (idx, r) in runs.iter().enumerate() {
            let df = (r[0].u - end.u).hypot(r[0].v - end.v);
            let dl = (r[r.len() - 1].u - end.u).hypot(r[r.len() - 1].v - end.v);
            let (rev, d) = if df <= dl { (false, df) } else { (true, dl) };
            if best.is_none() || d < best.unwrap().2 {
                best = Some((idx, rev, d));
            }
        }
        let (idx, rev, d) = best.unwrap();
        if d > tol {
            break; // disconnected leftover (shouldn't happen on dense loops)
        }
        let mut r = runs.swap_remove(idx);
        if rev {
            r.reverse();
        }
        contour.extend(r);
    }
    Some(contour)
}

/// Merge every cluster of overlapping closed openings into its envelope.
/// Returns the new hole list plus the merged index pairs (for warnings).
fn merge_overlapping_holes(
    mut holes: Vec<HoleResult>,
    circ: f64,
) -> (Vec<HoleResult>, Vec<(usize, usize)>) {
    let mut merged_pairs = Vec::new();
    'outer: loop {
        for i in 0..holes.len() {
            for j in (i + 1)..holes.len() {
                if !(holes[i].closed && holes[j].closed) {
                    continue;
                }
                let (Some(bi), Some(bj)) = (holes[i].bbox, holes[j].bbox) else {
                    continue;
                };
                if !bboxes_overlap_on_tube(&bi, &bj, circ) {
                    continue;
                }
                // Bring j into i's unwrap period before the real overlap test.
                let shift = (((bi.u_min + bi.u_max) - (bj.u_min + bj.u_max)) / 2.0 / circ).round()
                    * circ;
                let pts_j: Vec<DevPoint> = holes[j]
                    .pts
                    .iter()
                    .map(|p| DevPoint { theta: p.theta, u: p.u + shift, v: p.v })
                    .collect();
                let really_overlaps = pts_j
                    .iter()
                    .any(|p| point_in_loop(p.u, p.v, &holes[i].pts))
                    || holes[i]
                        .pts
                        .iter()
                        .any(|p| point_in_loop(p.u, p.v, &pts_j));
                if !really_overlaps {
                    continue;
                }
                if let Some(contour) = union_two_loops(&holes[i].pts, &pts_j) {
                    merged_pairs.push((holes[i].branch, holes[j].branch));
                    let branch = holes[i].branch.min(holes[j].branch);
                    let bbox = BBox2::from_points(contour.iter().map(|p| (p.u, p.v)));
                    holes[i] = HoleResult { branch, pts: contour, closed: true, bbox };
                    holes.remove(j);
                    continue 'outer;
                }
            }
        }
        break;
    }
    (holes, merged_pairs)
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

    /// Shoelace area of a closed developed loop.
    fn loop_area(pts: &[DevPoint]) -> f64 {
        let n = pts.len();
        let mut s = 0.0;
        for i in 0..n {
            let a = &pts[i];
            let b = &pts[(i + 1) % n];
            s += a.u * b.v - b.u * a.v;
        }
        (s / 2.0).abs()
    }

    #[test]
    fn overlapping_holes_are_merged_into_their_envelope() {
        // Two branches at the same height, azimuths 0 and 0.5 rad: their
        // openings overlap — the main-tube sheet must carry ONE envelope
        // contour, not two crossing loops.
        let mk = |psi: f64| MultiBranchSpec { r: 25.0, z: 0.0, phi: 1.2, psi };
        let res = multi(&MultiInput {
            r1: 40.0,
            branches: vec![mk(0.0), mk(0.5)],
            n_samples: 720,
        });
        assert_eq!(res.holes.len(), 1, "enveloppe unique attendue");
        assert!(res.holes[0].closed);
        assert!(
            res.warnings.iter().any(|w| w.contains("enveloppe")),
            "warnings: {:?}",
            res.warnings
        );

        // Envelope sanity vs the two isolated openings: same overall bbox,
        // and an area strictly between max(A, B) and A + B (the lens of the
        // overlap is counted once, not twice).
        let iso = |psi: f64| {
            multi(&MultiInput { r1: 40.0, branches: vec![mk(psi)], n_samples: 720 }).holes[0]
                .clone()
        };
        let (a, b) = (iso(0.0), iso(0.5));
        let (ba, bb) = (a.bbox.unwrap(), b.bbox.unwrap());
        let be = res.holes[0].bbox.unwrap();
        assert!((be.u_min - ba.u_min.min(bb.u_min)).abs() < 1e-6);
        assert!((be.u_max - ba.u_max.max(bb.u_max)).abs() < 1e-6);
        assert!((be.v_min - ba.v_min.min(bb.v_min)).abs() < 1e-6);
        assert!((be.v_max - ba.v_max.max(bb.v_max)).abs() < 1e-6);
        let (area_a, area_b) = (loop_area(&a.pts), loop_area(&b.pts));
        let area_e = loop_area(&res.holes[0].pts);
        assert!(area_e > area_a.max(area_b) * 1.01, "{area_e} vs {area_a}/{area_b}");
        assert!(area_e < (area_a + area_b) * 0.999, "la lentille doit être comptée une fois");
    }

    #[test]
    fn disjoint_holes_are_left_untouched() {
        let mk = |psi: f64| MultiBranchSpec { r: 15.0, z: 0.0, phi: 1.2, psi };
        let res = multi(&MultiInput {
            r1: 45.0,
            branches: vec![mk(0.0), mk(std::f64::consts::PI)],
            n_samples: 720,
        });
        assert_eq!(res.holes.len(), 2);
        assert!(res.warnings.is_empty());
    }

    /// Distance from a point to a closed polyline (segment-wise).
    fn dist_to_polyline(u: f64, v: f64, pts: &[DevPoint]) -> f64 {
        let n = pts.len();
        let mut best = f64::INFINITY;
        for i in 0..n {
            let a = &pts[i];
            let b = &pts[(i + 1) % n];
            let (dx, dy) = (b.u - a.u, b.v - a.v);
            let len2 = dx * dx + dy * dy;
            let t = if len2 > 0.0 {
                (((u - a.u) * dx + (v - a.v) * dy) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let (px, py) = (a.u + t * dx, a.v + t * dy);
            best = best.min((u - px).hypot(v - py));
        }
        best
    }

    #[test]
    fn rim_on_main_lies_on_the_envelope_contour() {
        // What the eye checks in the 3D view: where a branch lands on the
        // main tube, its rim must sit ON the opening contour of the main
        // sheet — no wall may land inside the cut-away region.
        // Sane V-node: axes crossing INSIDE the main tube (|z| < r1/tanφ·…),
        // so both tubes land on the wall around their mutual seam.
        let phi = std::f64::consts::FRAC_PI_4;
        let node = multi(&MultiInput {
            r1: 50.0,
            branches: vec![
                MultiBranchSpec { r: 30.0, z: -40.0, phi, psi: 0.0 },
                MultiBranchSpec { r: 30.0, z: 40.0, phi: std::f64::consts::PI - phi, psi: 0.0 },
            ],
            n_samples: 720,
        });
        assert_eq!(node.holes.len(), 1, "les deux lumières fusionnent en enveloppe");
        let env = &node.holes[0].pts;
        for br in &node.branches {
            let mut checked = 0usize;
            let mut far = 0usize;
            for p in &br.curve3d {
                let dist_main = (p.x * p.x + p.y * p.y).sqrt();
                if (dist_main - 50.0).abs() > 1e-6 {
                    continue; // seam point, off the wall
                }
                let u = 50.0 * p.y.atan2(p.x);
                if dist_to_polyline(u, p.z, env) > 0.05 {
                    far += 1;
                }
                checked += 1;
            }
            assert!(checked > 100, "trop peu de points de rive sur le tube ({checked})");
            // Tolerate a handful of crotch landings inside the envelope
            // (wall already cut away by the neighbour's bore right there).
            assert!(
                far * 20 < checked,
                "{far}/{checked} points de rive hors du contour d'enveloppe"
            );
        }
    }

    #[test]
    fn every_rim_point_lies_exactly_on_its_cutting_surface() {
        // Robustness of the matrix pipeline: every 3D rim point of a V-node
        // template must sit either ON the main cylinder (distance to Oz
        // = r1) or ON the neighbour cylinder (distance to its axis = r_j),
        // to numerical precision.
        let phi = std::f64::consts::FRAC_PI_4;
        let specs = [
            MultiBranchSpec { r: 25.0, z: -45.0, phi, psi: 0.3 },
            MultiBranchSpec { r: 25.0, z: 45.0, phi: std::f64::consts::PI - phi, psi: 0.3 },
        ];
        let node = multi(&MultiInput { r1: 40.0, branches: specs.to_vec(), n_samples: 720 });

        let frame = |s: &MultiBranchSpec| {
            let m = crate::geometry::rot_z(s.psi) * crate::geometry::rot_x(s.phi);
            (nalgebra::Vector3::new(0.0, 0.0, s.z), m * nalgebra::Vector3::z())
        };
        for (i, br) in node.branches.iter().enumerate() {
            let other = &specs[1 - i];
            let (c_j, d_j) = frame(other);
            for p in &br.curve3d {
                let v = nalgebra::Vector3::new(p.x, p.y, p.z);
                let dist_main = (v.x * v.x + v.y * v.y).sqrt();
                let rel = v - c_j;
                let dist_axis_j = (rel - d_j * rel.dot(&d_j)).norm();
                let on_main = (dist_main - 40.0).abs() < 1e-9;
                let on_neighbor = (dist_axis_j - other.r).abs() < 1e-9;
                assert!(
                    on_main || on_neighbor,
                    "point hors surface : d_main = {dist_main}, d_axe_voisin = {dist_axis_j}"
                );
            }
        }
    }
}
