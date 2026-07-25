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
//! Every tube ENDS at the node, arriving from the outside.  Each
//! generatrix keeps its outermost free segment: the lower lip is the LAST
//! obstacle met coming in — `t_lip(θ) = max(t_wall, max_j hi_j)`, with
//! `t_wall` the landing on the main tube (`t⁺` of the line–cylinder
//! quadratic, closed form) and `hi_j` the exit root of a higher-priority
//! neighbour's quadratic.  Priority is the list order: branch 1 is only
//! cut by the main tube, branch 2 rides in a saddle over branch 1's back
//! where they overlap (the through member of a K/KT overlap joint), and
//! so on — the lip always rests on material that exists, so the joint
//! closes without any tube running through another.
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
///
/// The template is bounded below by ONE continuous lip (`dev`): the
/// landing sinusoid where the tube reaches the main wall, lifted into a
/// saddle where it rides over a higher-priority branch — every tube ends
/// at the node, and the lip always rests on an intact wall.
#[derive(Debug, Clone, Serialize)]
pub struct MultiBranchResult {
    pub r: f64,
    pub z: f64,
    pub phi: f64,
    pub psi: f64,
    /// Lower lip of the tube: `u = r·θ`, `v = t_lip(θ)` — open curve over
    /// one period, `t_lip = max(t_wall, exit of prioritized neighbours)`.
    pub dev: Vec<DevPoint>,
    /// Always empty for branches (kept for payload stability — openings
    /// only exist on the main tube).
    pub holes: Vec<HoleResult>,
    /// 3D rim of the lip (world coordinates, mm).
    pub curve3d: Vec<Point3>,
    pub bbox: Option<BBox2>,
    pub circumference: f64,
    /// True when at least one neighbour crosses this branch.
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

/// Fabrication metrics of a **coplanar pair** of branches (K / N / X joint):
/// where their axes really cross, and how their footprints sit relative to
/// each other on the chord.  These are the numbers a design office checks
/// (EN 1993-1-8): eccentricity, gap, overlap ratio.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct NodePair {
    pub i: usize,
    pub j: usize,
    /// `true` when both branches leave on the same side of the chord (K / N
    /// joint); `false` when they are opposed through it (X joint).
    pub same_side: bool,
    /// Signed distance from the chord axis to the point where the two brace
    /// axes cross, mm.  `0` = concurrent axes (no eccentricity); positive on
    /// the side the branches come from.  `None` when the axes are parallel.
    pub eccentricity: Option<f64>,
    /// Gap between the two footprints, measured along the generatrix of the
    /// joint plane (mm).  `None` when they overlap or are opposed.
    pub gap: Option<f64>,
    /// Overlap ratio λov = q/p, %, measured in the joint plane on the
    /// footprint of the overlapping (lower-priority) branch taken alone.
    /// `None` when there is a gap.
    pub overlap_pct: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiPayload {
    pub mode: &'static str,
    pub r1: f64,
    pub circumference_main: f64,
    pub branches: Vec<MultiBranchResult>,
    /// All openings of the main tube, positioned in its unrolled plane.
    pub holes: Vec<HoleResult>,
    /// Joint metrics of every coplanar pair of branches.
    pub pairs: Vec<NodePair>,
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

/// Interval `t ∈ (t_lo, t_hi)` where the generatrix line `B + t·d` runs
/// INSIDE the (infinite) cylinder of branch `other` — both roots of the
/// distance-to-axis quadratic.  `None` when the line misses the cylinder.
fn branch_interval(
    b: &Vector3<f64>,
    d: &Vector3<f64>,
    other: &BranchFrame,
) -> Option<(f64, f64)> {
    let m0 = b - other.c;
    let e = d - other.d * d.dot(&other.d);
    let m = m0 - other.d * m0.dot(&other.d);
    let (a, bb, c) = (e.dot(&e), 2.0 * m.dot(&e), m.dot(&m) - other.r * other.r);
    if a.abs() < 1e-14 {
        return None;
    }
    let disc = bb * bb - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    Some(((-bb - sq) / (2.0 * a), (-bb + sq) / (2.0 * a)))
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
    // Axial ceiling: the far end of the tube as modelled (and exported) —
    // beyond it the neighbour has no material either.
    let ceilings: Vec<f64> = frames
        .iter()
        .map(|f| {
            let t_hi = (0..n)
                .filter_map(|k| {
                    entry_into_main(&f.base(tau * (k as f64) / (n as f64)), &f.d, r1)
                })
                .fold(f64::NEG_INFINITY, f64::max);
            t_hi + (3.0 * f.r).max(90.0)
        })
        .collect();
    let within_material = |j: usize, p: &Vector3<f64>| -> bool {
        let fr = &frames[j];
        let s = (p - fr.c).dot(&fr.d);
        s >= floors[j] - 1e-6 && s <= ceilings[j] + 1e-6
    };

    // --- Every tube ends at the node, arriving FROM the outside: each
    // generatrix keeps only its outermost free segment, so its lower lip
    // is the LAST obstacle met coming in — the main-tube wall, or the
    // BACK of a higher-priority branch it overlaps (the through member of
    // a K/KT overlap joint).  t_lip(θ) = max(t_wall, max_j hi_j) is one
    // continuous curve: landing sinusoid where the tube reaches the wall,
    // saddle where it rides on a prioritized neighbour.  The lip always
    // rests ON material that exists, so the joint closes, and nothing
    // ever continues through a neighbour.
    let mut branches = Vec::with_capacity(frames.len());
    let mut neighbor_pairs: Vec<(usize, usize)> = Vec::new();

    for (i, f) in frames.iter().enumerate() {
        let mut dev = Vec::with_capacity(n);
        let mut curve3d = Vec::with_capacity(n);
        let mut cut = false;

        for k in 0..n {
            let theta = tau * (k as f64) / (n as f64);
            let base = f.base(theta);
            let Some(t_wall) = entry_into_main(&base, &f.d, r1) else {
                continue;
            };
            // Priority: the list order decides who rides on whom — branch i
            // only stops on HIGHER-priority branches (j < i), whose walls
            // are never carved by i.  Taking the max exit point makes the
            // pass order-independent: the lip settles on the farthest
            // obstacle, i.e. the first one met coming from outside.
            let mut t_lip = t_wall;
            let mut on: Option<usize> = None;
            for (j, other) in frames.iter().enumerate().take(i) {
                if let Some((_, hi)) = branch_interval(&base, &f.d, other) {
                    if hi <= t_lip + 1e-9 {
                        continue; // crossing wholly below the current lip
                    }
                    let exit = base + hi * f.d;
                    if !within_material(j, &exit) {
                        continue; // phantom infinite-cylinder extension
                    }
                    t_lip = hi;
                    on = Some(j);
                }
            }
            if let Some(j) = on {
                cut = true;
                if !neighbor_pairs.contains(&(j, i)) {
                    neighbor_pairs.push((j, i));
                }
            }
            curve3d.push(Point3::from(base + t_lip * f.d));
            dev.push(DevPoint { theta, u: f.r * theta, v: t_lip });
        }

        let bbox = BBox2::from_points(dev.iter().map(|p| (p.u, p.v)));
        let spec = &input.branches[i];
        branches.push(MultiBranchResult {
            r: spec.r,
            z: spec.z,
            phi: spec.phi,
            psi: spec.psi,
            dev,
            cut_by_neighbor: cut,
            holes: Vec::new(),
            curve3d,
            bbox,
            circumference: std::f64::consts::TAU * spec.r,
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

    // Joint metrics are read on the ISOLATED footprints (λov is defined
    // "in the absence of the overlapped brace"), so before merging.
    let circ = std::f64::consts::TAU * r1;
    let pairs = node_pairs(input, &holes, r1, circ);

    // Overlapping openings are merged into their envelope: the template of
    // the main tube must only ever show the actual cut contour.
    let (holes, merged) = merge_overlapping_holes(holes, circ);

    // Warnings: mutual seams and merged openings.
    let mut warnings = Vec::new();
    for &(through, rider) in &neighbor_pairs {
        warnings.push(format!(
            "Le piquage {} chevauche le piquage {} (prioritaire) — sa lèvre s'arrête en selle sur son dos, la couture épouse son flanc.",
            rider + 1,
            through + 1
        ));
    }
    for (i, j) in merged {
        warnings.push(format!(
            "Les lumières des piquages {} et {} se chevauchent — la feuille du tube principal porte leur contour d'enveloppe.",
            i + 1,
            j + 1
        ));
    }

    // Design checks on the joints.
    let d0 = 2.0 * r1;
    for p in &pairs {
        if let Some(e) = p.eccentricity {
            if e.abs() > 0.25 * d0 {
                warnings.push(format!(
                    "Piquages {} et {} : excentrement e = {:.1} mm, au-delà de 0,25·Ø₁ ({:.1} mm) — hors du domaine courant de l'EN 1993-1-8, le moment secondaire doit être repris par le calcul.",
                    p.i + 1,
                    p.j + 1,
                    e,
                    0.25 * d0
                ));
            }
        }
        if let Some(ov) = p.overlap_pct {
            if ov < 25.0 {
                warnings.push(format!(
                    "Piquages {} et {} : recouvrement λov = {:.0} %, sous le minimum de 25 % (EN 1993-1-8) — soit augmenter le recouvrement, soit passer à un nœud avec jeu.",
                    p.i + 1,
                    p.j + 1,
                    ov
                ));
            }
        }
    }

    MultiPayload {
        mode: "multi",
        r1,
        circumference_main: circ,
        branches,
        holes,
        pairs,
        warnings,
    }
}

/// Axial extent `[v_min, v_max]` of a closed developed loop on the generatrix
/// of abscissa `u0` (periodic, `circ`).  `None` when the loop misses it.
fn span_at_u(pts: &[DevPoint], u0: f64, circ: f64) -> Option<(f64, f64)> {
    for k in -2i32..=2 {
        let u = u0 + f64::from(k) * circ;
        let mut vs: Vec<f64> = Vec::new();
        let m = pts.len();
        for idx in 0..m {
            let a = &pts[idx];
            let b = &pts[(idx + 1) % m];
            // Half-open rule: a vertex sitting exactly ON the generatrix
            // (θ = 0 does) is counted once, never twice, never zero times.
            if (a.u <= u) != (b.u <= u) {
                let s = (u - a.u) / (b.u - a.u);
                vs.push(a.v + s * (b.v - a.v));
            }
        }
        if vs.len() >= 2 {
            return Some((
                vs.iter().cloned().fold(f64::INFINITY, f64::min),
                vs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            ));
        }
    }
    None
}

/// Joint metrics of every coplanar pair of branches.
///
/// Two branches sharing an azimuth (or opposed by π) lie in one plane with
/// the chord axis: their axes then really cross, at a distance from the
/// chord axis that closed form gives as
/// `e = (z_j − z_i)·sinφ_i·sinφ_j / sin(φ_j − φ_i)` on the same side, the
/// denominator becoming `sin(φ_i + φ_j)` when the branches are opposed: this
/// is the eccentricity of the node, zero when they are concurrent.
///
/// Non-coplanar (spatial) pairs are skipped, the plane framework of the
/// standard defining no gap for them.
fn node_pairs(input: &MultiInput, holes: &[HoleResult], r1: f64, circ: f64) -> Vec<NodePair> {
    let tau = std::f64::consts::TAU;
    let mut out = Vec::new();
    let n = input.branches.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (&input.branches[i], &input.branches[j]);
            let dpsi = (b.psi - a.psi).rem_euclid(tau);
            let same_side = dpsi < 1e-6 || (tau - dpsi) < 1e-6;
            let opposed = (dpsi - std::f64::consts::PI).abs() < 1e-6;
            if !same_side && !opposed {
                continue; // spatial joint: no plane definition of the gap
            }
            let denom = if same_side { (b.phi - a.phi).sin() } else { (a.phi + b.phi).sin() };
            let eccentricity = if denom.abs() < 1e-9 {
                None // parallel axes: they never cross
            } else {
                Some((b.z - a.z) * a.phi.sin() * b.phi.sin() / denom)
            };

            // Gap / overlap along the generatrix of the joint plane.  Only
            // meaningful on the same side: opposed branches land on opposite
            // faces of the chord and can never touch.
            let (mut gap, mut overlap_pct) = (None, None);
            if same_side {
                // Crown generatrix of the joint plane: the branch leans along
                // n(ψ) = (sinψ, −cosψ, 0), i.e. azimuth α = ψ − π/2, so its
                // footprint straddles u = r1·(ψ − π/2) on the development.
                let u0 = r1 * (a.psi - std::f64::consts::FRAC_PI_2);
                let find = |idx: usize| holes.iter().find(|h| h.branch == idx);
                if let (Some(hi), Some(hj)) = (find(i), find(j)) {
                    if let (Some(si), Some(sj)) =
                        (span_at_u(&hi.pts, u0, circ), span_at_u(&hj.pts, u0, circ))
                    {
                        let q = si.1.min(sj.1) - si.0.max(sj.0);
                        if q > 0.0 {
                            // λov is read on the overlapping brace, i.e. the
                            // lower-priority one (j), taken alone.
                            let p = sj.1 - sj.0;
                            if p > 1e-9 {
                                overlap_pct = Some(100.0 * q / p);
                            }
                        } else {
                            gap = Some(-q);
                        }
                    }
                }
            }
            out.push(NodePair { i, j, same_side, eccentricity, gap, overlap_pct });
        }
    }
    out
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

    /// Map a developed point of a branch back to 3D world coordinates.
    fn to_world(spec: &MultiBranchSpec, p: &DevPoint) -> nalgebra::Vector3<f64> {
        let m = crate::geometry::rot_z(spec.psi) * crate::geometry::rot_x(spec.phi);
        let theta = p.u / spec.r;
        let base = nalgebra::Vector3::new(0.0, 0.0, spec.z)
            + spec.r
                * (theta.cos() * (m * nalgebra::Vector3::x())
                    + theta.sin() * (m * nalgebra::Vector3::y()));
        base + p.v * (m * nalgebra::Vector3::z())
    }

    fn v_specs() -> Vec<MultiBranchSpec> {
        let phi = std::f64::consts::FRAC_PI_4;
        vec![
            MultiBranchSpec { r: 25.0, z: -45.0, phi, psi: 0.0 },
            MultiBranchSpec { r: 25.0, z: 45.0, phi: std::f64::consts::PI - phi, psi: 0.0 },
        ]
    }

    /// Distance from a world point to the axis of `spec`.
    fn dist_to_axis(spec: &MultiBranchSpec, w: &nalgebra::Vector3<f64>) -> f64 {
        let m = crate::geometry::rot_z(spec.psi) * crate::geometry::rot_x(spec.phi);
        let d = m * nalgebra::Vector3::z();
        let rel = w - nalgebra::Vector3::new(0.0, 0.0, spec.z);
        (rel - d * rel.dot(&d)).norm()
    }

    #[test]
    fn crossing_branches_carve_each_other_and_keep_their_landing() {
        // Overlap model: branch 1 (through member) is untouched, branch 2
        // keeps its full circumference but its lip lifts into a saddle over
        // branch 1's back where they overlap.  Wherever the lip is NOT
        // lifted, the landing is IDENTICAL to the isolated computation.
        let specs = v_specs();
        let node = multi(&MultiInput { r1: 40.0, branches: specs.clone(), n_samples: 720 });
        assert!(node.branches[1].cut_by_neighbor && !node.branches[0].cut_by_neighbor);
        assert!(!node.warnings.is_empty());
        assert_eq!(node.branches[0].dev.len(), 720, "le prioritaire garde tout son pourtour");
        assert_eq!(node.branches[1].dev.len(), 720, "le chevauchant garde tout son pourtour");

        let alone = multi(&MultiInput { r1: 40.0, branches: vec![specs[1]], n_samples: 720 });
        let iso = &alone.branches[0].dev;
        let mut lifted = 0usize;
        for (p, q) in node.branches[1].dev.iter().zip(iso.iter()) {
            assert!((p.theta - q.theta).abs() < 1e-12);
            if (p.v - q.v).abs() < 1e-9 {
                continue; // free sector: landing untouched
            }
            // Lifted sector: strictly above the wall, exactly ON the back
            // of the through member.
            assert!(p.v > q.v, "la lèvre ne descend jamais sous l'atterrissage");
            let w = to_world(&specs[1], p);
            assert!(
                (dist_to_axis(&specs[0], &w) - specs[0].r).abs() < 1e-9,
                "selle hors du dos du prioritaire"
            );
            lifted += 1;
        }
        assert!(lifted > 30, "trop peu de génératrices en selle ({lifted})");
    }

    #[test]
    fn saddle_lip_rests_on_the_through_member_wall() {
        // Closure: every lifted lip point of branch 2 lies ON the surface
        // of branch 1 (distance to its axis = r) and outside the main tube
        // — the saddle rests on a wall that really exists.
        let specs = v_specs();
        let node = multi(&MultiInput { r1: 40.0, branches: specs.clone(), n_samples: 720 });
        let alone = multi(&MultiInput { r1: 40.0, branches: vec![specs[1]], n_samples: 720 });
        let mut checked = 0usize;
        for (p, q) in node.branches[1].dev.iter().zip(alone.branches[0].dev.iter()) {
            if (p.v - q.v).abs() < 1e-9 {
                continue;
            }
            let w = to_world(&specs[1], p);
            let dist_axis = dist_to_axis(&specs[0], &w);
            assert!(
                (dist_axis - specs[0].r).abs() < 1e-9,
                "selle hors de la paroi maîtresse : {dist_axis}"
            );
            let dist_main = (w.x * w.x + w.y * w.y).sqrt();
            assert!(dist_main >= 40.0 - 1e-6, "selle sous le tube principal");
            checked += 1;
        }
        assert!(checked > 30, "trop peu de points de selle ({checked})");
    }

    #[test]
    fn concurrent_axes_truss_node_lifts_the_overlapping_lips() {
        // Three tubes fanning onto the main tube, ALL axes through its
        // centre (z = 0, ψ = 0) — the truss-node layout.  Their footprints
        // overlap on the wall, so the lower-priority tubes MUST ride over
        // the higher ones: the lip lifts off the wall onto the back of the
        // through member.  This is the case where the first entry into the
        // neighbour's cylinder sits INSIDE the main tube — the obstacle
        // that counts is the EXIT point, probed for real material.
        let specs = vec![
            MultiBranchSpec { r: 30.0, z: 0.0, phi: 45f64.to_radians(), psi: 0.0 },
            MultiBranchSpec { r: 25.0, z: 0.0, phi: 90f64.to_radians(), psi: 0.0 },
            MultiBranchSpec { r: 22.5, z: 0.0, phi: 135f64.to_radians(), psi: 0.0 },
        ];
        let node = multi(&MultiInput { r1: 50.0, branches: specs.clone(), n_samples: 720 });

        // The through member keeps its full, untouched landing.
        assert_eq!(node.branches[0].dev.len(), 720, "le prioritaire reste entier");
        assert!(!node.branches[0].cut_by_neighbor);

        // 2 and 3 keep their full circumference but ride over the others.
        for bi in [1usize, 2] {
            let b = &node.branches[bi];
            assert_eq!(b.dev.len(), 720, "P{} garde tout son pourtour", bi + 1);
            assert!(b.cut_by_neighbor, "P{} devrait chevaucher un prioritaire", bi + 1);
            let alone =
                multi(&MultiInput { r1: 50.0, branches: vec![specs[bi]], n_samples: 720 });
            let mut lifted = 0usize;
            for (p, q) in b.dev.iter().zip(alone.branches[0].dev.iter()) {
                if (p.v - q.v).abs() < 1e-9 {
                    continue;
                }
                assert!(p.v > q.v, "la lèvre ne descend jamais sous l'atterrissage");
                // The lifted lip rests exactly on SOME higher-priority wall,
                // outside the main tube: the joint closes on real material.
                let w = to_world(&specs[bi], p);
                let on_prior = specs[..bi]
                    .iter()
                    .any(|t| (dist_to_axis(t, &w) - t.r).abs() < 1e-9);
                assert!(on_prior, "selle P{} hors de toute paroi prioritaire", bi + 1);
                let dist_main = (w.x * w.x + w.y * w.y).sqrt();
                assert!(dist_main >= 50.0 - 1e-6, "selle sous le tube principal");
                lifted += 1;
            }
            assert!(lifted > 50, "P{} : trop peu de selle ({lifted})", bi + 1);
        }
    }

    /// Closest point between two axes (lines), as a mid-point of the
    /// closest-approach segment — an INDEPENDENT check of the closed form.
    fn axes_crossing(a: &MultiBranchSpec, b: &MultiBranchSpec) -> nalgebra::Vector3<f64> {
        let dir = |s: &MultiBranchSpec| {
            (crate::geometry::rot_z(s.psi) * crate::geometry::rot_x(s.phi))
                * nalgebra::Vector3::z()
        };
        let (p1, d1) = (nalgebra::Vector3::new(0.0, 0.0, a.z), dir(a));
        let (p2, d2) = (nalgebra::Vector3::new(0.0, 0.0, b.z), dir(b));
        let r = p1 - p2;
        let (a11, b12, c22) = (d1.dot(&d1), d1.dot(&d2), d2.dot(&d2));
        let (d1r, d2r) = (d1.dot(&r), d2.dot(&r));
        let den = a11 * c22 - b12 * b12;
        let t = (b12 * d2r - c22 * d1r) / den;
        let u = (a11 * d2r - b12 * d1r) / den;
        ((p1 + t * d1) + (p2 + u * d2)) * 0.5
    }

    #[test]
    fn concurrent_axes_have_no_eccentricity() {
        // The application default: every axis through the centre, e = 0.
        let specs = vec![
            MultiBranchSpec { r: 30.0, z: 0.0, phi: 45f64.to_radians(), psi: 0.0 },
            MultiBranchSpec { r: 25.0, z: 0.0, phi: 135f64.to_radians(), psi: 0.0 },
        ];
        let node = multi(&MultiInput { r1: 50.0, branches: specs, n_samples: 720 });
        assert_eq!(node.pairs.len(), 1, "une paire coplanaire attendue");
        let p = node.pairs[0];
        assert!(p.same_side);
        assert!(p.eccentricity.unwrap().abs() < 1e-9, "axes concourants : e = 0");
    }

    #[test]
    fn eccentricity_matches_the_axes_crossing_point() {
        // K joint with axial offsets: the closed form must land on the point
        // where the two axes really cross, measured from the chord axis.
        let specs = vec![
            MultiBranchSpec { r: 30.0, z: -30.0, phi: 45f64.to_radians(), psi: 0.0 },
            MultiBranchSpec { r: 25.0, z: 30.0, phi: 135f64.to_radians(), psi: 0.0 },
        ];
        let node = multi(&MultiInput { r1: 50.0, branches: specs.clone(), n_samples: 720 });
        let e = node.pairs[0].eccentricity.unwrap();
        // Closed form: (z2 − z1)·sinφ1·sinφ2 / sin(φ2 − φ1) = 60·0.5/1 = 30.
        assert!((e - 30.0).abs() < 1e-9, "forme close inattendue : {e}");
        // Independent: distance from the chord axis to the crossing point,
        // signed along the outward radial direction of the joint plane.
        let w = axes_crossing(&specs[0], &specs[1]);
        let n = nalgebra::Vector3::new(specs[0].psi.sin(), -specs[0].psi.cos(), 0.0);
        assert!(
            (w - n * w.dot(&n)).norm() < 1e-9,
            "le point de croisement doit être dans le plan du nœud"
        );
        assert!((w.dot(&n) - e).abs() < 1e-9, "excentrement ≠ point de croisement");

        // Opposed branches (X joint) use sin(φ1 + φ2) — check the other sign.
        let opp = vec![
            MultiBranchSpec { r: 30.0, z: -20.0, phi: 60f64.to_radians(), psi: 0.0 },
            MultiBranchSpec {
                r: 25.0,
                z: 20.0,
                phi: 60f64.to_radians(),
                psi: std::f64::consts::PI,
            },
        ];
        let nx = multi(&MultiInput { r1: 50.0, branches: opp.clone(), n_samples: 720 });
        let ex = nx.pairs[0].eccentricity.unwrap();
        assert!(!nx.pairs[0].same_side, "piquages opposés");
        let wx = axes_crossing(&opp[0], &opp[1]);
        let nvec = nalgebra::Vector3::new(opp[0].psi.sin(), -opp[0].psi.cos(), 0.0);
        assert!((wx.dot(&nvec) - ex).abs() < 1e-9, "excentrement X ≠ croisement");
    }

    #[test]
    fn gap_and_overlap_are_read_in_the_joint_plane() {
        // A brace pierces the chord wall where its own axis exits it, so with
        // CONCURRENT axes two opposed braces necessarily land far apart: the
        // node has a gap, and eccentricity is what closes it.  Both regimes
        // must be measured in the plane of the joint.
        let mk = |z: f64, phi_deg: f64| MultiBranchSpec {
            r: 25.0,
            z,
            phi: phi_deg.to_radians(),
            psi: 0.0,
        };
        let circ = std::f64::consts::TAU * 50.0;
        let u_crown = 50.0 * -std::f64::consts::FRAC_PI_2;
        // Isolated footprint of one brace, measured on the crown generatrix.
        // Recomputed alone on purpose: as soon as two openings overlap the
        // payload only carries their merged envelope.
        let span = |spec: MultiBranchSpec| {
            let solo = multi(&MultiInput { r1: 50.0, branches: vec![spec], n_samples: 720 });
            super::span_at_u(&solo.holes[0].pts, u_crown, circ).unwrap()
        };

        // Concurrent axes (e = 0): a real gap between the two footprints.
        let concurrent = multi(&MultiInput {
            r1: 50.0,
            branches: vec![mk(0.0, 45.0), mk(0.0, 135.0)],
            n_samples: 720,
        });
        let p = concurrent.pairs[0];
        assert!(p.eccentricity.unwrap().abs() < 1e-9);
        assert!(p.overlap_pct.is_none(), "axes concourants : aucun recouvrement");
        let g = p.gap.expect("jeu attendu");
        let (s0, s1) = (span(mk(0.0, 45.0)), span(mk(0.0, 135.0)));
        // Footprint 1 sits BELOW footprint 0 here (each brace pierces on its
        // own side), so the gap is the distance between the facing lips.
        let expected = (s1.0 - s0.1).max(s0.0 - s1.1);
        assert!((g - expected).abs() < 1e-9, "jeu {g} ≠ écart des empreintes {expected}");
        assert!(g > 1.0, "un vrai jeu est attendu : {g}");

        // Offset the braces towards each other: the footprints overlap and
        // λov appears, read on the lower-priority (overlapping) brace.
        let offset = multi(&MultiInput {
            r1: 50.0,
            branches: vec![mk(-30.0, 45.0), mk(30.0, 135.0)],
            n_samples: 720,
        });
        let q = offset.pairs[0];
        assert!(q.gap.is_none(), "recouvrement : pas de jeu");
        let ov = q.overlap_pct.expect("λov attendu");
        let (t0, t1) = (span(mk(-30.0, 45.0)), span(mk(30.0, 135.0)));
        let qlen = t0.1.min(t1.1) - t0.0.max(t1.0);
        let plen = t1.1 - t1.0;
        assert!(
            (ov - 100.0 * qlen / plen).abs() < 1e-9,
            "λov {ov} ≠ q/p mesuré ({} / {})",
            qlen,
            plen
        );
        assert!(ov > 20.0 && ov < 100.0, "λov hors du plausible : {ov}");
        // And the eccentricity that produced it is the closed form.
        assert!((q.eccentricity.unwrap() - 30.0).abs() < 1e-9);
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
        let specs = vec![
            MultiBranchSpec { r: 30.0, z: -40.0, phi, psi: 0.0 },
            MultiBranchSpec { r: 30.0, z: 40.0, phi: std::f64::consts::PI - phi, psi: 0.0 },
        ];
        let node = multi(&MultiInput { r1: 50.0, branches: specs.clone(), n_samples: 720 });
        assert_eq!(node.holes.len(), 1, "les deux lumières fusionnent en enveloppe");
        let env = &node.holes[0].pts;
        let frames: Vec<BranchFrame> = specs.iter().map(BranchFrame::new).collect();
        for (i, br) in node.branches.iter().enumerate() {
            let mut checked = 0usize;
            for p in &br.curve3d {
                let dist_main = (p.x * p.x + p.y * p.y).sqrt();
                if (dist_main - 50.0).abs() > 1e-6 {
                    continue; // lifted lip riding a neighbour — not on the wall
                }
                // Landing points passing UNDER a neighbour are covered by
                // that neighbour's opening — the visible rim is the rest.
                let v = Vector3::new(p.x, p.y, p.z);
                let inside_neighbor = frames.iter().enumerate().any(|(j, fr)| {
                    if j == i {
                        return false;
                    }
                    let rel = v - fr.c;
                    (rel - fr.d * rel.dot(&fr.d)).norm() < fr.r - 1e-6
                });
                if inside_neighbor {
                    continue;
                }
                let u = 50.0 * p.y.atan2(p.x);
                let d = dist_to_polyline(u, p.z, env);
                assert!(d < 0.05, "rive visible à {d:.4} mm du contour d'enveloppe");
                checked += 1;
            }
            assert!(checked > 100, "trop peu de points de rive visibles ({checked})");
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
