//! Binary STL export — the exact 3D assembly as a surface mesh, in mm,
//! ready for Fusion 360 / FreeCAD / any CAD that imports mesh bodies.
//!
//! The mesh is generated from the same closed-form cuts as the templates:
//! * main tube wall between `z_lo` and `z_hi`, minus every opening
//!   (per-column solid spans, the complement of the developed hole loops);
//! * each branch tube meshed generator by generator from its developed cut
//!   `t_cut(θ)` to its far end — mutual seams included;
//! * plane mode: the tube truncated at the exact sinusoidal cut.
//!
//! Surfaces are zero-thickness (cut templates have no wall thickness to
//! guess); Fusion imports them as mesh bodies.

use nalgebra::Vector3;

use super::{ExportDocument, ExportError, SourceSpec};
use crate::geometry::{rot_x, rot_z};
use crate::intersection::{cyl_cyl, cyl_plane, DevPoint};
use crate::multi::multi;

type Tri = [Vector3<f64>; 3];

/// Wall of a cylinder of radius `r` around `Oz`, spanning `z ∈ [z_lo, z_hi]`,
/// minus the openings given as developed loops `(u, v)` (u periodic 2πr).
fn main_wall(
    r: f64,
    z_lo: f64,
    z_hi: f64,
    holes: &[(Vec<(f64, f64)>, bool)],
    out: &mut Vec<Tri>,
) {
    let circ = std::f64::consts::TAU * r;
    let spans_at = |u0: f64| -> Vec<(f64, f64)> {
        let mut cuts: Vec<(f64, f64)> = Vec::new();
        for (loop_pts, _) in holes {
            'shift: for kk in -2i32..=2 {
                let u = u0 + f64::from(kk) * circ;
                let mut vs: Vec<f64> = Vec::new();
                let m = loop_pts.len();
                for idx in 0..m {
                    let a = loop_pts[idx];
                    let b = loop_pts[(idx + 1) % m];
                    if (a.0 - u) * (b.0 - u) < 0.0 {
                        let s = (u - a.0) / (b.0 - a.0);
                        vs.push(a.1 + s * (b.1 - a.1));
                    }
                }
                if vs.len() >= 2 {
                    let lo = vs.iter().cloned().fold(f64::INFINITY, f64::min);
                    let hi = vs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    cuts.push((lo, hi));
                    break 'shift;
                }
            }
        }
        cuts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut spans = Vec::new();
        let mut lo = z_lo;
        for &(h_lo, h_hi) in &cuts {
            if h_lo > lo {
                spans.push((lo, h_lo));
            }
            lo = lo.max(h_hi);
        }
        if lo < z_hi {
            spans.push((lo, z_hi));
        }
        spans
    };

    let n_a = 720usize;
    let point = |alpha: f64, z: f64| {
        Vector3::new(r * alpha.cos(), r * alpha.sin(), z)
    };
    let mut prev_alpha = 0.0f64;
    let mut prev_spans = spans_at(0.0);
    for i in 1..=n_a {
        let alpha = std::f64::consts::TAU * (i as f64) / (n_a as f64);
        let spans = spans_at(r * alpha);
        if spans.len() == prev_spans.len() {
            for (sa, sb) in prev_spans.iter().zip(spans.iter()) {
                quad(
                    point(prev_alpha, sa.0),
                    point(prev_alpha, sa.1),
                    point(alpha, sb.0),
                    point(alpha, sb.1),
                    out,
                );
            }
        } else {
            // An opening starts/ends inside this hair-thin column.
            let dominant = if spans.len() > prev_spans.len() { &spans } else { &prev_spans };
            for s in dominant {
                quad(
                    point(prev_alpha, s.0),
                    point(prev_alpha, s.1),
                    point(alpha, s.0),
                    point(alpha, s.1),
                    out,
                );
            }
        }
        prev_alpha = alpha;
        prev_spans = spans;
    }
}

/// Local frame of a meshed tube: centre, axis and radial basis.
struct TubeFrame {
    c: Vector3<f64>,
    d: Vector3<f64>,
    u: Vector3<f64>,
    w: Vector3<f64>,
    r: f64,
}

/// Tube along `frame`, generators spanning `t ∈ [dev.v, t_end]` (or
/// `[t_end, dev.v]` when `t_end` is below), minus the crossing contours
/// given as closed developed loops (neighbour branches passing through).
fn branch_wall(
    frame: &TubeFrame,
    dev: &[DevPoint],
    holes: &[(Vec<(f64, f64)>, bool)],
    t_end: f64,
    out: &mut Vec<Tri>,
) {
    let TubeFrame { c, d, u, w, r } = *frame;
    if dev.len() < 3 {
        return;
    }
    let point = |theta: f64, t: f64| {
        let (st, ct) = theta.sin_cos();
        c + r * (ct * u + st * w) + t * d
    };
    let circ = std::f64::consts::TAU * r;
    let below = t_end < dev[0].v;

    // Kept spans of one generatrix: [t_wall, t_end] minus hole intervals.
    let spans_at = |k: usize| -> Vec<(f64, f64)> {
        let (lo0, hi0) = if below { (t_end, dev[k].v) } else { (dev[k].v, t_end) };
        let mut cuts: Vec<(f64, f64)> = Vec::new();
        let mut cap = f64::INFINITY;
        for (loop_pts, closed_flag) in holes {
            'shift: for kk in -2i32..=2 {
                let uu = dev[k].u + f64::from(kk) * circ;
                let mut vs: Vec<f64> = Vec::new();
                let m = loop_pts.len();
                for idx in 0..m {
                    let a = loop_pts[idx];
                    let b = loop_pts[(idx + 1) % m];
                    if (a.0 - uu) * (b.0 - uu) < 0.0 {
                        let s = (uu - a.0) / (b.0 - a.0);
                        vs.push(a.1 + s * (b.1 - a.1));
                    }
                }
                if !vs.is_empty() {
                    if *closed_flag && vs.len() >= 2 {
                        let lo = vs.iter().cloned().fold(f64::INFINITY, f64::min);
                        let hi = vs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                        cuts.push((lo, hi));
                    } else if !*closed_flag {
                        // Open saddle arc: the tube ENDS there (cap).
                        let v = vs.iter().cloned().fold(f64::INFINITY, f64::min);
                        cap = cap.min(v);
                    }
                    break 'shift;
                }
            }
        }
        let hi_eff = hi0.min(cap);
        cuts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut spans = Vec::new();
        let mut lo = lo0;
        for &(h_lo, h_hi) in &cuts {
            if h_lo > lo {
                spans.push((lo, h_lo.min(hi_eff)));
            }
            lo = lo.max(h_hi);
        }
        if lo < hi_eff {
            spans.push((lo, hi_eff));
        }
        spans.retain(|(a, b)| b > a);
        spans
    };

    let mean_step = std::f64::consts::TAU / dev.len() as f64;
    let mut emit = |ka: usize, theta_a: f64, kb: usize, theta_b: f64| {
        let sa = spans_at(ka);
        let sb = spans_at(kb);
        if sa.len() == sb.len() {
            for (a, b) in sa.iter().zip(sb.iter()) {
                quad(
                    point(theta_a, a.0),
                    point(theta_a, a.1),
                    point(theta_b, b.0),
                    point(theta_b, b.1),
                    out,
                );
            }
        } else {
            let dominant = if sa.len() > sb.len() { &sa } else { &sb };
            for s in dominant {
                quad(
                    point(theta_a, s.0),
                    point(theta_a, s.1),
                    point(theta_b, s.0),
                    point(theta_b, s.1),
                    out,
                );
            }
        }
    };
    for k in 0..dev.len() - 1 {
        if (dev[k + 1].theta - dev[k].theta).abs() < 4.0 * mean_step {
            emit(k, dev[k].theta, k + 1, dev[k + 1].theta);
        }
    }
    let last = dev.len() - 1;
    if dev[0].theta + std::f64::consts::TAU - dev[last].theta < 4.0 * mean_step {
        emit(last, dev[last].theta, 0, dev[0].theta + std::f64::consts::TAU);
    }
}

/// Two triangles for the quad (a0→a1) × (b0→b1).
fn quad(a0: Vector3<f64>, a1: Vector3<f64>, b0: Vector3<f64>, b1: Vector3<f64>, out: &mut Vec<Tri>) {
    out.push([a0, b0, a1]);
    out.push([b0, b1, a1]);
}

/// Build the whole assembly for any source mode.
fn build_triangles(doc: &ExportDocument) -> Result<Vec<Tri>, ExportError> {
    let mut tris: Vec<Tri> = Vec::new();
    match &doc.source {
        SourceSpec::Multi(input) => {
            let payload = multi(input);
            if payload.holes.is_empty() {
                return Err(ExportError::EmptyPattern);
            }
            let (mut v_lo, mut v_hi) = (0.0f64, 0.0f64);
            for h in &payload.holes {
                if let Some(b) = h.bbox {
                    v_lo = v_lo.min(b.v_min);
                    v_hi = v_hi.max(b.v_max);
                }
            }
            let pad = (1.6 * payload.r1).max(70.0);
            let loops: Vec<(Vec<(f64, f64)>, bool)> = payload
                .holes
                .iter()
                .map(|h| (h.pts.iter().map(|p| (p.u, p.v)).collect(), h.closed))
                .collect();
            main_wall(payload.r1, v_lo - pad, v_hi + pad, &loops, &mut tris);

            for br in &payload.branches {
                let m = rot_z(br.psi) * rot_x(br.phi);
                let t_hi = br
                    .dev
                    .iter()
                    .map(|p| p.v)
                    .chain(br.holes.iter().flat_map(|h| h.pts.iter().map(|p| p.v)))
                    .fold(f64::NEG_INFINITY, f64::max);
                let branch_holes: Vec<(Vec<(f64, f64)>, bool)> = br
                    .holes
                    .iter()
                    .map(|h| (h.pts.iter().map(|p| (p.u, p.v)).collect(), h.closed))
                    .collect();
                branch_wall(
                    &TubeFrame {
                        c: Vector3::new(0.0, 0.0, br.z),
                        d: m * Vector3::z(),
                        u: m * Vector3::x(),
                        w: m * Vector3::y(),
                        r: br.r,
                    },
                    &br.dev,
                    &branch_holes,
                    t_hi + (3.0 * br.r).max(90.0),
                    &mut tris,
                );
            }
        }
        SourceSpec::CylCyl(input) => {
            let payload = cyl_cyl(*input);
            if payload.dev_branch.is_empty() {
                return Err(ExportError::EmptyPattern);
            }
            let r1 = payload.r1;
            let r2 = payload.r2.unwrap_or(r1);
            let max_abs_z = payload
                .curve3d
                .iter()
                .map(|p| p.z.abs())
                .fold(0.0f64, f64::max);
            let h = (4.0 * r1).max(2.4 * max_abs_z + 2.0 * r2).max(200.0) / 2.0;
            let hole: Vec<(Vec<(f64, f64)>, bool)> = match (&payload.dev_main, payload.dev_main_closed) {
                (Some(pts), true) => {
                    vec![(pts.iter().map(|p| (p.u, p.v)).collect(), true)]
                }
                _ => Vec::new(),
            };
            main_wall(r1, -h, h, &hole, &mut tris);

            let outer = !matches!(input.branch, crate::geometry::Branch::Inner);
            let m = rot_x(payload.phi);
            let (t_lo, t_hi) = payload.dev_branch.iter().fold(
                (f64::INFINITY, f64::NEG_INFINITY),
                |(lo, hi), p| (lo.min(p.v), hi.max(p.v)),
            );
            let ext = (3.0 * r2).max(90.0);
            let t_end = if outer { t_hi + ext } else { t_lo - ext };
            branch_wall(
                &TubeFrame {
                    c: Vector3::zeros(),
                    d: m * Vector3::z(),
                    u: m * Vector3::x(),
                    w: m * Vector3::y(),
                    r: r2,
                },
                &payload.dev_branch,
                &[],
                t_end,
                &mut tris,
            );
        }
        SourceSpec::CylPlane(input) => {
            let payload = cyl_plane(*input);
            if payload.dev_branch.is_empty() {
                return Err(ExportError::EmptyPattern);
            }
            let r1 = payload.r1;
            let amp = payload.phi.tan().hypot(payload.phi_y.tan()) * r1;
            let bottom = input.z0 - amp - (2.4 * r1).max(140.0);
            // The cut profile z(θ) is the developed curve itself.
            branch_wall(
                &TubeFrame {
                    c: Vector3::zeros(),
                    d: Vector3::z(),
                    u: Vector3::x(),
                    w: Vector3::y(),
                    r: r1,
                },
                &payload.dev_branch,
                &[],
                bottom,
                &mut tris,
            );
        }
    }
    if tris.is_empty() {
        return Err(ExportError::EmptyPattern);
    }
    Ok(tris)
}

/// Render the document's geometry as a binary STL byte stream (mm).
pub fn render_stl(doc: &ExportDocument) -> Result<Vec<u8>, ExportError> {
    let tris = build_triangles(doc)?;

    let mut out = Vec::with_capacity(84 + 50 * tris.len());
    let mut header = [0u8; 80];
    let tag = b"Cylix - exact cylinder intersections - units: mm";
    header[..tag.len()].copy_from_slice(tag);
    out.extend_from_slice(&header);
    out.extend_from_slice(&u32::try_from(tris.len()).unwrap_or(u32::MAX).to_le_bytes());

    for t in &tris {
        let n = (t[1] - t[0]).cross(&(t[2] - t[0]));
        let n = if n.norm() > 1e-12 { n.normalize() } else { Vector3::zeros() };
        for v in [n, t[0], t[1], t[2]] {
            for c in [v.x, v.y, v.z] {
                out.extend_from_slice(&(c as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&0u16.to_le_bytes());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::super::{
        Layers, Orientation, PageFormat, PageSpec, ScaleMode, TitleBlock,
    };
    use super::*;
    use crate::geometry::Branch;
    use crate::intersection::CylCylInput;
    use crate::multi::{MultiBranchSpec, MultiInput};

    fn doc_with(source: SourceSpec) -> ExportDocument {
        ExportDocument {
            source,
            page: PageSpec {
                format: PageFormat::A4,
                orientation: Orientation::Landscape,
                margin_mm: 10.0,
            },
            scale: ScaleMode::OneToOne,
            layers: Layers::default(),
            title_block: TitleBlock::default(),
            annotations: vec![],
            patterns: vec![super::super::PatternKind::Branch],
            cut_width_mm: 0.35,
        }
    }

    fn parse(stl: &[u8]) -> (u32, Vec<[f32; 3]>) {
        let count = u32::from_le_bytes(stl[80..84].try_into().unwrap());
        let mut verts = Vec::new();
        for i in 0..count as usize {
            let base = 84 + 50 * i;
            for v in 0..3 {
                let off = base + 12 + 12 * v;
                let x = f32::from_le_bytes(stl[off..off + 4].try_into().unwrap());
                let y = f32::from_le_bytes(stl[off + 4..off + 8].try_into().unwrap());
                let z = f32::from_le_bytes(stl[off + 8..off + 12].try_into().unwrap());
                verts.push([x, y, z]);
            }
        }
        (count, verts)
    }

    #[test]
    fn stl_is_structurally_sound_and_on_surface() {
        let doc = doc_with(SourceSpec::CylCyl(CylCylInput {
            r1: 50.0,
            r2: 35.0,
            phi: 1.0,
            n_samples: 720,
            branch: Branch::Outer,
        }));
        let stl = render_stl(&doc).unwrap();
        let (count, verts) = parse(&stl);
        assert_eq!(stl.len(), 84 + 50 * count as usize);
        assert!(count > 2000, "count = {count}");
        // Every vertex is finite and lies on one of the two cylinders.
        let m = crate::geometry::rot_x(1.0);
        let d2 = m * nalgebra::Vector3::z();
        let mut on1 = 0usize;
        for v in &verts {
            assert!(v.iter().all(|c| c.is_finite()));
            let p = nalgebra::Vector3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2]));
            let dist1 = (p.x * p.x + p.y * p.y).sqrt();
            let dist2 = (p - d2 * p.dot(&d2)).norm();
            let ok1 = (dist1 - 50.0).abs() < 1e-3;
            let ok2 = (dist2 - 35.0).abs() < 1e-3;
            assert!(ok1 || ok2, "sommet hors surfaces : {dist1} / {dist2}");
            if ok1 {
                on1 += 1;
            }
        }
        assert!(on1 > 0 && on1 < verts.len(), "les deux tubes sont maillés");
    }

    #[test]
    fn stl_covers_all_three_modes() {
        let multi_doc = doc_with(SourceSpec::Multi(MultiInput {
            r1: 40.0,
            branches: vec![
                MultiBranchSpec { r: 25.0, z: -45.0, phi: std::f64::consts::FRAC_PI_4, psi: 0.0 },
                MultiBranchSpec {
                    r: 25.0,
                    z: 45.0,
                    phi: std::f64::consts::PI - std::f64::consts::FRAC_PI_4,
                    psi: 0.0,
                },
            ],
            n_samples: 720,
        }));
        let plane_doc = doc_with(SourceSpec::CylPlane(crate::intersection::CylPlaneInput {
            r1: 40.0,
            phi: 0.6,
            phi_y: 0.2,
            z0: 10.0,
            n_samples: 720,
        }));
        for d in [&multi_doc, &plane_doc] {
            let stl = render_stl(d).unwrap();
            let (count, _) = parse(&stl);
            assert!(count > 1000, "count = {count}");
        }
    }
}
