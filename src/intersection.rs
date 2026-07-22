//! Cylinder–cylinder and cylinder–plane intersection algebra plus the matching
//! "developed" (unrolled) patterns used to lay out the cuts at scale 1:1.

use serde::{Deserialize, Serialize};

use crate::geometry::{cyl2_point, unwrap_angles, BBox2, Branch, Point3};

/// Inputs for a cylinder–cylinder intersection.
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct CylCylInput {
    /// Radius of the main vertical cylinder, in millimetres.
    pub r1: f64,
    /// Radius of the inclined branch cylinder, in millimetres.
    pub r2: f64,
    /// Angle between the two cylinder axes, in radians (rotation around `Ox`).
    pub phi: f64,
    /// Number of samples along `θ ∈ [0, 2π)`.
    #[serde(default = "default_samples")]
    pub n_samples: usize,
    /// Which root of the quadratic to follow.
    #[serde(default)]
    pub branch: Branch,
}

/// Inputs for a cylinder–plane intersection (mitered cut on a single tube).
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct CylPlaneInput {
    /// Radius of the cylinder being cut, in millimetres.
    pub r1: f64,
    /// Tilt of the plane around `Ox`, in radians (`0` is horizontal).
    pub phi: f64,
    /// Additional tilt of the plane around `Oy`, in radians — a fully
    /// oriented cutting plane `z = z0 − y·tan(phi) − x·tan(phi_y)`.
    #[serde(default)]
    pub phi_y: f64,
    /// Vertical offset of the plane along `Oz`, in millimetres.
    #[serde(default)]
    pub z0: f64,
    /// Number of samples along `θ ∈ [0, 2π)`.
    #[serde(default = "default_samples")]
    pub n_samples: usize,
}

fn default_samples() -> usize {
    1440
}

/// 2D developed point in millimetres, plus the originating `θ` angle.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct DevPoint {
    pub theta: f64,
    pub u: f64,
    pub v: f64,
}

/// Final response payload returned to the frontend.
#[derive(Debug, Clone, Serialize)]
pub struct IntersectionPayload {
    pub mode: &'static str,
    pub r1: f64,
    pub r2: Option<f64>,
    pub phi: f64,
    /// Second plane tilt (around `Oy`); `0` outside the oriented-plane mode.
    pub phi_y: f64,
    pub branch: Option<Branch>,
    /// Intersection curve in 3D space (mm).
    pub curve3d: Vec<Point3>,
    /// Developed pattern on cylinder 2 (the branch cut).  Unwrap length on
    /// the abscissa, axial coordinate on the ordinate, both in millimetres.
    pub dev_branch: Vec<DevPoint>,
    /// "Gueule de loup" pattern on cylinder 1.  Always present when both
    /// cylinders are involved; `None` for the cyl/plane mode.  Points keep
    /// the θ traversal order so the polyline follows the actual contour.
    pub dev_main: Option<Vec<DevPoint>>,
    /// Whether `dev_main` is a closed loop.  True for a full penetration
    /// (`r2 ≤ r1`, the opening is a closed "egg"); false when the curve
    /// wraps around the whole main cylinder (`r2 > r1`).
    pub dev_main_closed: bool,
    /// Bounding boxes (mm) for layout & SVG sizing.
    pub bbox_branch: Option<BBox2>,
    pub bbox_main: Option<BBox2>,
    /// Total branch unwrap length (`2π·r2`) and main unwrap length (`2π·r1`).
    pub circumference_branch: Option<f64>,
    pub circumference_main: f64,
}

/// Compute the intersection of two cylinders and the two associated developed
/// patterns (branch + "gueule de loup" on the main cylinder).
pub fn cyl_cyl(input: CylCylInput) -> IntersectionPayload {
    let CylCylInput { r1, r2, phi, n_samples, branch } = input;
    let n = n_samples.max(64);

    let (sp, cp) = phi.sin_cos();
    let a = sp * sp;

    let mut curve3d = Vec::with_capacity(n);
    let mut dev_branch = Vec::with_capacity(n);
    let mut alphas = Vec::with_capacity(n);
    let mut zs = Vec::with_capacity(n);

    for k in 0..n {
        let theta = std::f64::consts::TAU * (k as f64) / (n as f64);
        let (st, ct) = theta.sin_cos();
        let b = -2.0 * r2 * st * cp * sp;
        let c0 = r2 * r2 * (ct * ct + st * st * cp * cp) - r1 * r1;

        let t = if a.abs() < 1e-14 {
            // Degenerate case — coaxial cylinders. The intersection only
            // exists on the full circle and only when r1 == r2; we just skip.
            f64::NAN
        } else {
            let disc = b * b - 4.0 * a * c0;
            if disc < 0.0 {
                f64::NAN
            } else {
                let sq = disc.sqrt();
                let denom = 2.0 * a;
                match branch {
                    Branch::Outer => (-b + sq) / denom,
                    Branch::Inner => (-b - sq) / denom,
                }
            }
        };

        if t.is_finite() {
            let p = cyl2_point(r2, phi, theta, t);
            curve3d.push(Point3::from(p));
            dev_branch.push(DevPoint { theta, u: r2 * theta, v: t });

            // Stash for the main cylinder unwrap.
            alphas.push(p.y.atan2(p.x));
            zs.push(p.z);
        }
    }

    // "Gueule de loup": unwrap α on cylinder 1 and keep the θ traversal
    // order — the samples then follow the physical contour of the opening
    // (sorting by `u` would interleave the two lips into a zigzag).
    let (dev_main, dev_main_closed) = if alphas.is_empty() {
        (None, false)
    } else {
        let alpha_unwrapped = unwrap_angles(&alphas);
        let pts = alpha_unwrapped
            .into_iter()
            .zip(zs.iter().copied())
            .map(|(alpha, z)| DevPoint { theta: alpha, u: r1 * alpha, v: z })
            .collect::<Vec<_>>();
        let closed = polyline_is_loop(&pts_uv(&pts));
        (Some(pts), closed)
    };

    let bbox_branch = BBox2::from_points(dev_branch.iter().map(|p| (p.u, p.v)));
    let bbox_main = dev_main
        .as_ref()
        .and_then(|pts| BBox2::from_points(pts.iter().map(|p| (p.u, p.v))));

    IntersectionPayload {
        mode: "cyl_cyl",
        r1,
        r2: Some(r2),
        phi,
        phi_y: 0.0,
        branch: Some(branch),
        curve3d,
        dev_branch,
        dev_main,
        dev_main_closed,
        bbox_branch,
        bbox_main,
        circumference_branch: Some(std::f64::consts::TAU * r2),
        circumference_main: std::f64::consts::TAU * r1,
    }
}

fn pts_uv(pts: &[DevPoint]) -> Vec<(f64, f64)> {
    pts.iter().map(|p| (p.u, p.v)).collect()
}

/// Whether a developed polyline in traversal order forms a closed loop:
/// the gap between its endpoints must be comparable to the sampling step,
/// not to the size of the pattern.
fn polyline_is_loop(pts: &[(f64, f64)]) -> bool {
    if pts.len() < 8 {
        return false;
    }
    let (fu, fv) = pts[0];
    let (lu, lv) = pts[pts.len() - 1];
    let gap = ((lu - fu).powi(2) + (lv - fv).powi(2)).sqrt();
    let total: f64 = pts
        .windows(2)
        .map(|w| ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt())
        .sum();
    let mean_step = total / (pts.len() - 1) as f64;
    gap <= 5.0 * mean_step.max(1e-9)
}

/// Cylinder–plane intersection: a single tube cut by a fully oriented
/// plane `z = z0 − y·tan(phi) − x·tan(phi_y)` (normal proportional to
/// `(tan(phi_y), tan(phi), 1)`, passing through `(0, 0, z0)`).
///
/// The developed cut stays an exact sinusoid: with `a = tan(phi_y)` and
/// `b = tan(phi)`, `v(u) = z0 − R·√(a²+b²)·sin(u/R + ψ)`, `ψ = atan2(a, b)`.
pub fn cyl_plane(input: CylPlaneInput) -> IntersectionPayload {
    let CylPlaneInput { r1, phi, phi_y, z0, n_samples } = input;
    let n = n_samples.max(64);
    let tan_phi = phi.tan();
    let tan_phi_y = phi_y.tan();

    let mut curve3d = Vec::with_capacity(n);
    let mut dev_branch = Vec::with_capacity(n);

    for k in 0..n {
        let theta = std::f64::consts::TAU * (k as f64) / (n as f64);
        let (st, ct) = theta.sin_cos();
        let x = r1 * ct;
        let y = r1 * st;
        let z = z0 - y * tan_phi - x * tan_phi_y;
        curve3d.push(Point3 { x, y, z });
        dev_branch.push(DevPoint { theta, u: r1 * theta, v: z });
    }

    let bbox_branch = BBox2::from_points(dev_branch.iter().map(|p| (p.u, p.v)));

    IntersectionPayload {
        mode: "cyl_plane",
        r1,
        r2: None,
        phi,
        phi_y,
        branch: None,
        curve3d,
        dev_branch,
        dev_main: None,
        dev_main_closed: false,
        bbox_branch,
        bbox_main: None,
        circumference_branch: Some(std::f64::consts::TAU * r1),
        circumference_main: std::f64::consts::TAU * r1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perpendicular_equal_radii_produces_two_arches() {
        // The textbook case: two equal cylinders crossing at right angle.
        // In that situation the intersection is a pair of ellipses and the
        // outer branch should reach a maximum at θ = π/2.
        let res = cyl_cyl(CylCylInput {
            r1: 50.0,
            r2: 50.0,
            phi: std::f64::consts::FRAC_PI_2,
            n_samples: 720,
            branch: Branch::Outer,
        });
        assert!(!res.dev_branch.is_empty());
        let v_max = res.dev_branch.iter().map(|p| p.v).fold(f64::MIN, f64::max);
        assert!((v_max - 50.0).abs() < 1e-6, "v_max = {v_max}");
    }

    #[test]
    fn gueule_de_loup_is_a_closed_contour_for_full_penetration() {
        // r2 < r1: the opening is a closed loop and consecutive samples must
        // stay close to each other (no zigzag between the two lips).
        let res = cyl_cyl(CylCylInput {
            r1: 50.0,
            r2: 35.0,
            phi: std::f64::consts::FRAC_PI_4,
            n_samples: 1440,
            branch: Branch::Outer,
        });
        assert!(res.dev_main_closed);
        let pts = res.dev_main.unwrap();
        let (mut max_seg, mut span_u) = (0.0f64, 0.0f64);
        for w in pts.windows(2) {
            let d = ((w[1].u - w[0].u).powi(2) + (w[1].v - w[0].v).powi(2)).sqrt();
            max_seg = max_seg.max(d);
        }
        for p in &pts {
            span_u = span_u.max(p.u.abs());
        }
        // A contour in traversal order has tiny segments; the sorted zigzag
        // of the old implementation produced segments spanning the height.
        assert!(max_seg < 2.0, "max segment = {max_seg} mm");
        assert!(span_u > 10.0);
    }

    #[test]
    fn gueule_de_loup_wraps_open_when_branch_is_larger() {
        // r2 > r1: the intersection circles the main cylinder — not a loop
        // in the developed plane.
        let res = cyl_cyl(CylCylInput {
            r1: 30.0,
            r2: 45.0,
            phi: std::f64::consts::FRAC_PI_3,
            n_samples: 1440,
            branch: Branch::Outer,
        });
        assert!(!res.dev_main_closed);
    }

    #[test]
    fn oriented_plane_is_a_phase_shifted_sinusoid() {
        // Two tilts: v(u) = z0 − R·√(a²+b²)·sin(u/R + ψ) with a = tan(φy),
        // b = tan(φx).  Check amplitude and the phase via the maximum point.
        let (r1, phix, phiy) = (30.0, 0.5f64, 0.35f64);
        let res = cyl_plane(CylPlaneInput {
            r1,
            phi: phix,
            phi_y: phiy,
            z0: 10.0,
            n_samples: 4096,
        });
        let amp = r1 * (phix.tan().hypot(phiy.tan()));
        let v_max = res.dev_branch.iter().map(|p| p.v).fold(f64::MIN, f64::max);
        let v_min = res.dev_branch.iter().map(|p| p.v).fold(f64::MAX, f64::min);
        assert!((v_max - (10.0 + amp)).abs() < 1e-3, "v_max = {v_max}");
        assert!((v_min - (10.0 - amp)).abs() < 1e-3, "v_min = {v_min}");
        // Maximum at θ* = −ψ + 3π/2 (mod 2π), ψ = atan2(a, b).
        let psi = phiy.tan().atan2(phix.tan());
        let theta_star = (-psi + 1.5 * std::f64::consts::PI).rem_euclid(std::f64::consts::TAU);
        let best = res
            .dev_branch
            .iter()
            .max_by(|a, b| a.v.partial_cmp(&b.v).unwrap())
            .unwrap();
        let dtheta = (best.theta - theta_star).abs().min(std::f64::consts::TAU - (best.theta - theta_star).abs());
        assert!(dtheta < 0.01, "θ_max = {}, attendu {}", best.theta, theta_star);
    }

    #[test]
    fn plane_cut_returns_full_period() {
        let res = cyl_plane(CylPlaneInput {
            r1: 30.0,
            phi: std::f64::consts::FRAC_PI_4,
            phi_y: 0.0,
            z0: 0.0,
            n_samples: 360,
        });
        assert_eq!(res.dev_branch.len(), 360);
        let bbox = res.bbox_branch.unwrap();
        // The samples cover [0, 2π); at 360 samples the last one sits one step
        // before 2π, hence the spread approaches but does not reach 2π·R.
        let circ = std::f64::consts::TAU * 30.0;
        let span = bbox.u_max - bbox.u_min;
        assert!(span > 0.99 * circ * (359.0 / 360.0), "span = {span}, circ = {circ}");
    }
}
