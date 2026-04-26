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
    pub branch: Option<Branch>,
    /// Intersection curve in 3D space (mm).
    pub curve3d: Vec<Point3>,
    /// Developed pattern on cylinder 2 (the branch cut).  Unwrap length on
    /// the abscissa, axial coordinate on the ordinate, both in millimetres.
    pub dev_branch: Vec<DevPoint>,
    /// "Gueule de loup" pattern on cylinder 1.  Always present when both
    /// cylinders are involved; `None` for the cyl/plane mode.
    pub dev_main: Option<Vec<DevPoint>>,
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

    // "Gueule de loup": unwrap α on cylinder 1, then sort by `u1 = R1·α`.
    let dev_main = if alphas.is_empty() {
        None
    } else {
        let alpha_unwrapped = unwrap_angles(&alphas);
        let mut tmp: Vec<(f64, f64)> = alpha_unwrapped
            .into_iter()
            .zip(zs.iter().copied())
            .map(|(alpha, z)| (r1 * alpha, z))
            .collect();
        tmp.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Some(
            tmp.into_iter()
                .map(|(u, v)| DevPoint { theta: u / r1, u, v })
                .collect::<Vec<_>>(),
        )
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
        branch: Some(branch),
        curve3d,
        dev_branch,
        dev_main,
        bbox_branch,
        bbox_main,
        circumference_branch: Some(std::f64::consts::TAU * r2),
        circumference_main: std::f64::consts::TAU * r1,
    }
}

/// Cylinder–plane intersection: a single tube cut by an inclined plane.
///
/// The plane equation, in millimetres, is `z = z0 - y·tan(phi)` — equivalent
/// to a plane normal `(0, sin phi, cos phi)` passing through `(0, 0, z0)`.
pub fn cyl_plane(input: CylPlaneInput) -> IntersectionPayload {
    let CylPlaneInput { r1, phi, z0, n_samples } = input;
    let n = n_samples.max(64);
    let tan_phi = phi.tan();

    let mut curve3d = Vec::with_capacity(n);
    let mut dev_branch = Vec::with_capacity(n);

    for k in 0..n {
        let theta = std::f64::consts::TAU * (k as f64) / (n as f64);
        let (st, ct) = theta.sin_cos();
        let x = r1 * ct;
        let y = r1 * st;
        let z = z0 - y * tan_phi;
        curve3d.push(Point3 { x, y, z });
        dev_branch.push(DevPoint { theta, u: r1 * theta, v: z });
    }

    let bbox_branch = BBox2::from_points(dev_branch.iter().map(|p| (p.u, p.v)));

    IntersectionPayload {
        mode: "cyl_plane",
        r1,
        r2: None,
        phi,
        branch: None,
        curve3d,
        dev_branch,
        dev_main: None,
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
    fn plane_cut_returns_full_period() {
        let res = cyl_plane(CylPlaneInput {
            r1: 30.0,
            phi: std::f64::consts::FRAC_PI_4,
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
