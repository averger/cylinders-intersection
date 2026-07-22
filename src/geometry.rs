//! Core geometric primitives shared by the intersection routines.
//!
//! Conventions
//! -----------
//! * Cylinder 1 (the "main", or "vertical") has its axis on `Oz`, radius `r1`,
//!   and equation `x² + y² = r1²`.
//! * Cylinder 2 (the "branch", or "inclined") is obtained by rotating an
//!   axis-aligned cylinder around the `Ox` axis by `phi` radians.  Its
//!   parameterisation, before rotation, is `(r2·cosθ, r2·sinθ, t)`.
//! * A plane is described by an angle `phi` measured around the `Ox` axis
//!   relative to the horizontal plane `z = 0`, and a vertical offset `z0`.
//!
//! All lengths are expressed in millimetres so that the SVG export can be
//! rendered at print scale 1:1.

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Branch {
    /// `+sqrt(Δ)` root — the outer (entry side) curve.
    #[default]
    Outer,
    /// `-sqrt(Δ)` root — the inner (exit side) curve.
    Inner,
}

/// 3D point used for plotting and SVG export.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl From<Vector3<f64>> for Point3 {
    fn from(v: Vector3<f64>) -> Self {
        Point3 { x: v.x, y: v.y, z: v.z }
    }
}

/// Rotation matrix around the `Ox` axis by `phi` radians.
#[inline]
pub fn rot_x(phi: f64) -> Matrix3<f64> {
    let (s, c) = phi.sin_cos();
    Matrix3::new(
        1.0, 0.0, 0.0,
        0.0,   c,  -s,
        0.0,   s,   c,
    )
}

/// Forward parameterisation of cylinder 2 (radius `r2`, tilted by `phi`):
/// the point is the matrix product `Rx(φ) · P₀(θ, t)` of equation (1) in
/// `docs/THEORY.md` — the code mirrors the derivation literally.
#[inline]
pub fn cyl2_point(r2: f64, phi: f64, theta: f64, t: f64) -> Vector3<f64> {
    let (st, ct) = theta.sin_cos();
    let p0 = Vector3::new(r2 * ct, r2 * st, t);
    rot_x(phi) * p0
}

/// Standardised summary of the bounding box of a list of 2D points.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct BBox2 {
    pub u_min: f64,
    pub u_max: f64,
    pub v_min: f64,
    pub v_max: f64,
}

impl BBox2 {
    pub fn from_points<I: IntoIterator<Item = (f64, f64)>>(pts: I) -> Option<Self> {
        let mut iter = pts.into_iter();
        let (u, v) = iter.next()?;
        let (mut u_min, mut u_max, mut v_min, mut v_max) = (u, u, v, v);
        for (u, v) in iter {
            if u < u_min { u_min = u; }
            if u > u_max { u_max = u; }
            if v < v_min { v_min = v; }
            if v > v_max { v_max = v; }
        }
        Some(BBox2 { u_min, u_max, v_min, v_max })
    }

    pub fn width(&self) -> f64 { self.u_max - self.u_min }
    pub fn height(&self) -> f64 { self.v_max - self.v_min }
}

/// Numerically-robust `unwrap` of an angle sequence (analogous to `numpy.unwrap`).
pub fn unwrap_angles(alphas: &[f64]) -> Vec<f64> {
    if alphas.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(alphas.len());
    out.push(alphas[0]);
    let mut prev = alphas[0];
    let two_pi = std::f64::consts::TAU;
    let mut offset = 0.0;
    for &a in &alphas[1..] {
        let mut diff = a + offset - prev;
        while diff > std::f64::consts::PI {
            offset -= two_pi;
            diff -= two_pi;
        }
        while diff < -std::f64::consts::PI {
            offset += two_pi;
            diff += two_pi;
        }
        let v = a + offset;
        out.push(v);
        prev = v;
    }
    out
}
