//! Tauri entry point — exposes the cylinder-intersection routines as
//! `invoke`-able commands consumed by the Svelte frontend.

pub mod geometry;
pub mod intersection;

use intersection::{cyl_cyl, cyl_plane, CylCylInput, CylPlaneInput, IntersectionPayload};

#[tauri::command]
fn intersect_cyl_cyl(input: CylCylInput) -> Result<IntersectionPayload, String> {
    if input.r1 <= 0.0 || input.r2 <= 0.0 {
        return Err("Radii must be strictly positive.".into());
    }
    if !input.phi.is_finite() {
        return Err("phi must be finite.".into());
    }
    Ok(cyl_cyl(input))
}

#[tauri::command]
fn intersect_cyl_plane(input: CylPlaneInput) -> Result<IntersectionPayload, String> {
    if input.r1 <= 0.0 {
        return Err("r1 must be strictly positive.".into());
    }
    let half_pi = std::f64::consts::FRAC_PI_2;
    if input.phi.abs() >= half_pi - 1e-3 {
        return Err(
            "phi must lie in (-π/2, π/2): the plane cannot become parallel to the cylinder axis."
                .into(),
        );
    }
    Ok(cyl_plane(input))
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            intersect_cyl_cyl,
            intersect_cyl_plane,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
