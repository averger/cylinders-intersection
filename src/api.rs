//! HTTP layer: exposes the intersection routines as JSON endpoints and serves
//! the Svelte single-page application bundled into the binary.

use std::sync::Arc;

use axum::{
    extract::State,
    http::{header, StatusCode, Uri},
    response::{IntoResponse, Response},
    routing::{get, post},
    Json, Router,
};
use serde::Serialize;

use crate::assets::Assets;
use crate::export::{render_dxf, render_pdf, ExportDocument, ExportError};
use crate::intersection::{cyl_cyl, cyl_plane, CylCylInput, CylPlaneInput, IntersectionPayload};
use crate::multi::{multi, MultiInput, MultiPayload};

#[derive(Clone)]
pub struct AppState {
    pub started_at: Arc<std::time::Instant>,
}

impl AppState {
    pub fn new() -> Self {
        Self {
            started_at: Arc::new(std::time::Instant::now()),
        }
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Serialize)]
struct Health {
    name: &'static str,
    version: &'static str,
    uptime_ms: u128,
}

pub fn router() -> Router {
    let state = AppState::new();
    Router::new()
        .route("/api/health", get(health))
        .route("/api/intersect/cyl-cyl", post(api_cyl_cyl))
        .route("/api/intersect/cyl-plane", post(api_cyl_plane))
        .route("/api/intersect/multi", post(api_multi))
        .route("/api/export/pdf", post(api_export_pdf))
        .route("/api/export/dxf", post(api_export_dxf))
        .fallback(static_fallback)
        .with_state(state)
}

async fn health(State(state): State<AppState>) -> Json<Health> {
    Json(Health {
        name: env!("CARGO_PKG_NAME"),
        version: env!("CARGO_PKG_VERSION"),
        uptime_ms: state.started_at.elapsed().as_millis(),
    })
}

async fn api_cyl_cyl(Json(input): Json<CylCylInput>) -> Result<Json<IntersectionPayload>, ApiError> {
    if input.r1 <= 0.0 || input.r2 <= 0.0 {
        return Err(ApiError::bad_request("Radii must be strictly positive."));
    }
    if !input.phi.is_finite() {
        return Err(ApiError::bad_request("phi must be finite."));
    }
    Ok(Json(cyl_cyl(input)))
}

async fn api_cyl_plane(Json(input): Json<CylPlaneInput>) -> Result<Json<IntersectionPayload>, ApiError> {
    if input.r1 <= 0.0 {
        return Err(ApiError::bad_request("r1 must be strictly positive."));
    }
    let half_pi = std::f64::consts::FRAC_PI_2;
    if input.phi.abs() >= half_pi - 1e-3 || input.phi_y.abs() >= half_pi - 1e-3 {
        return Err(ApiError::bad_request(
            "phi and phi_y must lie in (-π/2, π/2): the plane cannot become parallel to the cylinder axis.",
        ));
    }
    Ok(Json(cyl_plane(input)))
}

async fn api_multi(Json(input): Json<MultiInput>) -> Result<Json<MultiPayload>, ApiError> {
    if input.r1 <= 0.0 {
        return Err(ApiError::bad_request("r1 must be strictly positive."));
    }
    if input.branches.is_empty() || input.branches.len() > 8 {
        return Err(ApiError::bad_request("The node accepts between 1 and 8 branches."));
    }
    for (i, b) in input.branches.iter().enumerate() {
        if b.r <= 0.0 || b.r > input.r1 {
            return Err(ApiError::bad_request(format!(
                "Branch {}: radius must lie in (0, r1].",
                i + 1
            )));
        }
        if !(b.phi > 1e-3 && b.phi < std::f64::consts::PI - 1e-3) {
            return Err(ApiError::bad_request(format!(
                "Branch {}: phi must lie strictly between 0 and π (no branch parallel to the main axis).",
                i + 1
            )));
        }
        if !b.z.is_finite() || !b.psi.is_finite() {
            return Err(ApiError::bad_request(format!("Branch {}: invalid z or psi.", i + 1)));
        }
    }
    Ok(Json(multi(&input)))
}

async fn api_export_pdf(Json(doc): Json<ExportDocument>) -> Result<Response, ApiError> {
    let bytes = render_pdf(&doc).map_err(ApiError::from)?;
    Ok((
        StatusCode::OK,
        [
            (header::CONTENT_TYPE, "application/pdf".to_string()),
            (
                header::CONTENT_DISPOSITION,
                "attachment; filename=\"gabarits.pdf\"".to_string(),
            ),
        ],
        bytes,
    )
        .into_response())
}

async fn api_export_dxf(Json(doc): Json<ExportDocument>) -> Result<Response, ApiError> {
    let dxf = render_dxf(&doc).map_err(ApiError::from)?;
    Ok((
        StatusCode::OK,
        [
            (header::CONTENT_TYPE, "application/dxf".to_string()),
            (
                header::CONTENT_DISPOSITION,
                "attachment; filename=\"gabarits.dxf\"".to_string(),
            ),
        ],
        dxf,
    )
        .into_response())
}

impl From<ExportError> for ApiError {
    fn from(e: ExportError) -> Self {
        ApiError::bad_request(e.to_string())
    }
}

#[derive(Debug)]
pub struct ApiError {
    status: StatusCode,
    message: String,
}

impl ApiError {
    fn bad_request(msg: impl Into<String>) -> Self {
        ApiError {
            status: StatusCode::BAD_REQUEST,
            message: msg.into(),
        }
    }
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        let body = serde_json::json!({ "error": self.message });
        (self.status, Json(body)).into_response()
    }
}

/// Serve embedded SPA assets, falling back to `index.html` for any unknown
/// path so that the Svelte router can take over.
async fn static_fallback(uri: Uri) -> Response {
    let path = uri.path().trim_start_matches('/');
    let candidate = if path.is_empty() { "index.html" } else { path };

    if let Some(content) = Assets::get(candidate) {
        let mime = mime_guess::from_path(candidate).first_or_octet_stream();
        return (
            StatusCode::OK,
            [(header::CONTENT_TYPE, mime.as_ref().to_string())],
            content.data.into_owned(),
        )
            .into_response();
    }

    if let Some(content) = Assets::get("index.html") {
        return (
            StatusCode::OK,
            [(header::CONTENT_TYPE, "text/html; charset=utf-8".to_string())],
            content.data.into_owned(),
        )
            .into_response();
    }

    (
        StatusCode::NOT_FOUND,
        [(header::CONTENT_TYPE, "text/plain; charset=utf-8".to_string())],
        b"Frontend assets are not bundled. Run `npm run build` inside ./web.".to_vec(),
    )
        .into_response()
}
