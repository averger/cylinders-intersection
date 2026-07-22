//! Export pipeline: a single, backend-owned document model rendered to
//! print-accurate PDF (vector, mm-exact) and CAD-ready DXF (R12).
//!
//! The frontend editor manipulates exactly this JSON model; the backend
//! recomputes the geometry from the source parameters (single source of
//! truth) and lays the patterns out according to the page specification.

mod dxf;
mod pdf;

pub use dxf::render_dxf;
pub use pdf::render_pdf;

use serde::Deserialize;

use crate::intersection::{cyl_cyl, cyl_plane, CylCylInput, CylPlaneInput, IntersectionPayload};

/// Which geometry feeds the document.  Mirrors the two compute endpoints.
#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case")]
pub enum SourceSpec {
    CylCyl(CylCylInput),
    CylPlane(CylPlaneInput),
}

impl SourceSpec {
    pub fn compute(&self) -> IntersectionPayload {
        match *self {
            SourceSpec::CylCyl(input) => cyl_cyl(input),
            SourceSpec::CylPlane(input) => cyl_plane(input),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PageFormat {
    A4,
    A3,
    A2,
}

impl PageFormat {
    /// Portrait dimensions in mm.
    pub fn dims_mm(self) -> (f64, f64) {
        match self {
            PageFormat::A4 => (210.0, 297.0),
            PageFormat::A3 => (297.0, 420.0),
            PageFormat::A2 => (420.0, 594.0),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Orientation {
    Portrait,
    Landscape,
}

#[derive(Debug, Clone, Copy, Deserialize)]
pub struct PageSpec {
    pub format: PageFormat,
    pub orientation: Orientation,
    /// Printable margin on all four sides, mm.
    #[serde(default = "default_margin")]
    pub margin_mm: f64,
}

fn default_margin() -> f64 {
    10.0
}

impl PageSpec {
    /// Page size in mm, orientation applied.
    pub fn size_mm(&self) -> (f64, f64) {
        let (w, h) = self.format.dims_mm();
        match self.orientation {
            Orientation::Portrait => (w, h),
            Orientation::Landscape => (h, w),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ScaleMode {
    /// Fit the pattern on a single page; the effective scale is printed on
    /// the sheet.  Not intended for direct tracing.
    Fit,
    /// True 1:1 scale, tiled over as many pages as needed with assembly
    /// marks and a glue overlap.
    OneToOne,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(default)]
pub struct Layers {
    pub grid: bool,
    pub frame: bool,
    pub axis: bool,
    pub labels: bool,
    pub scale_bar: bool,
}

impl Default for Layers {
    fn default() -> Self {
        Layers { grid: true, frame: true, axis: true, labels: true, scale_bar: true }
    }
}

#[derive(Debug, Clone, Deserialize, Default)]
#[serde(default)]
pub struct TitleBlock {
    pub title: String,
    pub project: String,
    pub author: String,
    pub date: String,
    pub notes: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PatternKind {
    /// Cut template of the inclined tube (or of the single tube in the
    /// cylinder/plane mode).
    Branch,
    /// "Gueule de loup" opening on the main cylinder.
    Main,
}

/// Free text annotation anchored in pattern coordinates (mm).
#[derive(Debug, Clone, Deserialize)]
pub struct Annotation {
    pub pattern: PatternKind,
    pub u: f64,
    pub v: f64,
    pub text: String,
    #[serde(default = "default_annotation_size")]
    pub size_mm: f64,
}

fn default_annotation_size() -> f64 {
    4.0
}

/// The complete export document — what the SVG editor produces.
#[derive(Debug, Clone, Deserialize)]
pub struct ExportDocument {
    pub source: SourceSpec,
    pub page: PageSpec,
    pub scale: ScaleMode,
    #[serde(default)]
    pub layers: Layers,
    #[serde(default)]
    pub title_block: TitleBlock,
    #[serde(default)]
    pub annotations: Vec<Annotation>,
    /// Which patterns to export, in order.
    pub patterns: Vec<PatternKind>,
    /// Cut line width, mm.
    #[serde(default = "default_cut_width")]
    pub cut_width_mm: f64,
}

fn default_cut_width() -> f64 {
    0.35
}

/// A pattern fully laid out in its own (u, v) mm coordinate system, ready
/// for a renderer.  `v` grows upwards.
#[derive(Debug, Clone)]
pub struct Sheet {
    pub kind: PatternKind,
    pub name: String,
    pub meta: String,
    /// The cut polyline, in mm.
    pub cut: Vec<(f64, f64)>,
    pub closed: bool,
    /// Unwrapped cylinder footprint: (u0, v0, width, height), mm.
    pub frame: (f64, f64, f64, f64),
    /// `θ = 0` generator abscissa, mm.
    pub axis_u: f64,
    /// Unwrapped circumference of the tube this pattern wraps onto, mm.
    pub circumference: f64,
    /// Overall bounding box of everything drawable: (u_min, v_min, u_max, v_max).
    pub bbox: (f64, f64, f64, f64),
    /// Bounding box of what must actually reach the paper at 1:1 — the cut
    /// curve and the annotations.  The 1:1 tiling covers THIS box, so no
    /// page is wasted on the empty part of the unwrapped frame.
    pub cut_bbox: (f64, f64, f64, f64),
    pub annotations: Vec<Annotation>,
}

impl Sheet {
    /// Whether the rectangle `(u0..u0+w, v0..v0+h)`, grown by `margin`,
    /// contains any part of the cut curve or an annotation anchor.
    pub fn rect_has_content(&self, u0: f64, v0: f64, w: f64, h: f64, margin: f64) -> bool {
        let (lo_u, hi_u) = (u0 - margin, u0 + w + margin);
        let (lo_v, hi_v) = (v0 - margin, v0 + h + margin);
        let inside = |u: f64, v: f64| u >= lo_u && u <= hi_u && v >= lo_v && v <= hi_v;
        self.cut.iter().any(|&(u, v)| inside(u, v))
            || self.annotations.iter().any(|a| inside(a.u, a.v))
    }

    /// True tube generatrices (multiples of a quarter turn, `u = 0 ⇔ θ = 0`)
    /// falling within the drawing extents, with their angle in degrees.
    pub fn generatrices(&self) -> Vec<(f64, u32)> {
        let (u0, _, u1, _) = self.bbox;
        let step = self.circumference / 4.0;
        if step <= 0.0 {
            return Vec::new();
        }
        let k0 = (u0 / step).floor() as i64;
        let mut out = Vec::new();
        let mut k = k0;
        while (k as f64) * step <= u1 + 1e-9 {
            let u = (k as f64) * step;
            if u >= u0 - 1e-9 {
                out.push((u, (((k * 90) % 360 + 360) % 360) as u32));
            }
            k += 1;
        }
        out
    }

    /// `v` values where the cut polyline crosses the vertical line `u = g` —
    /// the alignment tick positions on a generatrix.
    pub fn curve_crossings(&self, g: f64) -> Vec<f64> {
        let mut out = Vec::new();
        let n = self.cut.len();
        if n < 2 {
            return out;
        }
        let last = if self.closed { n } else { n - 1 };
        for i in 0..last {
            let (ua, va) = self.cut[i];
            let (ub, vb) = self.cut[(i + 1) % n];
            if (ua - g) * (ub - g) < 0.0 {
                let s = (g - ua) / (ub - ua);
                out.push(va + s * (vb - va));
            }
        }
        out
    }
}

/// Errors surfaced to the API layer.
#[derive(Debug, thiserror::Error)]
pub enum ExportError {
    #[error("le calcul ne produit aucune courbe : paramètres hors domaine")]
    EmptyPattern,
    #[error("motif demandé indisponible pour ce mode ({0:?})")]
    MissingPattern(PatternKind),
    #[error("{0}")]
    Invalid(String),
}

/// Compute the geometry and lay out every requested pattern.
pub fn build_sheets(doc: &ExportDocument) -> Result<Vec<Sheet>, ExportError> {
    if doc.patterns.is_empty() {
        return Err(ExportError::Invalid("aucun motif sélectionné".into()));
    }
    if !(doc.cut_width_mm > 0.0 && doc.cut_width_mm <= 5.0) {
        return Err(ExportError::Invalid("largeur de trait hors limites".into()));
    }
    let payload = doc.source.compute();
    if payload.dev_branch.is_empty() {
        return Err(ExportError::EmptyPattern);
    }

    let mut sheets = Vec::with_capacity(doc.patterns.len());
    for &kind in &doc.patterns {
        sheets.push(layout_sheet(doc, &payload, kind)?);
    }
    Ok(sheets)
}

fn layout_sheet(
    doc: &ExportDocument,
    payload: &IntersectionPayload,
    kind: PatternKind,
) -> Result<Sheet, ExportError> {
    // The branch development is an OPEN curve: u = 0 and u = 2πR coincide
    // once the sheet is rolled, so no closing chord must ever be drawn.
    let (points, closed, circumference, diameter, name) = match kind {
        PatternKind::Branch => (
            &payload.dev_branch,
            false,
            payload
                .circumference_branch
                .unwrap_or(payload.circumference_main),
            payload.r2.unwrap_or(payload.r1) * 2.0,
            if payload.mode == "cyl_cyl" {
                "Gabarit tube incliné"
            } else {
                "Gabarit coupe en sifflet"
            },
        ),
        PatternKind::Main => {
            let pts = payload
                .dev_main
                .as_ref()
                .ok_or(ExportError::MissingPattern(kind))?;
            (
                pts,
                payload.dev_main_closed,
                payload.circumference_main,
                payload.r1 * 2.0,
                "Gueule de loup — cylindre principal",
            )
        }
    };
    if points.is_empty() {
        return Err(ExportError::EmptyPattern);
    }

    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    let cut: Vec<(f64, f64)> = points.iter().map(|p| (p.u, p.v)).collect();
    for &(u, v) in &cut {
        u_min = u_min.min(u);
        u_max = u_max.max(u);
        v_min = v_min.min(v);
        v_max = v_max.max(v);
    }

    // The frame is the full unwrapped footprint of the tube — circumference
    // wide, aligned to the unwrap period containing the curve (after the
    // angle unwrap, `u` may live in any 2πR-period).
    let center = (u_min + u_max) / 2.0;
    let frame_start = (center / circumference).floor() * circumference;
    let frame = (frame_start, v_min, circumference, (v_max - v_min).max(0.1));

    let bbox = (
        u_min.min(frame_start),
        v_min,
        u_max.max(frame_start + circumference),
        v_max,
    );

    let annotations: Vec<Annotation> = doc
        .annotations
        .iter()
        .filter(|a| a.pattern == kind)
        .cloned()
        .collect();

    // What must land on paper: the curve itself plus annotation anchors.
    let mut cut_bbox = (u_min, v_min, u_max, v_max);
    for a in &annotations {
        cut_bbox.0 = cut_bbox.0.min(a.u);
        cut_bbox.1 = cut_bbox.1.min(a.v);
        cut_bbox.2 = cut_bbox.2.max(a.u);
        cut_bbox.3 = cut_bbox.3.max(a.v);
    }

    let phi_deg = payload.phi.to_degrees();
    let meta = match kind {
        PatternKind::Branch => {
            let mut m = format!(
                "Ø {diameter:.1} mm — périmètre {circumference:.1} mm — phi = {phi_deg:.1}°"
            );
            if payload.phi_y.abs() > 1e-9 {
                m.push_str(&format!(" — phi_y = {:.1}°", payload.phi_y.to_degrees()));
            }
            m
        }
        PatternKind::Main => format!(
            "Ø {diameter:.1} mm — lumière de piquage — phi = {phi_deg:.1}°"
        ),
    };

    Ok(Sheet {
        kind,
        name: name.to_string(),
        meta,
        cut,
        closed,
        frame,
        axis_u: 0.0,
        circumference,
        bbox,
        cut_bbox,
        annotations,
    })
}

/// Tiling plan for 1:1 output: which (col, row) tiles cover the sheet.
#[derive(Debug, Clone, Copy)]
pub struct TilePlan {
    pub cols: usize,
    pub rows: usize,
    /// Drawing step between consecutive tiles (printable size minus overlap), mm.
    pub step_x: f64,
    pub step_y: f64,
    /// Printable area of one page, mm.
    pub view_w: f64,
    pub view_h: f64,
    /// Glue overlap between adjacent tiles, mm.
    pub overlap: f64,
}

pub const TILE_OVERLAP_MM: f64 = 12.0;
/// Vertical strip reserved for the title block at the bottom of each page, mm.
pub const TITLE_BLOCK_H_MM: f64 = 24.0;

/// Compute the tile grid needed to cover `(w × h)` mm of drawing at 1:1.
pub fn plan_tiles(page: &PageSpec, w: f64, h: f64) -> TilePlan {
    let (pw, ph) = page.size_mm();
    let view_w = (pw - 2.0 * page.margin_mm).max(20.0);
    let view_h = (ph - 2.0 * page.margin_mm - TITLE_BLOCK_H_MM).max(20.0);
    let step_x = (view_w - TILE_OVERLAP_MM).max(10.0);
    let step_y = (view_h - TILE_OVERLAP_MM).max(10.0);
    let cols = if w <= view_w { 1 } else { 1 + ((w - view_w) / step_x).ceil() as usize };
    let rows = if h <= view_h { 1 } else { 1 + ((h - view_h) / step_y).ceil() as usize };
    TilePlan { cols, rows, step_x, step_y, view_w, view_h, overlap: TILE_OVERLAP_MM }
}

/// Human tile label: columns as letters, rows as numbers ("B3").
pub fn tile_label(col: usize, row: usize) -> String {
    let mut letters = String::new();
    let mut c = col;
    loop {
        letters.insert(0, (b'A' + (c % 26) as u8) as char);
        if c < 26 {
            break;
        }
        c = c / 26 - 1;
    }
    format!("{letters}{}", row + 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Branch;

    fn doc() -> ExportDocument {
        ExportDocument {
            source: SourceSpec::CylCyl(CylCylInput {
                r1: 50.0,
                r2: 35.0,
                phi: 1.0,
                n_samples: 720,
                branch: Branch::Outer,
            }),
            page: PageSpec {
                format: PageFormat::A4,
                orientation: Orientation::Landscape,
                margin_mm: 10.0,
            },
            scale: ScaleMode::OneToOne,
            layers: Layers::default(),
            title_block: TitleBlock::default(),
            annotations: vec![],
            patterns: vec![PatternKind::Branch, PatternKind::Main],
            cut_width_mm: 0.35,
        }
    }

    #[test]
    fn sheets_cover_the_circumference() {
        let sheets = build_sheets(&doc()).unwrap();
        assert_eq!(sheets.len(), 2);
        let branch = &sheets[0];
        // Frame must span the full unwrapped circumference 2π·r2.
        assert!((branch.frame.2 - std::f64::consts::TAU * 35.0).abs() < 1e-6);
        // Branch development: open curve (u = 0 and u = 2πR meet on the tube).
        assert!(!branch.closed);
        // r2 < r1: the gueule de loup is a closed opening.
        assert!(sheets[1].closed);
    }

    #[test]
    fn tile_plan_is_minimal_and_sufficient() {
        let page = PageSpec {
            format: PageFormat::A4,
            orientation: Orientation::Landscape,
            margin_mm: 10.0,
        };
        // A4 landscape printable: 277 × 163 (with 24 mm title strip).
        let plan = plan_tiles(&page, 500.0, 120.0);
        assert_eq!(plan.rows, 1);
        assert!(plan.cols >= 2);
        // Coverage: last tile must reach past the drawing width.
        let covered = plan.view_w + (plan.cols - 1) as f64 * plan.step_x;
        assert!(covered >= 500.0);
        // Minimality: one tile less must not suffice.
        let covered_less = plan.view_w + (plan.cols - 2) as f64 * plan.step_x;
        assert!(covered_less < 500.0);
    }

    #[test]
    fn tile_labels_are_spreadsheet_style() {
        assert_eq!(tile_label(0, 0), "A1");
        assert_eq!(tile_label(2, 1), "C2");
        assert_eq!(tile_label(26, 0), "AA1");
    }

    #[test]
    fn missing_main_pattern_is_a_clean_error() {
        let mut d = doc();
        d.source = SourceSpec::CylPlane(CylPlaneInput {
            r1: 40.0,
            phi: 0.5,
            phi_y: 0.0,
            z0: 0.0,
            n_samples: 360,
        });
        d.patterns = vec![PatternKind::Main];
        assert!(matches!(
            build_sheets(&d),
            Err(ExportError::MissingPattern(PatternKind::Main))
        ));
    }
}
