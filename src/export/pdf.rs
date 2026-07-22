//! Minimal, dependency-free vector PDF writer specialised for the export
//! documents of this app.
//!
//! Design notes
//! ------------
//! * Every page installs a `mm → pt` CTM as its first operation, so all
//!   subsequent coordinates and line widths are expressed in millimetres.
//!   Combined with the page `MediaBox` derived from the ISO format, printing
//!   at "actual size" reproduces the patterns exactly at scale 1:1.
//! * Text uses the built-in Helvetica fonts with `WinAnsiEncoding`, which
//!   covers French accents and the drafting symbols we need (°, Ø, —).
//! * The writer emits PDF 1.4 with a classic xref table — readable by every
//!   viewer and CAD/print pipeline.

use super::{
    build_sheets, plan_tiles, tile_label, ExportDocument, ExportError, ScaleMode, Sheet,
    TITLE_BLOCK_H_MM,
};

const MM_TO_PT: f64 = 72.0 / 25.4;

// ---------------------------------------------------------------------------
// Content stream builder (all coordinates in mm, origin bottom-left, y up)
// ---------------------------------------------------------------------------

struct Cs {
    s: String,
}

impl Cs {
    fn new() -> Self {
        let mut s = String::with_capacity(16 * 1024);
        // Everything after this operates in millimetres.
        s.push_str(&format!("{MM_TO_PT:.10} 0 0 {MM_TO_PT:.10} 0 0 cm\n"));
        Cs { s }
    }

    fn push(&mut self, op: &str) {
        self.s.push_str(op);
        self.s.push('\n');
    }

    fn save(&mut self) {
        self.push("q");
    }
    fn restore(&mut self) {
        self.push("Q");
    }

    fn clip_rect(&mut self, x: f64, y: f64, w: f64, h: f64) {
        self.push(&format!("{x:.3} {y:.3} {w:.3} {h:.3} re W n"));
    }

    fn line_width(&mut self, mm: f64) {
        self.push(&format!("{mm:.4} w"));
    }

    fn stroke_gray(&mut self, g: f64) {
        self.push(&format!("{g:.3} G"));
    }

    fn stroke_rgb(&mut self, r: f64, g: f64, b: f64) {
        self.push(&format!("{r:.3} {g:.3} {b:.3} RG"));
    }

    fn fill_rgb(&mut self, r: f64, g: f64, b: f64) {
        self.push(&format!("{r:.3} {g:.3} {b:.3} rg"));
    }

    fn fill_gray(&mut self, g: f64) {
        self.push(&format!("{g:.3} g"));
    }

    fn dash(&mut self, on: f64, off: f64) {
        self.push(&format!("[{on:.3} {off:.3}] 0 d"));
    }

    fn solid(&mut self) {
        self.push("[] 0 d");
    }

    fn segment(&mut self, x1: f64, y1: f64, x2: f64, y2: f64) {
        self.push(&format!("{x1:.3} {y1:.3} m {x2:.3} {y2:.3} l S"));
    }

    fn polyline(&mut self, pts: &[(f64, f64)], closed: bool) {
        if pts.len() < 2 {
            return;
        }
        self.push(&format!("{:.3} {:.3} m", pts[0].0, pts[0].1));
        for &(x, y) in &pts[1..] {
            self.push(&format!("{x:.3} {y:.3} l"));
        }
        self.push(if closed { "s" } else { "S" });
    }

    fn rect_stroke(&mut self, x: f64, y: f64, w: f64, h: f64) {
        self.push(&format!("{x:.3} {y:.3} {w:.3} {h:.3} re S"));
    }

    /// `font`: 1 = Helvetica, 2 = Helvetica-Bold.  `size` in mm.
    fn text(&mut self, font: u8, size: f64, x: f64, y: f64, s: &str) {
        let encoded = encode_winansi(s);
        self.push(&format!(
            "BT /F{font} {size:.3} Tf {x:.3} {y:.3} Td ({encoded}) Tj ET"
        ));
    }
}

/// Encode a string as escaped WinAnsi bytes for a PDF literal string.
fn encode_winansi(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for ch in s.chars() {
        let byte: u8 = match ch {
            '\\' | '(' | ')' => {
                out.push('\\');
                ch as u8
            }
            c if (c as u32) < 0x80 => c as u8,
            // WinAnsi (CP-1252) codepoints we actually use.
            'à' => 0xE0, 'â' => 0xE2, 'ç' => 0xE7, 'è' => 0xE8, 'é' => 0xE9,
            'ê' => 0xEA, 'ë' => 0xEB, 'î' => 0xEE, 'ï' => 0xEF, 'ô' => 0xF4,
            'ö' => 0xF6, 'ù' => 0xF9, 'û' => 0xFB, 'ü' => 0xFC,
            'À' => 0xC0, 'É' => 0xC9, 'È' => 0xC8, 'Ê' => 0xCA,
            '°' => 0xB0, '²' => 0xB2, '·' => 0xB7, 'Ø' => 0xD8, 'ø' => 0xF8,
            '–' => 0x96, '—' => 0x97, '±' => 0xB1, '×' => 0xD7,
            _ => b'?',
        };
        if byte < 0x80 && !matches!(ch, '\\' | '(' | ')') {
            out.push(byte as char);
        } else if matches!(ch, '\\' | '(' | ')') {
            out.push(ch);
        } else {
            out.push_str(&format!("\\{byte:03o}"));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Document assembly
// ---------------------------------------------------------------------------

struct PageDef {
    width_mm: f64,
    height_mm: f64,
    content: String,
}

fn assemble(pages: &[PageDef]) -> Vec<u8> {
    let n_pages = pages.len();
    // Object numbering: 1 catalog, 2 pages, 3 F1, 4 F2, then per page
    // (page object, content stream).
    let total_objs = 4 + 2 * n_pages;
    let mut bodies: Vec<Vec<u8>> = Vec::with_capacity(total_objs);

    let kids: Vec<String> = (0..n_pages).map(|i| format!("{} 0 R", 5 + 2 * i)).collect();
    bodies.push(b"<< /Type /Catalog /Pages 2 0 R >>".to_vec());
    bodies.push(
        format!(
            "<< /Type /Pages /Kids [{}] /Count {} >>",
            kids.join(" "),
            n_pages
        )
        .into_bytes(),
    );
    bodies.push(
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>"
            .to_vec(),
    );
    bodies.push(
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>"
            .to_vec(),
    );

    for (i, page) in pages.iter().enumerate() {
        let content_obj = 6 + 2 * i;
        let w_pt = page.width_mm * MM_TO_PT;
        let h_pt = page.height_mm * MM_TO_PT;
        bodies.push(
            format!(
                "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {w_pt:.3} {h_pt:.3}] \
                 /Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> /Contents {content_obj} 0 R >>"
            )
            .into_bytes(),
        );
        let stream = page.content.as_bytes();
        let mut body = format!("<< /Length {} >>\nstream\n", stream.len()).into_bytes();
        body.extend_from_slice(stream);
        body.extend_from_slice(b"\nendstream");
        bodies.push(body);
    }

    let mut out: Vec<u8> = Vec::with_capacity(64 * 1024);
    out.extend_from_slice(b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");
    let mut offsets = Vec::with_capacity(bodies.len());
    for (i, body) in bodies.iter().enumerate() {
        offsets.push(out.len());
        out.extend_from_slice(format!("{} 0 obj\n", i + 1).as_bytes());
        out.extend_from_slice(body);
        out.extend_from_slice(b"\nendobj\n");
    }
    let xref_pos = out.len();
    out.extend_from_slice(format!("xref\n0 {}\n", bodies.len() + 1).as_bytes());
    out.extend_from_slice(b"0000000000 65535 f \n");
    for off in &offsets {
        out.extend_from_slice(format!("{off:010} 00000 n \n").as_bytes());
    }
    out.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{}\n%%EOF\n",
            bodies.len() + 1,
            xref_pos
        )
        .as_bytes(),
    );
    out
}

// ---------------------------------------------------------------------------
// Page rendering
// ---------------------------------------------------------------------------

/// Draw the model-space layers of a sheet (grid, frame, axes, cut,
/// annotations).  The caller has already installed the model transform.
fn draw_model(cs: &mut Cs, doc: &ExportDocument, sheet: &Sheet) {
    let (u0, v0, u1, v1) = sheet.bbox;

    if doc.layers.grid {
        cs.save();
        cs.solid();
        let start_u = (u0 / 10.0).floor() * 10.0;
        let start_v = (v0 / 10.0).floor() * 10.0;
        let mut u = start_u;
        while u <= u1 + 1e-9 {
            let major = (u / 50.0 - (u / 50.0).round()).abs() < 1e-9;
            cs.line_width(if major { 0.09 } else { 0.05 });
            cs.stroke_gray(if major { 0.72 } else { 0.86 });
            cs.segment(u, v0, u, v1);
            u += 10.0;
        }
        let mut v = start_v;
        while v <= v1 + 1e-9 {
            let major = (v / 50.0 - (v / 50.0).round()).abs() < 1e-9;
            cs.line_width(if major { 0.09 } else { 0.05 });
            cs.stroke_gray(if major { 0.72 } else { 0.86 });
            cs.segment(u0, v, u1, v);
            v += 10.0;
        }
        cs.restore();
    }

    if doc.layers.frame {
        cs.save();
        cs.line_width(0.15);
        cs.stroke_gray(0.45);
        cs.dash(1.6, 1.4);
        let (fx, fy, fw, fh) = sheet.frame;
        cs.rect_stroke(fx, fy, fw, fh);
        cs.restore();
    }

    if doc.layers.axis {
        // True tube generatrices (multiples of 90°, u = 0 ⇔ θ = 0): the
        // wrap-alignment marks — match them with lines traced on the tube.
        cs.save();
        cs.line_width(0.12);
        cs.stroke_rgb(0.20, 0.55, 0.70);
        cs.dash(3.0, 1.2);
        let (_fx, fy, _fw, fh) = sheet.frame;
        let gens = sheet.generatrices();
        for &(u, _) in &gens {
            cs.segment(u, fy, u, fy + fh);
        }

        // Axis-plane datum v = 0: longitudinal positioning reference.
        let (_, v0, _, v1) = sheet.bbox;
        if v0 < 0.0 && v1 > 0.0 {
            cs.stroke_gray(0.25);
            cs.dash(5.0, 1.6);
            cs.segment(sheet.bbox.0, 0.0, sheet.bbox.2, 0.0);
        }

        // Alignment ticks where each generatrix crosses the cut line.
        cs.solid();
        cs.line_width(0.3);
        cs.stroke_gray(0.0);
        for &(u, _) in &gens {
            for v in sheet.curve_crossings(u) {
                cs.segment(u - 2.5, v, u + 2.5, v);
            }
        }

        if doc.layers.labels {
            cs.fill_rgb(0.20, 0.55, 0.70);
            for &(u, deg) in &gens {
                cs.text(1, 2.4, u + 0.8, fy + 0.9, &format!("{deg}°"));
            }
            if v0 < 0.0 && v1 > 0.0 {
                cs.fill_gray(0.25);
                cs.text(1, 2.0, sheet.bbox.0 + 1.0, 0.7, "réf. plan des axes");
            }
        }
        cs.restore();
    }

    // The cut line itself — the pattern's accent colour, mirroring the 2D
    // view (kept dark enough to trace and cut confidently).
    cs.save();
    cs.solid();
    cs.line_width(doc.cut_width_mm);
    match sheet.kind {
        super::PatternKind::Branch => cs.stroke_rgb(0.80, 0.27, 0.04),
        super::PatternKind::Main => cs.stroke_rgb(0.03, 0.45, 0.62),
    }
    cs.polyline(&sheet.cut, sheet.closed);
    cs.restore();

    // Annotations.
    cs.save();
    cs.fill_gray(0.1);
    for a in &sheet.annotations {
        cs.text(1, a.size_mm.clamp(1.5, 20.0), a.u, a.v, &a.text);
    }
    cs.restore();
}

/// Bottom strip: identity, metadata, scale statement, tile info.
#[allow(clippy::too_many_arguments)]
fn draw_title_block(
    cs: &mut Cs,
    doc: &ExportDocument,
    sheet: &Sheet,
    page_w: f64,
    margin: f64,
    scale_note: &str,
    tile_note: &str,
    page_index: usize,
    page_count: usize,
) {
    let y0 = margin;
    let h = TITLE_BLOCK_H_MM - 4.0;
    let w = page_w - 2.0 * margin;

    cs.save();
    cs.solid();
    cs.line_width(0.25);
    cs.stroke_gray(0.0);
    cs.rect_stroke(margin, y0, w, h);

    // Vertical separators: identity | metadata | scale/tile.
    let x1 = margin + w * 0.42;
    let x2 = margin + w * 0.74;
    cs.line_width(0.12);
    cs.segment(x1, y0, x1, y0 + h);
    cs.segment(x2, y0, x2, y0 + h);

    cs.fill_gray(0.0);
    cs.text(2, 3.4, margin + 3.0, y0 + h - 5.2, &sheet.name);
    cs.fill_gray(0.25);
    cs.text(1, 2.6, margin + 3.0, y0 + h - 9.6, &sheet.meta);
    if !doc.title_block.notes.is_empty() {
        // Kept above the control ruler strip to avoid any overlap.
        cs.text(1, 2.4, margin + 3.0, y0 + h - 13.8, &doc.title_block.notes);
    }

    let mut ty = y0 + h - 5.2;
    for (label, value) in [
        ("Projet", &doc.title_block.project),
        ("Titre", &doc.title_block.title),
        ("Auteur", &doc.title_block.author),
        ("Date", &doc.title_block.date),
    ] {
        if value.is_empty() {
            continue;
        }
        cs.fill_gray(0.45);
        cs.text(1, 2.2, x1 + 3.0, ty, label);
        cs.fill_gray(0.0);
        cs.text(1, 2.6, x1 + 17.0, ty, value);
        ty -= 4.4;
    }

    cs.fill_gray(0.0);
    cs.text(2, 3.2, x2 + 3.0, y0 + h - 5.2, scale_note);
    if !tile_note.is_empty() {
        cs.text(2, 4.6, x2 + 3.0, y0 + h - 12.0, tile_note);
    }
    cs.fill_gray(0.35);
    cs.text(
        1,
        2.2,
        x2 + 3.0,
        y0 + 2.2,
        &format!("page {page_index}/{page_count} — cylinders-intersection"),
    );
    cs.restore();

    // 100 mm control ruler (only meaningful at 1:1, where it is exact).
    if doc.layers.scale_bar && matches!(doc.scale, ScaleMode::OneToOne) {
        let rx = margin + 3.0;
        let ry = y0 + 3.4;
        cs.save();
        cs.line_width(0.3);
        cs.stroke_gray(0.0);
        cs.segment(rx, ry, rx + 100.0, ry);
        for k in 0..=10 {
            let x = rx + 10.0 * k as f64;
            let tick = if k % 5 == 0 { 1.6 } else { 1.0 };
            cs.segment(x, ry, x, ry + tick);
        }
        cs.fill_gray(0.25);
        cs.text(1, 2.0, rx + 102.0, ry - 0.6, "100 mm");
        cs.restore();
    }
}

fn scale_note_for(scale: f64) -> String {
    if (scale - 1.0).abs() < 1e-9 {
        "Échelle 1:1".to_string()
    } else if scale < 1.0 {
        format!("Échelle 1:{:.2} — aperçu, ne pas tracer", 1.0 / scale)
    } else {
        format!("Échelle {:.2}:1 — aperçu, ne pas tracer", scale)
    }
}

fn render_fit_page(doc: &ExportDocument, sheet: &Sheet, idx: usize, count: usize) -> PageDef {
    let (pw, ph) = doc.page.size_mm();
    let margin = doc.page.margin_mm;
    let view_w = pw - 2.0 * margin;
    let view_h = ph - 2.0 * margin - TITLE_BLOCK_H_MM;

    // Fit the useful drawing (cut + annotations), not the empty frame.
    let (u0, v0, u1, v1) = sheet.cut_bbox;
    let w = (u1 - u0).max(1e-6);
    let h = (v1 - v0).max(1e-6);
    let s = (view_w / w).min(view_h / h).min(1.0);

    let mut cs = Cs::new();
    // Center the drawing in the viewport.
    let ox = margin + (view_w - w * s) / 2.0 - u0 * s;
    let oy = margin + TITLE_BLOCK_H_MM + (view_h - h * s) / 2.0 - v0 * s;
    cs.save();
    cs.clip_rect(margin, margin + TITLE_BLOCK_H_MM, view_w, view_h);
    cs.push(&format!("{s:.6} 0 0 {s:.6} {ox:.3} {oy:.3} cm"));
    // Compensate line widths so they stay visually constant on paper.
    let mut scaled = doc.clone();
    scaled.cut_width_mm = doc.cut_width_mm / s;
    draw_model(&mut cs, &scaled, sheet);
    cs.restore();

    draw_title_block(
        &mut cs, doc, sheet, pw, margin, &scale_note_for(s), "", idx, count,
    );

    PageDef { width_mm: pw, height_mm: ph, content: cs.s }
}

#[allow(clippy::too_many_arguments)]
fn render_tile_page(
    doc: &ExportDocument,
    sheet: &Sheet,
    col: usize,
    row: usize,
    plan: super::TilePlan,
    idx: usize,
    count: usize,
) -> PageDef {
    let (pw, ph) = doc.page.size_mm();
    let margin = doc.page.margin_mm;
    let (u0, _v0, _u1, v1) = sheet.cut_bbox;

    // Row 0 is the top row so labels read naturally on the assembled wall.
    let tile_u = u0 + col as f64 * plan.step_x;
    let tile_v_top = v1 - row as f64 * plan.step_y;
    let tile_v = tile_v_top - plan.view_h;

    let view_x = margin;
    let view_y = margin + TITLE_BLOCK_H_MM;

    let mut cs = Cs::new();
    cs.save();
    cs.clip_rect(view_x, view_y, plan.view_w, plan.view_h);
    cs.push(&format!(
        "1 0 0 1 {:.3} {:.3} cm",
        view_x - tile_u,
        view_y - tile_v
    ));
    draw_model(&mut cs, doc, sheet);

    // Glue overlap indicators, drawn in model space inside the clip.
    cs.line_width(0.12);
    cs.stroke_gray(0.55);
    cs.dash(4.0, 2.0);
    if col + 1 < plan.cols {
        let x = tile_u + plan.step_x;
        cs.segment(x, tile_v, x, tile_v_top);
    }
    if row + 1 < plan.rows {
        let y = tile_v + plan.overlap;
        cs.segment(tile_u, y, tile_u + plan.view_w, y);
    }
    cs.restore();

    // Crop marks at the viewport corners.
    cs.save();
    cs.solid();
    cs.line_width(0.15);
    cs.stroke_gray(0.0);
    let corners = [
        (view_x, view_y),
        (view_x + plan.view_w, view_y),
        (view_x, view_y + plan.view_h),
        (view_x + plan.view_w, view_y + plan.view_h),
    ];
    for (cx, cy) in corners {
        cs.segment(cx - 4.0, cy, cx + 4.0, cy);
        cs.segment(cx, cy - 4.0, cx, cy + 4.0);
    }
    cs.restore();

    let label = tile_label(col, row);
    let neighbors = {
        let mut parts: Vec<String> = Vec::new();
        if col + 1 < plan.cols {
            parts.push(format!("droite: {}", tile_label(col + 1, row)));
        }
        if row + 1 < plan.rows {
            parts.push(format!("dessous: {}", tile_label(col, row + 1)));
        }
        if parts.is_empty() {
            String::new()
        } else {
            format!("raccords — {}", parts.join(" · "))
        }
    };
    // Neighbor note just above the title strip, useful while gluing.
    if !neighbors.is_empty() {
        cs.fill_gray(0.4);
        cs.text(1, 2.2, view_x + 1.0, view_y + 1.2, &neighbors);
    }

    draw_title_block(
        &mut cs,
        doc,
        sheet,
        pw,
        margin,
        "Échelle 1:1",
        &format!("Tuile {label}"),
        idx,
        count,
    );

    PageDef { width_mm: pw, height_mm: ph, content: cs.s }
}

/// Sheet index plus its optional tile placement `(col, row, plan)`.
type PageLayout = (usize, Option<(usize, usize, super::TilePlan)>);

/// Render the whole export document to PDF bytes.
pub fn render_pdf(doc: &ExportDocument) -> Result<Vec<u8>, ExportError> {
    let sheets = build_sheets(doc)?;

    // First pass: count pages for the "page i/n" footer.
    let mut layouts: Vec<PageLayout> = Vec::new();
    let mut count = 0usize;
    for (si, sheet) in sheets.iter().enumerate() {
        match doc.scale {
            ScaleMode::Fit => {
                count += 1;
                layouts.push((si, None));
            }
            ScaleMode::OneToOne => {
                // Tile over the cut extents only, and drop tiles where no
                // cut line nor annotation lands: nobody glues empty grids.
                let (u0, v0, u1, v1) = sheet.cut_bbox;
                let plan = plan_tiles(&doc.page, u1 - u0, v1 - v0);
                for row in 0..plan.rows {
                    for col in 0..plan.cols {
                        let tile_u = u0 + col as f64 * plan.step_x;
                        let tile_v_top = v1 - row as f64 * plan.step_y;
                        if !sheet.rect_has_content(
                            tile_u,
                            tile_v_top - plan.view_h,
                            plan.view_w,
                            plan.view_h,
                            plan.overlap,
                        ) {
                            continue;
                        }
                        count += 1;
                        layouts.push((si, Some((col, row, plan))));
                    }
                }
            }
        }
    }

    let mut pages = Vec::with_capacity(count);
    for (idx, (si, tile)) in layouts.iter().enumerate() {
        let sheet = &sheets[*si];
        let page = match tile {
            None => render_fit_page(doc, sheet, idx + 1, count),
            Some((col, row, plan)) => {
                render_tile_page(doc, sheet, *col, *row, *plan, idx + 1, count)
            }
        };
        pages.push(page);
    }

    Ok(assemble(&pages))
}

#[cfg(test)]
mod tests {
    use super::super::{
        Layers, Orientation, PageFormat, PageSpec, PatternKind, SourceSpec, TitleBlock,
    };
    use super::*;
    use crate::geometry::Branch;
    use crate::intersection::CylCylInput;

    fn doc(scale: ScaleMode) -> ExportDocument {
        ExportDocument {
            source: SourceSpec::CylCyl(CylCylInput {
                r1: 80.0,
                r2: 60.0,
                phi: 0.9,
                n_samples: 720,
                branch: Branch::Outer,
            }),
            page: PageSpec {
                format: PageFormat::A4,
                orientation: Orientation::Landscape,
                margin_mm: 10.0,
            },
            scale,
            layers: Layers::default(),
            title_block: TitleBlock {
                title: "Piquage Ø160/Ø120".into(),
                project: "Démo".into(),
                author: "Atelier".into(),
                date: "2026-07-22".into(),
                notes: "Découpe plasma — vérifier la règle 100 mm.".into(),
            },
            annotations: vec![],
            patterns: vec![PatternKind::Branch, PatternKind::Main],
            cut_width_mm: 0.35,
        }
    }

    #[test]
    fn pdf_has_valid_skeleton() {
        let bytes = render_pdf(&doc(ScaleMode::Fit)).unwrap();
        assert!(bytes.starts_with(b"%PDF-1.4"));
        assert!(bytes.ends_with(b"%%EOF\n"));
        let text = String::from_utf8_lossy(&bytes);
        assert!(text.contains("/Type /Catalog"));
        assert!(text.contains("/Count 2"), "fit mode: one page per pattern");
    }

    #[test]
    fn one_to_one_produces_tiles() {
        // Ø160 branch unwraps to ~377 mm — wider than an A4 landscape
        // printable width, so at least 2 columns are required.
        let bytes = render_pdf(&doc(ScaleMode::OneToOne)).unwrap();
        let text = String::from_utf8_lossy(&bytes);
        let pages = text.matches("/Type /Page ").count();
        // Branch sinusoid spans the full 377 mm circumference → 2 columns;
        // the gueule de loup egg fits on a single sheet.  Every page carries
        // actual cut geometry.
        assert_eq!(pages, 3, "expected 3 useful pages, got {pages}");
        assert!(text.contains("(Tuile A1)"));
        assert!(text.contains("(Tuile B1)"));
    }

    #[test]
    fn xref_offsets_point_at_objects() {
        // Work on raw bytes: the binary header comment is not valid UTF-8,
        // so string conversions would shift every offset.
        let bytes = render_pdf(&doc(ScaleMode::Fit)).unwrap();
        let find = |needle: &[u8], from: usize| -> usize {
            bytes[from..]
                .windows(needle.len())
                .position(|w| w == needle)
                .map(|p| p + from)
                .unwrap()
        };
        let sx = find(b"startxref\n", 0);
        let tail = &bytes[sx + 10..];
        let line_end = tail.iter().position(|&b| b == b'\n').unwrap();
        let offset: usize = std::str::from_utf8(&tail[..line_end]).unwrap().parse().unwrap();
        assert!(bytes[offset..].starts_with(b"xref"));
        // First object offset must land exactly on "1 0 obj".
        let first = find(b"\n0000000000 65535 f \n", 0) + 21;
        let entry = std::str::from_utf8(&bytes[first..first + 10]).unwrap();
        let obj1: usize = entry.parse().unwrap();
        assert!(bytes[obj1..].starts_with(b"1 0 obj"));
    }

    #[test]
    fn tiling_covers_the_cut_not_the_empty_frame() {
        // Ø160 main / Ø60 branch: the gueule de loup egg occupies ~60 mm of
        // the 502 mm unwrapped circumference.  Tiling the frame would burn
        // 2+ pages on empty grid; tiling the cut box needs exactly one.
        let mut d = doc(ScaleMode::OneToOne);
        d.source = SourceSpec::CylCyl(CylCylInput {
            r1: 80.0,
            r2: 30.0,
            phi: 1.0,
            n_samples: 720,
            branch: Branch::Outer,
        });
        let bytes = render_pdf(&d).unwrap();
        let text = String::from_utf8_lossy(&bytes);
        let pages = text.matches("/Type /Page ").count();
        // Branch sinusoid (188 mm wide) → 1 page; egg (~60 mm) → 1 page.
        assert_eq!(pages, 2, "expected 2 useful pages, got {pages}");
    }

    #[test]
    fn winansi_escapes_accents_and_parens() {
        let e = encode_winansi("é (test) — 90°");
        assert!(e.contains("\\351"));
        assert!(e.contains("\\(test\\)"));
        assert!(e.contains("\\227"));
        assert!(e.contains("\\260"));
    }
}
