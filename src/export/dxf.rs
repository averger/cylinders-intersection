//! DXF R12 (AC1009) ASCII writer — the most widely importable flavour of
//! DXF, accepted by AutoCAD, QCAD, LibreCAD, laser/plasma CAM pipelines…
//!
//! Layout conventions
//! ------------------
//! * Model space is in millimetres, patterns are placed side by side with a
//!   50 mm gap, each with its own origin marker.
//! * Layers: `CUT` (the toolpath, colour white), `FRAME` (unwrapped tube
//!   footprint, green), `AXIS` (quarter generators, cyan), `TEXT` (labels,
//!   yellow), `ANNOT` (user annotations, magenta).
//! * The on-screen grid is intentionally not exported: CAD tools bring
//!   their own grids, and CAM pipelines must only ever see the `CUT` layer.

use super::{build_sheets, ExportDocument, ExportError, Sheet};

const GAP_MM: f64 = 50.0;

struct DxfWriter {
    out: String,
}

impl DxfWriter {
    fn new() -> Self {
        DxfWriter { out: String::with_capacity(64 * 1024) }
    }

    fn tag(&mut self, code: i32, value: &str) {
        self.out.push_str(&format!("{code}\n{value}\n"));
    }

    fn tag_f(&mut self, code: i32, value: f64) {
        self.tag(code, &format!("{value:.4}"));
    }

    fn header(&mut self) {
        self.tag(0, "SECTION");
        self.tag(2, "HEADER");
        self.tag(9, "$ACADVER");
        self.tag(1, "AC1009");
        // Units flag: metric.
        self.tag(9, "$MEASUREMENT");
        self.tag(70, "1");
        self.tag(0, "ENDSEC");
    }

    fn tables(&mut self, layers: &[(&str, i32)]) {
        self.tag(0, "SECTION");
        self.tag(2, "TABLES");

        self.tag(0, "TABLE");
        self.tag(2, "LTYPE");
        self.tag(70, "1");
        self.tag(0, "LTYPE");
        self.tag(2, "CONTINUOUS");
        self.tag(70, "0");
        self.tag(3, "Solid line");
        self.tag(72, "65");
        self.tag(73, "0");
        self.tag(40, "0.0");
        self.tag(0, "ENDTAB");

        self.tag(0, "TABLE");
        self.tag(2, "LAYER");
        self.tag(70, &layers.len().to_string());
        for (name, color) in layers {
            self.tag(0, "LAYER");
            self.tag(2, name);
            self.tag(70, "0");
            self.tag(62, &color.to_string());
            self.tag(6, "CONTINUOUS");
        }
        self.tag(0, "ENDTAB");

        self.tag(0, "ENDSEC");
    }

    fn begin_entities(&mut self) {
        self.tag(0, "SECTION");
        self.tag(2, "ENTITIES");
    }

    fn polyline(&mut self, layer: &str, pts: &[(f64, f64)], closed: bool) {
        if pts.len() < 2 {
            return;
        }
        self.tag(0, "POLYLINE");
        self.tag(8, layer);
        self.tag(66, "1");
        self.tag(70, if closed { "1" } else { "0" });
        for &(x, y) in pts {
            self.tag(0, "VERTEX");
            self.tag(8, layer);
            self.tag_f(10, x);
            self.tag_f(20, y);
            self.tag_f(30, 0.0);
        }
        self.tag(0, "SEQEND");
    }

    fn line(&mut self, layer: &str, x1: f64, y1: f64, x2: f64, y2: f64) {
        self.tag(0, "LINE");
        self.tag(8, layer);
        self.tag_f(10, x1);
        self.tag_f(20, y1);
        self.tag_f(30, 0.0);
        self.tag_f(11, x2);
        self.tag_f(21, y2);
        self.tag_f(31, 0.0);
    }

    fn rect(&mut self, layer: &str, x: f64, y: f64, w: f64, h: f64) {
        let pts = [(x, y), (x + w, y), (x + w, y + h), (x, y + h)];
        self.polyline(layer, &pts, true);
    }

    fn text(&mut self, layer: &str, x: f64, y: f64, height: f64, value: &str) {
        self.tag(0, "TEXT");
        self.tag(8, layer);
        self.tag_f(10, x);
        self.tag_f(20, y);
        self.tag_f(30, 0.0);
        self.tag_f(40, height);
        self.tag(1, &sanitize(value));
    }

    fn finish(mut self) -> String {
        self.tag(0, "ENDSEC");
        self.tag(0, "EOF");
        self.out
    }
}

/// R12 text is codepage-bound; keep exports portable by transliterating to
/// plain ASCII.
fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            'à' | 'â' => 'a',
            'ç' => 'c',
            'é' | 'è' | 'ê' | 'ë' => 'e',
            'î' | 'ï' => 'i',
            'ô' | 'ö' => 'o',
            'ù' | 'û' | 'ü' => 'u',
            'É' | 'È' | 'Ê' => 'E',
            'À' => 'A',
            'Ø' | 'ø' => 'D',
            '°' => 'd',
            '—' | '–' | '·' => '-',
            c if c.is_ascii() => c,
            _ => '?',
        })
        .collect()
}

fn write_sheet(w: &mut DxfWriter, sheet: &Sheet, offset_x: f64) {
    let (u0, v0, u1, v1) = sheet.bbox;
    let dx = offset_x - u0;

    let (fx, fy, fw, fh) = sheet.frame;
    w.rect("FRAME", fx + dx, fy, fw, fh);

    // True generatrices with angle labels, alignment ticks on the cut line,
    // and the axis-plane datum — the wrap-alignment marks.
    for (g, deg) in sheet.generatrices() {
        // θ = 0 goes to the dedicated REF layer (wrap-alignment datum).
        let layer = if deg == 0 { "REF" } else { "AXIS" };
        w.line(layer, g + dx, fy, g + dx, fy + fh);
        w.text(
            layer,
            g + dx + 0.8,
            fy + 0.9,
            2.4,
            &if deg == 0 { "0d REF".to_string() } else { format!("{deg}d") },
        );
        for v in sheet.curve_crossings(g) {
            w.line(layer, g + dx - 2.5, v, g + dx + 2.5, v);
        }
    }
    if sheet.bbox.1 < 0.0 && sheet.bbox.3 > 0.0 {
        w.line("AXIS", sheet.bbox.0 + dx, 0.0, sheet.bbox.2 + dx, 0.0);
        w.text("AXIS", sheet.bbox.0 + dx + 1.0, 0.7, 2.0, "ref plan des axes");
    }

    let cut: Vec<(f64, f64)> = sheet.cut.iter().map(|&(u, v)| (u + dx, v)).collect();
    w.polyline("CUT", &cut, sheet.closed);
    for (pts, closed) in &sheet.holes {
        let shifted: Vec<(f64, f64)> = pts.iter().map(|&(u, v)| (u + dx, v)).collect();
        w.polyline("CUT", &shifted, *closed);
    }

    w.text("TEXT", u0 + dx, v1 + 8.0, 5.0, &sheet.name);
    w.text("TEXT", u0 + dx, v1 + 2.0, 3.0, &sheet.meta);

    for a in &sheet.annotations {
        w.text("ANNOT", a.u + dx, a.v, a.size_mm.clamp(1.5, 20.0), &a.text);
    }

    // Origin marker of this pattern (u = 0, v = 0 in pattern coordinates).
    w.line("AXIS", dx - 3.0, 0.0, dx + 3.0, 0.0);
    w.line("AXIS", dx, -3.0, dx, 3.0);

    let _ = (v0, u1);
}

/// Render the export document as a DXF R12 string.
pub fn render_dxf(doc: &ExportDocument) -> Result<String, ExportError> {
    let sheets = build_sheets(doc)?;

    let mut w = DxfWriter::new();
    w.header();
    w.tables(&[
        ("CUT", 7),
        ("FRAME", 3),
        ("AXIS", 4),
        ("REF", 1),
        ("TEXT", 2),
        ("ANNOT", 6),
    ]);
    w.begin_entities();

    let mut offset_x = 0.0;
    for sheet in &sheets {
        write_sheet(&mut w, sheet, offset_x);
        let (u0, _, u1, _) = sheet.bbox;
        offset_x += (u1 - u0) + GAP_MM;
    }

    Ok(w.finish())
}

#[cfg(test)]
mod tests {
    use super::super::{
        Layers, Orientation, PageFormat, PageSpec, PatternKind, ScaleMode, SourceSpec, TitleBlock,
    };
    use super::*;
    use crate::geometry::Branch;
    use crate::intersection::CylCylInput;

    fn doc() -> ExportDocument {
        ExportDocument {
            source: SourceSpec::CylCyl(CylCylInput {
                r1: 50.0,
                r2: 35.0,
                phi: 1.1,
                n_samples: 360,
                branch: Branch::Outer,
            }),
            page: PageSpec {
                format: PageFormat::A3,
                orientation: Orientation::Landscape,
                margin_mm: 10.0,
            },
            scale: ScaleMode::OneToOne,
            layers: Layers::default(),
            title_block: TitleBlock::default(),
            annotations: vec![super::super::Annotation {
                pattern: PatternKind::Branch,
                u: 30.0,
                v: 10.0,
                text: "repère soudure".into(),
                size_mm: 4.0,
            }],
            patterns: vec![PatternKind::Branch, PatternKind::Main],
            cut_width_mm: 0.35,
        }
    }

    #[test]
    fn dxf_is_structurally_sound() {
        let dxf = render_dxf(&doc()).unwrap();
        assert!(dxf.starts_with("0\nSECTION\n2\nHEADER\n"));
        assert!(dxf.ends_with("0\nEOF\n"));
        assert_eq!(dxf.matches("\nPOLYLINE\n").count(), 4, "2 cuts + 2 frames");
        assert_eq!(
            dxf.matches("\nSEQEND\n").count(),
            dxf.matches("\nPOLYLINE\n").count()
        );
        for layer in ["CUT", "FRAME", "AXIS", "TEXT", "ANNOT"] {
            assert!(dxf.contains(&format!("\n{layer}\n")), "layer {layer}");
        }
    }

    #[test]
    fn dxf_text_is_ascii_only() {
        let dxf = render_dxf(&doc()).unwrap();
        assert!(dxf.is_ascii(), "R12 export must stay plain ASCII");
        assert!(dxf.contains("repere soudure"));
    }

    #[test]
    fn closed_cut_uses_closed_polyline_flag() {
        let dxf = render_dxf(&doc()).unwrap();
        // At least one POLYLINE with the closed flag (70 → 1) exists.
        assert!(dxf.contains("POLYLINE\n8\nCUT\n66\n1\n70\n1\n"));
    }
}
