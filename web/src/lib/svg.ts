import type { DevPoint } from "./api";

interface PatternMeta {
  title: string;
  subtitle: string;
  diameter: number;       // mm — for the legend
  circumference: number;  // mm — full unwrapped width
  closed: boolean;        // close back to the first point?
}

/**
 * Build a 1:1 scale SVG document of a developed cut pattern.  The SVG is
 * sized in mm so that any printer with the proper "Actual size" / "Échelle
 * 100%" setting reproduces the gabarit at scale 1:1.
 */
export function buildPatternSVG(
  points: DevPoint[],
  meta: PatternMeta,
): string {
  if (points.length === 0) {
    return `<svg xmlns="http://www.w3.org/2000/svg" width="0mm" height="0mm" viewBox="0 0 0 0"></svg>`;
  }

  const padding = 12; // mm of margin around the pattern
  const ruler = 10;   // mm — extra strip below for the "100 mm" reference

  let uMin = points[0].u, uMax = uMin;
  let vMin = points[0].v, vMax = vMin;
  for (const p of points) {
    if (p.u < uMin) uMin = p.u;
    if (p.u > uMax) uMax = p.u;
    if (p.v < vMin) vMin = p.v;
    if (p.v > vMax) vMax = p.v;
  }
  // Pattern goes over the full circumference width:
  if (meta.circumference > uMax - uMin) {
    uMax = uMin + meta.circumference;
  }

  const widthMm = uMax - uMin + padding * 2;
  const heightMm = vMax - vMin + padding * 2 + ruler + 6;

  // Build the polyline path (Y is inverted: SVG has Y-down, our v is Y-up).
  const pathD = points
    .map((p, i) => {
      const x = (p.u - uMin) + padding;
      const y = heightMm - ruler - 6 - ((p.v - vMin) + padding);
      return `${i === 0 ? "M" : "L"}${x.toFixed(3)} ${y.toFixed(3)}`;
    })
    .join(" ");
  const close = meta.closed ? " Z" : "";

  // Draw the rectangle representing the unwrapped cylinder footprint:
  const rectX = padding;
  const rectW = meta.circumference;
  const rectY = heightMm - ruler - 6 - ((vMax - vMin) + padding);
  const rectH = vMax - vMin + 0.0001;

  // Reference 100 mm scale bar at the bottom-left:
  const ref = 100;
  const refX1 = padding;
  const refX2 = padding + ref;
  const refY = heightMm - 4;

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="${widthMm.toFixed(3)}mm" height="${heightMm.toFixed(3)}mm"
     viewBox="0 0 ${widthMm.toFixed(3)} ${heightMm.toFixed(3)}"
     shape-rendering="geometricPrecision">
  <title>${escapeXml(meta.title)}</title>
  <desc>${escapeXml(meta.subtitle)}</desc>
  <style>
    .frame { fill: none; stroke: #999; stroke-width: 0.15; stroke-dasharray: 1.2 1.2; }
    .cut   { fill: none; stroke: #ff5b1a; stroke-width: 0.25; stroke-linejoin: round; stroke-linecap: round; }
    .axis  { fill: none; stroke: #555; stroke-width: 0.1; stroke-dasharray: 0.8 0.8; }
    .label { font: 3px "Inter", sans-serif; fill: #222; }
    .small { font: 2px "Inter", sans-serif; fill: #666; }
    .scale { fill: none; stroke: #111; stroke-width: 0.25; }
    .scaleTick { fill: none; stroke: #111; stroke-width: 0.25; }
  </style>

  <rect x="0" y="0" width="${widthMm.toFixed(3)}" height="${heightMm.toFixed(3)}" fill="#ffffff"/>

  <!-- unwrapped cylinder rectangle -->
  <rect class="frame" x="${rectX.toFixed(3)}" y="${rectY.toFixed(3)}" width="${rectW.toFixed(3)}" height="${rectH.toFixed(3)}"/>

  <!-- reference 0 axis (theta = 0) -->
  <line class="axis" x1="${rectX.toFixed(3)}" y1="${rectY.toFixed(3)}" x2="${rectX.toFixed(3)}" y2="${(rectY + rectH).toFixed(3)}"/>

  <!-- the cut curve -->
  <path class="cut" d="${pathD}${close}"/>

  <!-- legend -->
  <text class="label" x="${padding}" y="${(padding - 4).toFixed(1)}">
    ${escapeXml(meta.title)} — ${escapeXml(meta.subtitle)}
  </text>
  <text class="small" x="${padding}" y="${(padding - 1).toFixed(1)}">
    Ø ${meta.diameter.toFixed(2)} mm   |   périmètre = ${meta.circumference.toFixed(2)} mm   |   échelle 1:1
  </text>

  <!-- scale bar -->
  <line class="scale" x1="${refX1}" y1="${refY.toFixed(1)}" x2="${refX2}" y2="${refY.toFixed(1)}"/>
  <line class="scaleTick" x1="${refX1}" y1="${(refY - 1.2).toFixed(1)}" x2="${refX1}" y2="${(refY + 1.2).toFixed(1)}"/>
  <line class="scaleTick" x1="${refX2}" y1="${(refY - 1.2).toFixed(1)}" x2="${refX2}" y2="${(refY + 1.2).toFixed(1)}"/>
  <text class="small" x="${(refX1 + 2).toFixed(1)}" y="${(refY - 1.6).toFixed(1)}">100 mm de référence — vérifier après impression</text>
</svg>`;
}

export function downloadSVG(filename: string, svg: string) {
  const blob = new Blob([svg], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 0);
}

export function downloadCSV(filename: string, rows: string[][]) {
  const csv = rows.map((r) => r.map(csvEscape).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 0);
}

function csvEscape(s: string): string {
  if (/[",\n]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
  return s;
}

function escapeXml(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&apos;");
}
