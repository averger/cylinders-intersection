<script lang="ts">
  import type { DevPoint } from "../lib/api";
  import { buildPatternSVG, downloadSVG, downloadCSV } from "../lib/svg";

  interface Props {
    title: string;
    subtitle: string;
    points: DevPoint[];
    diameter: number;
    circumference: number;
    closed: boolean;
    accent?: string;
    filename: string;
  }

  let {
    title,
    subtitle,
    points,
    diameter,
    circumference,
    closed,
    accent = "#ff5b1a",
    filename,
  }: Props = $props();

  // 1 mm = 1 SVG user unit. We compute a viewBox in mm and let the page CSS
  // scale it to fit the on-screen card. The downloaded SVG keeps the mm sizing.
  let bbox = $derived.by(() => {
    if (points.length === 0)
      return { uMin: 0, uMax: circumference, vMin: 0, vMax: 1 };
    let uMin = points[0].u, uMax = uMin, vMin = points[0].v, vMax = vMin;
    for (const p of points) {
      if (p.u < uMin) uMin = p.u;
      if (p.u > uMax) uMax = p.u;
      if (p.v < vMin) vMin = p.v;
      if (p.v > vMax) vMax = p.v;
    }
    if (circumference > uMax - uMin) uMax = uMin + circumference;
    return { uMin, uMax, vMin, vMax };
  });

  let pad = 6; // mm of margin inside the visible viewport
  let width = $derived(bbox.uMax - bbox.uMin + pad * 2);
  let height = $derived(Math.max(bbox.vMax - bbox.vMin + pad * 2, 24));

  let path = $derived.by(() => {
    if (points.length === 0) return "";
    const d = points
      .map((p, i) => {
        const x = p.u - bbox.uMin + pad;
        const y = height - (p.v - bbox.vMin + pad);
        return `${i === 0 ? "M" : "L"}${x.toFixed(3)} ${y.toFixed(3)}`;
      })
      .join(" ");
    return closed ? `${d} Z` : d;
  });

  function exportSVG() {
    const svg = buildPatternSVG(points, {
      title,
      subtitle,
      diameter,
      circumference,
      closed,
    });
    downloadSVG(`${filename}.svg`, svg);
  }

  function exportCSV() {
    const rows = [["theta_rad", "u_mm", "v_mm"]];
    for (const p of points) {
      rows.push([p.theta.toFixed(6), p.u.toFixed(4), p.v.toFixed(4)]);
    }
    downloadCSV(`${filename}.csv`, rows);
  }

  // Smart tick generator: pick a step that yields ~6-12 ticks.
  function ticks(min: number, max: number): number[] {
    const span = max - min;
    if (span <= 0) return [];
    const target = 8;
    const raw = span / target;
    const pow10 = Math.pow(10, Math.floor(Math.log10(raw)));
    const candidates = [1, 2, 5, 10].map((m) => m * pow10);
    let step = candidates[0];
    for (const c of candidates) if (Math.abs(c - raw) < Math.abs(step - raw)) step = c;
    const out: number[] = [];
    const start = Math.ceil(min / step) * step;
    for (let v = start; v <= max + 1e-6; v += step) out.push(v);
    return out;
  }

  let uTicks = $derived(ticks(bbox.uMin, bbox.uMax));
  let vTicks = $derived(ticks(bbox.vMin, bbox.vMax));
</script>

<section class="glass card-hairline p-5 flex flex-col gap-3 min-h-0">
  <header class="flex items-start justify-between gap-3">
    <div>
      <span class="pill" style="border-color: {accent}55; color: {accent}">{title}</span>
      <h3 class="text-pearl text-lg font-semibold tracking-tight mt-2">{subtitle}</h3>
      <p class="text-[11px] text-ash mt-0.5 num">
        Ø {diameter.toFixed(2)} mm · périmètre {circumference.toFixed(2)} mm
        · largeur dévelop. {(bbox.uMax - bbox.uMin).toFixed(2)} mm
        · hauteur {(bbox.vMax - bbox.vMin).toFixed(2)} mm
      </p>
    </div>
    <div class="flex gap-2 no-print">
      <button class="btn-ghost" onclick={exportCSV} aria-label="Export CSV">
        <span class="text-[11px] uppercase tracking-[0.18em]">CSV</span>
      </button>
      <button class="btn-ghost" onclick={exportSVG} aria-label="Export SVG">
        <span class="text-[11px] uppercase tracking-[0.18em]">SVG 1:1</span>
      </button>
    </div>
  </header>

  <div class="relative flex-1 min-h-[180px] rounded-[12px] overflow-hidden bg-black/30 border border-white/5 grid-bg">
    <svg
      viewBox="0 0 {width} {height}"
      preserveAspectRatio="xMidYMid meet"
      class="w-full h-full"
      xmlns="http://www.w3.org/2000/svg"
    >
      <!-- bounding rectangle of the unwrapped cylinder -->
      <rect
        x={pad}
        y={height - (bbox.vMax - bbox.vMin) - pad}
        width={circumference}
        height={bbox.vMax - bbox.vMin + 0.0001}
        fill="rgba(255,255,255,0.018)"
        stroke="rgba(255,255,255,0.18)"
        stroke-width={width / 600}
        stroke-dasharray="{(width / 200).toFixed(2)} {(width / 200).toFixed(2)}"
      />

      <!-- u ticks -->
      {#each uTicks as t (t)}
        <g>
          <line
            x1={t - bbox.uMin + pad}
            y1={height - pad}
            x2={t - bbox.uMin + pad}
            y2={pad}
            stroke="rgba(255,255,255,0.04)"
            stroke-width={width / 1500}
          />
          <text
            x={t - bbox.uMin + pad}
            y={height - pad / 2}
            fill="rgba(184,184,192,0.8)"
            font-size={Math.max(2, height / 22)}
            text-anchor="middle"
            font-family="JetBrains Mono, monospace"
          >
            {t.toFixed(0)}
          </text>
        </g>
      {/each}

      <!-- v ticks -->
      {#each vTicks as t (t)}
        <g>
          <line
            x1={pad}
            y1={height - (t - bbox.vMin + pad)}
            x2={width - pad}
            y2={height - (t - bbox.vMin + pad)}
            stroke="rgba(255,255,255,0.04)"
            stroke-width={width / 1500}
          />
          <text
            x={pad / 2}
            y={height - (t - bbox.vMin + pad) + 1}
            fill="rgba(184,184,192,0.8)"
            font-size={Math.max(2, height / 22)}
            text-anchor="end"
            font-family="JetBrains Mono, monospace"
          >
            {t.toFixed(0)}
          </text>
        </g>
      {/each}

      <!-- the cut curve -->
      <path
        d={path}
        fill="none"
        stroke={accent}
        stroke-width={Math.max(0.45, width / 700)}
        stroke-linejoin="round"
        stroke-linecap="round"
      />

      <!-- Endpoints markers -->
      {#if points.length > 0 && !closed}
        <circle
          cx={points[0].u - bbox.uMin + pad}
          cy={height - (points[0].v - bbox.vMin + pad)}
          r={Math.max(0.6, width / 600)}
          fill={accent}
        />
        <circle
          cx={points[points.length - 1].u - bbox.uMin + pad}
          cy={height - (points[points.length - 1].v - bbox.vMin + pad)}
          r={Math.max(0.6, width / 600)}
          fill={accent}
        />
      {/if}
    </svg>

    <!-- corner labels -->
    <div class="absolute left-2 top-2 text-[10px] text-ash uppercase tracking-[0.18em]">
      v (mm)
    </div>
    <div class="absolute right-2 bottom-2 text-[10px] text-ash uppercase tracking-[0.18em]">
      u (mm)
    </div>
  </div>
</section>
