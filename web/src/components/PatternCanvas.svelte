<script lang="ts">
  // One developed pattern (mm coordinates) with its layers, 1:1 tile
  // overlay and draggable annotations.  Double-click adds an annotation.
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
  import { downloadSVG } from "../lib/svg";
  import type { PatternKind } from "../lib/export";
  import type { DevPoint } from "../lib/api";

  /** Explicit pattern data — used by the multi-branch mode where the sheets
   * do not map 1-to-1 onto the historical branch/main pair. */
  export interface PatternData {
    points: DevPoint[];
    closed: boolean;
    holes?: { pts: DevPoint[]; closed: boolean }[];
    title: string;
    diameter: number;
    circumference: number;
    accent?: string;
    /** Annotations are only editable on the sheet that owns `kind`. */
    annotable?: boolean;
    slugName?: string;
  }

  interface Props {
    kind: PatternKind;
    data?: PatternData;
  }
  let { kind, data }: Props = $props();

  const MARGIN = 16; // mm of breathing room around the drawing

  let points = $derived.by(() => {
    if (data) return data.points;
    const r = store.result;
    if (!r) return [];
    return kind === "branch" ? r.dev_branch : (r.dev_main ?? []);
  });
  // The branch development is always an open curve (u = 0 meets u = 2πR on
  // the rolled tube); the gueule de loup closes when the backend says so.
  let closed = $derived(
    data ? data.closed : kind === "branch" ? false : (store.result?.dev_main_closed ?? false),
  );
  let holes = $derived(data?.holes ?? []);
  let annotable = $derived(data ? data.annotable === true : true);
  // Viewport framed on the cut extents (+ margin) — the gueule de loup
  // sheet no longer spans the whole unwrapped circumference.
  let box = $derived.by(() => {
    if (data) {
      const all = [...data.points, ...holes.flatMap((h) => h.pts)];
      if (all.length === 0) return null;
      let uMin = all[0].u,
        uMax = uMin,
        vMin = all[0].v,
        vMax = vMin;
      for (const p of all) {
        uMin = Math.min(uMin, p.u);
        uMax = Math.max(uMax, p.u);
        vMin = Math.min(vMin, p.v);
        vMax = Math.max(vMax, p.v);
      }
      if (annotable) {
        for (const a of editor.annotations) {
          if (a.pattern !== kind) continue;
          uMin = Math.min(uMin, a.u);
          uMax = Math.max(uMax, a.u);
          vMin = Math.min(vMin, a.v);
          vMax = Math.max(vMax, a.v);
        }
      }
      return { uMin, vMin, w: Math.max(uMax - uMin, 1), h: Math.max(vMax - vMin, 0.1) };
    }
    const cb = editor.cutBox(kind);
    if (!cb) return null;
    return { uMin: cb.uMin, vMin: cb.vMin, w: Math.max(cb.uMax - cb.uMin, 1), h: Math.max(cb.vMax - cb.vMin, 0.1) };
  });

  let accent = $derived(data?.accent ?? (kind === "branch" ? "var(--ember)" : "var(--cyan)"));
  let title = $derived(
    data
      ? data.title
      : kind === "branch"
        ? store.result?.mode === "cyl_cyl"
          ? "Développé tube incliné — Ø₂"
          : "Développé du tube — coupe plane"
        : "Gueule de loup — Ø₁",
  );
  let diameter = $derived.by(() => {
    if (data) return data.diameter;
    const r = store.result;
    if (!r) return 0;
    return kind === "branch" ? (r.r2 ?? r.r1) * 2 : r.r1 * 2;
  });
  let circumference = $derived.by(() => {
    if (data) return data.circumference;
    const r = store.result;
    if (!r) return 0;
    return kind === "branch" ? (r.circumference_branch ?? r.circumference_main) : r.circumference_main;
  });

  // Unwrapped-period start containing the curve (u may live in any
  // 2πR-period after the angle unwrap).
  let frameStart = $derived.by(() => {
    if (!box || circumference <= 0) return 0;
    const center = box.uMin + box.w / 2;
    return Math.floor(center / circumference) * circumference;
  });

  // True tube generatrices (multiples of 90°, u = 0 ⇔ θ = 0) falling
  // inside the cropped viewport.
  let gens = $derived.by(() => {
    if (!box || circumference <= 0) return [];
    const step = circumference / 4;
    const out: { u: number; deg: number; ref: boolean }[] = [];
    const k0 = Math.floor((box.uMin - 1e-9) / step);
    for (let k = k0; k * step <= box.uMin + box.w + 1e-9; k++) {
      const u = k * step;
      if (u < box.uMin - 1e-9) continue;
      // k multiple of 4 ⇔ the θ = 0 generatrix — THE wrap-alignment datum,
      // essential when the pattern is phase-shifted (oriented plane).
      out.push({ u, deg: ((k * 90) % 360 + 360) % 360, ref: k % 4 === 0 });
    }
    return out;
  });

  // Alignment marks: where each generatrix crosses the cut lines (ticks to
  // match with lines traced on the tube) + the axis-plane datum v = 0.
  function crossingsOf(pts: DevPoint[], isClosed: boolean, g: number, out: number[]) {
    const n = pts.length;
    if (n < 2) return;
    const last = isClosed ? n : n - 1;
    for (let i = 0; i < last; i++) {
      const a = pts[i];
      const b = pts[(i + 1) % n];
      if ((a.u - g) * (b.u - g) < 0) {
        const t = (g - a.u) / (b.u - a.u);
        out.push(a.v + t * (b.v - a.v));
      }
    }
  }
  function crossings(g: number): number[] {
    const out: number[] = [];
    crossingsOf(points, closed, g, out);
    for (const h of holes) crossingsOf(h.pts, h.closed, g, out);
    return out;
  }
  let ticks = $derived(gens.flatMap((g) => crossings(g.u).map((v) => ({ u: g.u, v }))));
  let showDatum = $derived(!!box && box.vMin < 0 && box.vMin + box.h > 0);

  let viewW = $derived(box ? box.w + MARGIN * 2 : 100);
  let viewH = $derived(box ? box.h + MARGIN * 2 : 60);

  function X(u: number): number {
    return box ? u - box.uMin + MARGIN : 0;
  }
  function Y(v: number): number {
    return box ? box.h - (v - box.vMin) + MARGIN : 0;
  }

  function pathOf(pts: DevPoint[], isClosed: boolean): string {
    if (pts.length === 0) return "";
    const d = pts
      .map((p, i) => `${i === 0 ? "M" : "L"}${X(p.u).toFixed(3)} ${Y(p.v).toFixed(3)}`)
      .join(" ");
    return isClosed ? `${d} Z` : d;
  }
  let cutPath = $derived(pathOf(points, closed));
  let holePaths = $derived(holes.map((h) => pathOf(h.pts, h.closed)));

  function gridLines(min: number, span: number): { at: number; major: boolean }[] {
    const out: { at: number; major: boolean }[] = [];
    const start = Math.floor(min / 10) * 10;
    for (let g = start; g <= min + span + 1e-9; g += 10) {
      if (g < min - 1e-9) continue;
      out.push({ at: g, major: Math.abs(g / 50 - Math.round(g / 50)) < 1e-9 });
    }
    return out;
  }
  let uGrid = $derived(box ? gridLines(box.uMin, box.w) : []);
  let vGrid = $derived(box ? gridLines(box.vMin, box.h) : []);

  // Non-empty pages only, tiled over the cut extents — mirrors the backend.
  let tiles = $derived.by(() => {
    const raw = data
      ? editor.tilesFor(
          [...points, ...holes.flatMap((h) => h.pts)],
          annotable ? editor.annotations.filter((a) => a.pattern === kind) : [],
        )
      : editor.tiles(kind);
    return raw.map((t) => ({ x: X(t.u0), y: Y(t.vTop), w: t.w, h: t.h, label: t.label }));
  });

  let annotationsHere = $derived(
    annotable
      ? editor.annotations.map((a, index) => ({ a, index })).filter(({ a }) => a.pattern === kind)
      : [],
  );

  // ----- annotation dragging ----------------------------------------------
  let svgEl = $state<SVGSVGElement | undefined>(undefined);
  let dragIndex: number | null = null;

  function clientToMm(e: PointerEvent | MouseEvent): { u: number; v: number } | null {
    if (!svgEl || !box) return null;
    const ctm = svgEl.getScreenCTM();
    if (!ctm) return null;
    const local = new DOMPoint(e.clientX, e.clientY).matrixTransform(ctm.inverse());
    return { u: local.x - MARGIN + box.uMin, v: box.h - (local.y - MARGIN) + box.vMin };
  }

  function startDrag(e: PointerEvent, index: number) {
    e.stopPropagation();
    dragIndex = index;
    editor.selected = index;
    (e.currentTarget as Element).setPointerCapture(e.pointerId);
  }
  function moveDrag(e: PointerEvent) {
    if (dragIndex === null) return;
    const mm = clientToMm(e);
    if (!mm) return;
    editor.annotations[dragIndex].u = Math.round(mm.u * 10) / 10;
    editor.annotations[dragIndex].v = Math.round(mm.v * 10) / 10;
  }
  function endDrag() {
    if (dragIndex !== null) store.persist();
    dragIndex = null;
  }

  function onBackgroundDblClick(e: MouseEvent) {
    if (!annotable) return;
    const mm = clientToMm(e);
    if (mm) editor.addAnnotation(kind, Math.round(mm.u), Math.round(mm.v));
  }

  // ----- standalone SVG export (print colours, 1 mm = 1 unit) -------------
  export function exportSVG() {
    if (!box || points.length === 0) return;
    const esc = (s: string) =>
      s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
    const parts: string[] = [
      `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${viewW.toFixed(2)}mm" height="${viewH.toFixed(2)}mm" viewBox="0 0 ${viewW.toFixed(2)} ${viewH.toFixed(2)}" shape-rendering="geometricPrecision">`,
      `<rect width="100%" height="100%" fill="#ffffff"/>`,
    ];
    if (editor.layers.grid) {
      for (const g of uGrid)
        parts.push(
          `<line x1="${X(g.at).toFixed(2)}" y1="${Y(box.vMin).toFixed(2)}" x2="${X(g.at).toFixed(2)}" y2="${Y(box.vMin + box.h).toFixed(2)}" stroke="${g.major ? "#b9b9b9" : "#e2e2e2"}" stroke-width="${g.major ? 0.09 : 0.05}"/>`,
        );
      for (const g of vGrid)
        parts.push(
          `<line x1="${X(box.uMin).toFixed(2)}" y1="${Y(g.at).toFixed(2)}" x2="${X(box.uMin + box.w).toFixed(2)}" y2="${Y(g.at).toFixed(2)}" stroke="${g.major ? "#b9b9b9" : "#e2e2e2"}" stroke-width="${g.major ? 0.09 : 0.05}"/>`,
        );
    }
    if (editor.layers.frame) {
      parts.push(
        `<rect x="${X(frameStart).toFixed(2)}" y="${Y(box.vMin + box.h).toFixed(2)}" width="${circumference.toFixed(2)}" height="${box.h.toFixed(2)}" fill="none" stroke="#777" stroke-width="0.15" stroke-dasharray="1.6 1.4"/>`,
      );
    }
    if (editor.layers.axis) {
      for (const g of gens) {
        parts.push(
          g.ref
            ? `<line x1="${X(g.u).toFixed(2)}" y1="${Y(box.vMin).toFixed(2)}" x2="${X(g.u).toFixed(2)}" y2="${Y(box.vMin + box.h).toFixed(2)}" stroke="#cc4a0c" stroke-width="0.3"/>`
            : `<line x1="${X(g.u).toFixed(2)}" y1="${Y(box.vMin).toFixed(2)}" x2="${X(g.u).toFixed(2)}" y2="${Y(box.vMin + box.h).toFixed(2)}" stroke="#555" stroke-width="0.12" stroke-dasharray="3 1.2"/>`,
        );
        if (editor.layers.labels)
          parts.push(
            `<text x="${(X(g.u) + 0.8).toFixed(2)}" y="${(Y(box.vMin) - 0.9).toFixed(2)}" font-size="2.4" fill="${g.ref ? "#cc4a0c" : "#555"}" font-family="Helvetica, Arial, sans-serif">${g.ref ? `${g.deg}° réf` : `${g.deg}°`}</text>`,
          );
      }
      if (showDatum) {
        parts.push(
          `<line x1="${X(box.uMin).toFixed(2)}" y1="${Y(0).toFixed(2)}" x2="${X(box.uMin + box.w).toFixed(2)}" y2="${Y(0).toFixed(2)}" stroke="#333" stroke-width="0.1" stroke-dasharray="5 1.6"/>`,
          `<text x="${(X(box.uMin) + 1).toFixed(2)}" y="${(Y(0) - 0.7).toFixed(2)}" font-size="2" fill="#333" font-family="Helvetica, Arial, sans-serif">réf. plan des axes</text>`,
        );
      }
      for (const t of ticks) {
        parts.push(
          `<line x1="${(X(t.u) - 2.5).toFixed(2)}" y1="${Y(t.v).toFixed(2)}" x2="${(X(t.u) + 2.5).toFixed(2)}" y2="${Y(t.v).toFixed(2)}" stroke="#000" stroke-width="0.3"/>`,
        );
      }
    }
    parts.push(
      `<path d="${cutPath}" fill="none" stroke="#000" stroke-width="${editor.cutWidth}" stroke-linejoin="round" stroke-linecap="round"/>`,
    );
    for (const hp of holePaths) {
      parts.push(
        `<path d="${hp}" fill="none" stroke="#000" stroke-width="${editor.cutWidth}" stroke-linejoin="round" stroke-linecap="round"/>`,
      );
    }
    for (const { a } of annotationsHere) {
      parts.push(
        `<text x="${X(a.u).toFixed(2)}" y="${Y(a.v).toFixed(2)}" font-size="${a.size_mm}" fill="#111" font-family="Helvetica, Arial, sans-serif">${esc(a.text)}</text>`,
      );
    }
    parts.push(`</svg>`);
    const name =
      data?.slugName ?? (kind === "branch" ? "cylix-tube" : "cylix-gueule-de-loup");
    downloadSVG(`${name}.svg`, parts.join("\n"));
  }
</script>

<section class="glass card-hairline relative overflow-hidden grid-bg flex flex-col min-h-[260px]">
  <header class="flex items-center justify-between gap-3 px-4 pt-3 pb-1 z-10">
    <div class="flex items-center gap-3 min-w-0">
      <span class="pill shrink-0" style="border-color: color-mix(in srgb, {accent} 40%, transparent); color: {accent}">{title}</span>
      {#if points.length > 0}
        <span class="num text-[10px] text-ash truncate">
          Ø {diameter.toFixed(1)} mm · largeur {box ? box.w.toFixed(1) : "—"} mm · hauteur {box
            ? box.h.toFixed(1)
            : "—"} mm
        </span>
      {/if}
    </div>
    <button
      class="text-[10px] uppercase tracking-[0.16em] text-ash hover:text-pearl transition-colors shrink-0"
      onclick={() => exportSVG()}
    >
      svg 1:1 ↓
    </button>
  </header>

  {#if box && points.length > 0}
    <svg
      bind:this={svgEl}
      viewBox="0 0 {viewW} {viewH}"
      preserveAspectRatio="xMidYMid meet"
      class="w-full flex-1 min-h-0 px-3 pb-3 select-none"
      role="application"
      aria-label="{title} — double-clic pour annoter"
      ondblclick={onBackgroundDblClick}
      onpointermove={moveDrag}
      onpointerup={endDrag}
    >
      {#if editor.layers.grid}
        {#each uGrid as g (g.at)}
          <line
            x1={X(g.at)}
            y1={Y(box.vMin)}
            x2={X(g.at)}
            y2={Y(box.vMin + box.h)}
            style="stroke: {g.major ? 'var(--grid-major)' : 'var(--grid-line)'}"
            stroke-width={g.major ? viewW / 1600 : viewW / 2600}
          />
        {/each}
        {#each vGrid as g (g.at)}
          <line
            x1={X(box.uMin)}
            y1={Y(g.at)}
            x2={X(box.uMin + box.w)}
            y2={Y(g.at)}
            style="stroke: {g.major ? 'var(--grid-major)' : 'var(--grid-line)'}"
            stroke-width={g.major ? viewW / 1600 : viewW / 2600}
          />
        {/each}
      {/if}

      {#if editor.layers.frame}
        <rect
          x={X(frameStart)}
          y={Y(box.vMin + box.h)}
          width={circumference}
          height={box.h}
          style="fill: var(--frame-fill); stroke: var(--frame-line)"
          stroke-width={viewW / 900}
          stroke-dasharray="{viewW / 180} {viewW / 220}"
        />
      {/if}

      {#if editor.layers.axis}
        {#each gens as g (g.u)}
          <line
            x1={X(g.u)}
            y1={Y(box.vMin)}
            x2={X(g.u)}
            y2={Y(box.vMin + box.h)}
            style="stroke: {g.ref
              ? 'var(--ember)'
              : 'color-mix(in srgb, var(--cyan) 30%, transparent)'}"
            stroke-width={g.ref ? viewW / 600 : viewW / 1400}
            stroke-dasharray={g.ref ? "none" : `${viewW / 140} ${viewW / 300}`}
          />
          {#if editor.layers.labels}
            <text
              x={X(g.u) + viewW / 300}
              y={Y(box.vMin) - viewH / 90}
              font-size={Math.max(2, viewW / 90)}
              style="fill: {g.ref
                ? 'var(--ember)'
                : 'color-mix(in srgb, var(--cyan) 70%, transparent)'}"
              font-weight={g.ref ? "600" : "400"}
              font-family="JetBrains Mono, monospace"
            >
              {g.ref ? `${g.deg}° · réf` : `${g.deg}°`}
            </text>
          {/if}
        {/each}

        {#if showDatum && box}
          <line
            x1={X(box.uMin)}
            y1={Y(0)}
            x2={X(box.uMin + box.w)}
            y2={Y(0)}
            style="stroke: color-mix(in srgb, var(--text) 40%, transparent)"
            stroke-width={viewW / 1600}
            stroke-dasharray="{viewW / 110} {viewW / 260}"
          />
          {#if editor.layers.labels}
            <text
              x={X(box.uMin) + viewW / 250}
              y={Y(0) - viewH / 120}
              font-size={Math.max(1.8, viewW / 110)}
              style="fill: color-mix(in srgb, var(--text) 55%, transparent)"
              font-family="JetBrains Mono, monospace"
            >
              réf. plan des axes
            </text>
          {/if}
        {/if}

        {#each ticks as t (`${t.u}-${t.v}`)}
          <line
            x1={X(t.u) - viewW / 130}
            y1={Y(t.v)}
            x2={X(t.u) + viewW / 130}
            y2={Y(t.v)}
            style="stroke: var(--text)"
            stroke-width={viewW / 700}
          />
        {/each}
      {/if}

      {#each tiles as t (t.label)}
        <rect
          x={t.x}
          y={t.y}
          width={t.w}
          height={t.h}
          fill="none"
          style="stroke: color-mix(in srgb, var(--ember) 35%, transparent)"
          stroke-width={viewW / 1100}
        />
        <text
          x={t.x + viewW / 250}
          y={t.y + Math.max(2.6, viewW / 70)}
          font-size={Math.max(2.6, viewW / 80)}
          style="fill: color-mix(in srgb, var(--ember) 60%, transparent)"
          font-family="JetBrains Mono, monospace"
        >
          {t.label}
        </text>
      {/each}

      <path
        d={cutPath}
        fill="none"
        style="stroke: {accent}"
        stroke-width={Math.max(editor.cutWidth, viewW / 700)}
        stroke-linejoin="round"
        stroke-linecap="round"
      />
      {#each holePaths as hp, hi (hi)}
        <path
          d={hp}
          fill="none"
          style="stroke: {accent}"
          stroke-width={Math.max(editor.cutWidth, viewW / 700)}
          stroke-linejoin="round"
          stroke-linecap="round"
        />
      {/each}

      {#if !closed && points.length > 0}
        <circle cx={X(points[0].u)} cy={Y(points[0].v)} r={viewW / 400} style="fill: {accent}" />
        <circle
          cx={X(points[points.length - 1].u)}
          cy={Y(points[points.length - 1].v)}
          r={viewW / 400}
          style="fill: {accent}"
        />
      {/if}

      {#each annotationsHere as { a, index } (index)}
        <!-- svelte-ignore a11y_no_static_element_interactions -->
        <text
          x={X(a.u)}
          y={Y(a.v)}
          font-size={Math.max(a.size_mm, viewW / 110)}
          style="fill: {editor.selected === index ? 'var(--text)' : 'color-mix(in srgb, var(--text) 85%, transparent)'}"
          font-family="Inter, sans-serif"
          class="cursor-move"
          onpointerdown={(e) => startDrag(e, index)}
        >
          {a.text}
        </text>
        {#if editor.selected === index}
          <circle
            cx={X(a.u)}
            cy={Y(a.v)}
            r={viewW / 260}
            fill="none"
            style="stroke: var(--ember)"
            stroke-width={viewW / 1400}
          />
        {/if}
      {/each}
    </svg>
  {:else}
    <div class="flex-1 grid place-items-center text-ash text-sm">Aucune courbe à afficher.</div>
  {/if}
</section>
