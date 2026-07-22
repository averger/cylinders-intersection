<script lang="ts">
  import Segmented from "./Segmented.svelte";
  import Slider from "./Slider.svelte";
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
  import { downloadSVG } from "../lib/svg";
  import {
    planTiles,
    tileLabel,
    type PageFormat,
    type PageOrientation,
    type PatternKind,
    type ScaleMode,
  } from "../lib/export";

  const MARGIN = 16; // mm of breathing room around the drawing in the preview

  const formatOptions = [
    { label: "A4", value: "a4" as PageFormat },
    { label: "A3", value: "a3" as PageFormat },
    { label: "A2", value: "a2" as PageFormat },
  ];
  const orientationOptions = [
    { label: "Paysage", value: "landscape" as PageOrientation },
    { label: "Portrait", value: "portrait" as PageOrientation },
  ];
  const scaleOptions = [
    { label: "1:1 · tuiles", value: "one_to_one" as ScaleMode },
    { label: "Ajusté · 1 page", value: "fit" as ScaleMode },
  ];

  let patternOptions = $derived(
    editor.available.map((k) => ({
      label: k === "branch" ? "Tube incliné" : "Gueule de loup",
      value: k as PatternKind,
    })),
  );

  let points = $derived.by(() => {
    const r = store.result;
    if (!r) return [];
    return editor.active === "branch" ? r.dev_branch : (r.dev_main ?? []);
  });
  let closed = $derived(
    editor.active === "branch" ? true : (store.result?.dev_main_closed ?? false),
  );
  let box = $derived(editor.patternBox(editor.active));

  let viewW = $derived(box ? box.w + MARGIN * 2 : 100);
  let viewH = $derived(box ? box.h + MARGIN * 2 : 60);

  // mm → preview coordinates (SVG y grows downward).
  function X(u: number): number {
    return box ? u - box.uMin + MARGIN : 0;
  }
  function Y(v: number): number {
    return box ? box.h - (v - box.vMin) + MARGIN : 0;
  }

  let cutPath = $derived.by(() => {
    if (points.length === 0) return "";
    const d = points
      .map((p, i) => `${i === 0 ? "M" : "L"}${X(p.u).toFixed(3)} ${Y(p.v).toFixed(3)}`)
      .join(" ");
    return closed ? `${d} Z` : d;
  });

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

  let tiles = $derived.by(() => {
    if (!box || editor.scale !== "one_to_one") return [];
    const plan = planTiles(editor.page, box.w, box.h);
    const out: { x: number; y: number; w: number; h: number; label: string }[] = [];
    for (let row = 0; row < plan.rows; row++) {
      for (let col = 0; col < plan.cols; col++) {
        const u = box.uMin + col * plan.stepX;
        const vTop = box.vMin + box.h - row * plan.stepY;
        out.push({
          x: X(u),
          y: Y(vTop),
          w: plan.viewW,
          h: plan.viewH,
          label: tileLabel(col, row),
        });
      }
    }
    return out;
  });

  // Effective scale of the single-page "fit" mode (printable area over
  // drawing extents, never upscaled) — same rule as the backend renderer.
  let fitScale = $derived.by(() => {
    if (!box || editor.scale !== "fit") return 1;
    const plan = planTiles(editor.page, 0.0001, 0.0001);
    return Math.min(plan.viewW / box.w, plan.viewH / box.h, 1);
  });

  let annotationsHere = $derived(
    editor.annotations
      .map((a, index) => ({ a, index }))
      .filter(({ a }) => a.pattern === editor.active),
  );

  // ----- annotation dragging -------------------------------------------------
  let svgEl = $state<SVGSVGElement | undefined>(undefined);
  let dragIndex: number | null = null;

  function clientToMm(e: PointerEvent | MouseEvent): { u: number; v: number } | null {
    if (!svgEl || !box) return null;
    const pt = new DOMPoint(e.clientX, e.clientY);
    const ctm = svgEl.getScreenCTM();
    if (!ctm) return null;
    const local = pt.matrixTransform(ctm.inverse());
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
    dragIndex = null;
  }

  function onBackgroundDblClick(e: MouseEvent) {
    const mm = clientToMm(e);
    if (mm) editor.addAnnotation(Math.round(mm.u), Math.round(mm.v));
  }

  // ----- standalone SVG export (WYSIWYG, print colours, 1 mm = 1 unit) ------
  function exportEditedSVG() {
    if (!box || points.length === 0) return;
    const esc = (s: string) =>
      s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
    const parts: string[] = [];
    parts.push(
      `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${viewW.toFixed(2)}mm" height="${viewH.toFixed(2)}mm" viewBox="0 0 ${viewW.toFixed(2)} ${viewH.toFixed(2)}" shape-rendering="geometricPrecision">`,
      `<rect width="100%" height="100%" fill="#ffffff"/>`,
    );
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
        `<rect x="${X(box.uMin).toFixed(2)}" y="${Y(box.vMin + box.h).toFixed(2)}" width="${box.w.toFixed(2)}" height="${box.h.toFixed(2)}" fill="none" stroke="#777" stroke-width="0.15" stroke-dasharray="1.6 1.4"/>`,
      );
    }
    if (editor.layers.axis) {
      for (let q = 0; q <= 4; q++) {
        const u = box.uMin + (box.w * q) / 4;
        parts.push(
          `<line x1="${X(u).toFixed(2)}" y1="${Y(box.vMin).toFixed(2)}" x2="${X(u).toFixed(2)}" y2="${Y(box.vMin + box.h).toFixed(2)}" stroke="#555" stroke-width="0.12" stroke-dasharray="3 1.2"/>`,
        );
        if (editor.layers.labels)
          parts.push(
            `<text x="${(X(u) + 0.8).toFixed(2)}" y="${(Y(box.vMin) - 0.9).toFixed(2)}" font-size="2.4" fill="#555" font-family="Helvetica, Arial, sans-serif">${q * 90}°</text>`,
          );
      }
    }
    parts.push(
      `<path d="${cutPath}" fill="none" stroke="#000" stroke-width="${editor.cutWidth}" stroke-linejoin="round" stroke-linecap="round"/>`,
    );
    for (const { a } of annotationsHere) {
      parts.push(
        `<text x="${X(a.u).toFixed(2)}" y="${Y(a.v).toFixed(2)}" font-size="${a.size_mm}" fill="#111" font-family="Helvetica, Arial, sans-serif">${esc(a.text)}</text>`,
      );
    }
    parts.push(`</svg>`);
    const name = editor.active === "branch" ? "gabarit-tube" : "gueule-de-loup";
    downloadSVG(`${name}-edite.svg`, parts.join("\n"));
  }
</script>

<section id="studio" class="max-w-[1600px] mx-auto px-6 lg:px-10 pb-16 scroll-mt-24">
  <div class="flex flex-col lg:flex-row lg:items-end lg:justify-between gap-3 mb-6">
    <div>
      <span class="pill">atelier d'export</span>
      <h2 class="text-pearl text-2xl lg:text-3xl font-semibold tracking-tight mt-2">
        Éditez, annotez, exportez
      </h2>
      <p class="text-ash text-sm max-w-2xl mt-1">
        Composez la planche avant l'export&nbsp;: calques, annotations déplaçables à la
        souris (double-clic pour en ajouter), cartouche et mise en page. Le PDF est
        vectoriel et mm-exact&nbsp;— en 1:1 il est tuilé avec repères d'assemblage&nbsp;;
        le DXF part directement en CAO/CNC.
      </p>
    </div>
    <div class="text-right text-[11px] text-ash num leading-relaxed">
      {editor.pageCount()} page{editor.pageCount() > 1 ? "s" : ""}
      {editor.page.format.toUpperCase()}
      {editor.page.orientation === "landscape" ? "paysage" : "portrait"}<br />
      {editor.scale === "one_to_one"
        ? "échelle 1:1 — tuiles à recoller"
        : `aperçu ajusté ${fitScale < 1 ? `(1:${(1 / fitScale).toFixed(2)})` : "(1:1)"}`}
    </div>
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-[380px_1fr] gap-6 items-start">
    <!-- ============================ Settings ============================ -->
    <aside class="glass-strong card-hairline p-6 flex flex-col gap-6">
      {#if patternOptions.length > 1}
        <section class="space-y-3">
          <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Motif affiché</span>
          <Segmented
            options={patternOptions}
            value={editor.active}
            onchange={(v) => (editor.active = v)}
          />
        </section>
      {/if}

      <section class="space-y-3">
        <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Mise en page</span>
        <div class="flex flex-wrap gap-3">
          <Segmented
            options={formatOptions}
            value={editor.page.format}
            onchange={(v) => (editor.page = { ...editor.page, format: v })}
          />
          <Segmented
            options={orientationOptions}
            value={editor.page.orientation}
            onchange={(v) => (editor.page = { ...editor.page, orientation: v })}
          />
        </div>
        <Segmented options={scaleOptions} value={editor.scale} onchange={(v) => (editor.scale = v)} />
      </section>

      <section class="space-y-3">
        <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Calques</span>
        <div class="grid grid-cols-2 gap-x-4 gap-y-2 text-sm text-silver">
          {#each [["grid", "Grille 10 mm"], ["frame", "Emprise du tube"], ["axis", "Génératrices 90°"], ["labels", "Étiquettes"], ["scale_bar", "Règle 100 mm"]] as [key, label] (key)}
            <label class="flex items-center gap-2 cursor-pointer normal-case tracking-normal text-sm text-silver">
              <input
                type="checkbox"
                class="accent-[#ff5b1a] w-3.5 h-3.5"
                checked={editor.layers[key as keyof typeof editor.layers]}
                onchange={(e) =>
                  (editor.layers = {
                    ...editor.layers,
                    [key]: (e.currentTarget as HTMLInputElement).checked,
                  })}
              />
              {label}
            </label>
          {/each}
        </div>
      </section>

      <section class="space-y-3">
        <Slider
          label="Trait de coupe"
          unit="mm"
          min={0.1}
          max={1}
          step={0.05}
          value={editor.cutWidth}
          decimals={2}
          onchange={(v) => (editor.cutWidth = v)}
        />
      </section>

      <section class="space-y-3">
        <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Cartouche</span>
        <div class="grid grid-cols-2 gap-2">
          {#each [["title", "Titre"], ["project", "Projet"], ["author", "Auteur"], ["date", "Date"]] as [key, ph] (key)}
            <input
              type="text"
              placeholder={ph}
              class="bg-black/40 border border-mist/50 rounded-md px-2.5 py-1.5 text-sm text-pearl placeholder:text-ash/60 focus:outline-none focus:border-ember/70 transition-colors"
              value={editor.titleBlock[key as keyof typeof editor.titleBlock]}
              oninput={(e) =>
                (editor.titleBlock = {
                  ...editor.titleBlock,
                  [key]: (e.currentTarget as HTMLInputElement).value,
                })}
            />
          {/each}
        </div>
        <input
          type="text"
          placeholder="Notes (procédé, matière, consignes…)"
          class="w-full bg-black/40 border border-mist/50 rounded-md px-2.5 py-1.5 text-sm text-pearl placeholder:text-ash/60 focus:outline-none focus:border-ember/70 transition-colors"
          value={editor.titleBlock.notes}
          oninput={(e) =>
            (editor.titleBlock = {
              ...editor.titleBlock,
              notes: (e.currentTarget as HTMLInputElement).value,
            })}
        />
      </section>

      <section class="space-y-3">
        <div class="flex items-center justify-between">
          <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Annotations</span>
          <button
            class="text-[11px] text-ember hover:text-ember-soft transition-colors uppercase tracking-[0.14em]"
            onclick={() => {
              if (box) editor.addAnnotation(box.uMin + box.w / 2, box.vMin + box.h / 2);
            }}
          >
            + ajouter
          </button>
        </div>
        {#if annotationsHere.length === 0}
          <p class="text-[11px] text-ash/70">
            Aucune annotation sur ce motif — double-cliquez sur le plan pour en placer une.
          </p>
        {:else}
          <ul class="space-y-2">
            {#each annotationsHere as { a, index } (index)}
              <li class="flex items-center gap-2">
                <input
                  type="text"
                  class="flex-1 bg-black/40 border rounded-md px-2.5 py-1.5 text-sm text-pearl focus:outline-none transition-colors {editor.selected === index
                    ? 'border-ember/80'
                    : 'border-mist/50 focus:border-ember/70'}"
                  value={a.text}
                  onfocus={() => (editor.selected = index)}
                  oninput={(e) => (editor.annotations[index].text = (e.currentTarget as HTMLInputElement).value)}
                />
                <span class="num text-[10px] text-ash w-24 text-right"
                  >{a.u.toFixed(0)} ; {a.v.toFixed(0)} mm</span
                >
                <button
                  class="text-ash hover:text-ember transition-colors text-base leading-none px-1"
                  aria-label="Supprimer l'annotation"
                  onclick={() => editor.removeAnnotation(index)}
                >
                  ×
                </button>
              </li>
            {/each}
          </ul>
        {/if}
      </section>

      <section class="pt-2 border-t border-white/5 space-y-3">
        <div class="flex flex-wrap gap-2">
          <button
            class="btn-primary flex-1 min-w-[120px]"
            disabled={editor.busy !== null || !store.result}
            onclick={() => editor.exportPdf()}
          >
            {editor.busy === "pdf" ? "génération…" : "Exporter PDF"}
          </button>
          <button
            class="btn-ghost flex-1 min-w-[100px] justify-center"
            disabled={editor.busy !== null || !store.result}
            onclick={() => editor.exportDxf()}
          >
            {editor.busy === "dxf" ? "génération…" : "DXF · CAO"}
          </button>
          <button
            class="btn-ghost flex-1 min-w-[100px] justify-center"
            disabled={!store.result}
            onclick={exportEditedSVG}
          >
            SVG 1:1
          </button>
        </div>
        {#if editor.error}
          <p class="text-[12px] text-ember">{editor.error}</p>
        {/if}
        <p class="text-[10px] text-ash/70 leading-relaxed">
          PDF&nbsp;: toutes les planches sélectionnées ({editor.available.length} motif{editor
            .available.length > 1
            ? "s"
            : ""}), cartouche et repères de collage inclus. DXF&nbsp;R12&nbsp;: calques
          CUT / FRAME / AXIS / TEXT / ANNOT, unités mm.
        </p>
      </section>
    </aside>

    <!-- ============================ Preview ============================ -->
    <div
      class="glass-strong card-hairline rounded-[var(--radius-card)] overflow-hidden relative min-h-[420px] lg:min-h-[560px] grid-bg"
    >
      {#if box && points.length > 0}
        <svg
          bind:this={svgEl}
          viewBox="0 0 {viewW} {viewH}"
          preserveAspectRatio="xMidYMid meet"
          class="w-full h-full p-4 select-none"
          role="application"
          aria-label="Aperçu du gabarit — double-clic pour annoter"
          ondblclick={onBackgroundDblClick}
          onpointermove={moveDrag}
          onpointerup={endDrag}
        >
          <!-- grid -->
          {#if editor.layers.grid}
            {#each uGrid as g (g.at)}
              <line
                x1={X(g.at)}
                y1={Y(box.vMin)}
                x2={X(g.at)}
                y2={Y(box.vMin + box.h)}
                stroke="rgba(255,255,255,{g.major ? 0.1 : 0.045})"
                stroke-width={g.major ? viewW / 1600 : viewW / 2600}
              />
            {/each}
            {#each vGrid as g (g.at)}
              <line
                x1={X(box.uMin)}
                y1={Y(g.at)}
                x2={X(box.uMin + box.w)}
                y2={Y(g.at)}
                stroke="rgba(255,255,255,{g.major ? 0.1 : 0.045})"
                stroke-width={g.major ? viewW / 1600 : viewW / 2600}
              />
            {/each}
          {/if}

          <!-- unwrapped tube frame -->
          {#if editor.layers.frame}
            <rect
              x={X(box.uMin)}
              y={Y(box.vMin + box.h)}
              width={box.w}
              height={box.h}
              fill="rgba(255,255,255,0.015)"
              stroke="rgba(255,255,255,0.22)"
              stroke-width={viewW / 900}
              stroke-dasharray="{viewW / 180} {viewW / 220}"
            />
          {/if}

          <!-- quarter generators -->
          {#if editor.layers.axis}
            {#each [0, 1, 2, 3, 4] as q (q)}
              <line
                x1={X(box.uMin + (box.w * q) / 4)}
                y1={Y(box.vMin)}
                x2={X(box.uMin + (box.w * q) / 4)}
                y2={Y(box.vMin + box.h)}
                stroke="rgba(41,194,255,0.28)"
                stroke-width={viewW / 1400}
                stroke-dasharray="{viewW / 140} {viewW / 300}"
              />
              {#if editor.layers.labels}
                <text
                  x={X(box.uMin + (box.w * q) / 4) + viewW / 300}
                  y={Y(box.vMin) - viewH / 90}
                  font-size={Math.max(2, viewW / 90)}
                  fill="rgba(136,223,250,0.7)"
                  font-family="JetBrains Mono, monospace"
                >
                  {q * 90}°
                </text>
              {/if}
            {/each}
          {/if}

          <!-- page tiling overlay (1:1) -->
          {#each tiles as t (t.label)}
            <rect
              x={t.x}
              y={t.y}
              width={t.w}
              height={t.h}
              fill="none"
              stroke="rgba(255,91,26,0.35)"
              stroke-width={viewW / 1100}
            />
            <text
              x={t.x + viewW / 250}
              y={t.y + Math.max(2.6, viewW / 70)}
              font-size={Math.max(2.6, viewW / 80)}
              fill="rgba(255,91,26,0.6)"
              font-family="JetBrains Mono, monospace"
            >
              {t.label}
            </text>
          {/each}

          <!-- cut line -->
          <path
            d={cutPath}
            fill="none"
            stroke="#ff5b1a"
            stroke-width={Math.max(editor.cutWidth, viewW / 700)}
            stroke-linejoin="round"
            stroke-linecap="round"
          />

          <!-- annotations -->
          {#each annotationsHere as { a, index } (index)}
            <!-- svelte-ignore a11y_no_static_element_interactions -->
            <text
              x={X(a.u)}
              y={Y(a.v)}
              font-size={Math.max(a.size_mm, viewW / 110)}
              fill={editor.selected === index ? "#ffffff" : "rgba(245,245,247,0.85)"}
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
                stroke="#ff5b1a"
                stroke-width={viewW / 1400}
              />
            {/if}
          {/each}
        </svg>

        <div class="absolute left-4 top-4 pill pointer-events-none">
          {editor.active === "branch" ? "gabarit tube" : "gueule de loup"}
        </div>
        <div
          class="absolute right-4 bottom-3 text-[10px] text-ash/70 uppercase tracking-[0.18em] pointer-events-none"
        >
          double-clic&nbsp;: annoter · glisser&nbsp;: déplacer
        </div>
      {:else}
        <div class="absolute inset-0 grid place-items-center text-ash text-sm">
          Aucun motif à afficher — ajustez les paramètres.
        </div>
      {/if}
    </div>
  </div>
</section>
