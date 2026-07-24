<script lang="ts">
  // 2D workspace: BOTH developed patterns stacked (the intersection unrolled
  // on each cylinder), sharing the editor's layers / page setup.
  import PatternCanvas from "./PatternCanvas.svelte";
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
</script>

<div class="absolute inset-0 flex flex-col gap-3 p-3 overflow-y-auto">
  {#if store.params.mode === "multi"}
    {#if store.multiResult}
      {@const m = store.multiResult}
      <PatternCanvas
        kind="main"
        data={{
          points: m.holes[0]?.pts ?? [],
          closed: m.holes[0]?.closed ?? false,
          holes: m.holes.slice(1).map((h) => ({ pts: h.pts, closed: h.closed })),
          title: "Tube principal — lumières",
          diameter: m.r1 * 2,
          circumference: m.circumference_main,
          accent: "var(--cyan)",
          annotable: true,
          slugName: "cylix-tube-principal",
        }}
      />
      {#each m.branches as b, i (i)}
        <PatternCanvas
          kind="branch"
          data={{
            points: b.dev,
            closed: false,
            holes: b.holes.map((h) => ({ pts: h.pts, closed: h.closed })),
            title: `Piquage ${i + 1} — Ø ${(b.r * 2).toFixed(0)}${b.cut_by_neighbor ? " · selle sur piquage prioritaire" : ""}`,
            diameter: b.r * 2,
            circumference: b.circumference,
            accent: "var(--ember)",
            annotable: false,
            slugName: `cylix-piquage-${i + 1}`,
          }}
        />
      {/each}
      {#each m.warnings as w, wi (wi)}
        <p class="shrink-0 text-[11px] text-ember px-1">⚠ {w}</p>
      {/each}
    {/if}
  {:else}
    <PatternCanvas kind="branch" />
    {#if store.result?.dev_main}
      <PatternCanvas kind="main" />
    {/if}
  {/if}

  <div
    class="shrink-0 flex flex-wrap items-center justify-between gap-2 text-[10px] text-ash uppercase tracking-[0.16em] px-1 pb-1"
  >
    <span>double-clic : annoter · glisser : déplacer l'annotation</span>
    <span class="num normal-case tracking-normal">
      export : {editor.pageCount()} page{editor.pageCount() > 1 ? "s" : ""}
      {editor.page.format.toUpperCase()}
      {editor.page.orientation === "landscape" ? "paysage" : "portrait"} ·
      {editor.scale === "one_to_one" ? "échelle 1:1 tuilée" : "ajusté 1 page"}
    </span>
  </div>
</div>
