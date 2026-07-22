<script lang="ts">
  // 2D workspace: BOTH developed patterns stacked (the intersection unrolled
  // on each cylinder), sharing the editor's layers / page setup.
  import PatternCanvas from "./PatternCanvas.svelte";
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
</script>

<div class="absolute inset-0 flex flex-col gap-3 p-3 overflow-y-auto">
  <PatternCanvas kind="branch" />
  {#if store.result?.dev_main}
    <PatternCanvas kind="main" />
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
