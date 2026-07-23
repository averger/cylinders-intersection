<script lang="ts">
  // Annotations + export actions — shown in the sidebar when the 2D
  // workspace is active. Layout options live in the header Options menu.
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
</script>

<div class="flex flex-col gap-5">
  <section class="space-y-2">
    <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Annotations</span>
    {#if editor.annotations.length === 0}
      <p class="text-[11px] text-ash/70">
        Double-cliquez sur une mise à plat pour placer un repère texte.
      </p>
    {:else}
      <ul class="space-y-1.5">
        {#each editor.annotations as a, index (index)}
          <li class="flex items-center gap-1.5">
            <span
              class="shrink-0 w-1.5 h-1.5 rounded-full"
              style="background: {a.pattern === 'branch' ? 'var(--ember)' : 'var(--cyan)'}"
              title={a.pattern === "branch" ? "tube incliné" : "gueule de loup"}
            ></span>
            <input
              type="text"
              class="flex-1 min-w-0 panel-input px-2 py-1 text-[12px] {editor.selected ===
              index
                ? 'border-ember/80'
                : 'border-mist/50 focus:border-ember/70'}"
              value={a.text}
              onfocus={() => (editor.selected = index)}
              oninput={(e) => {
                editor.annotations[index].text = (e.currentTarget as HTMLInputElement).value;
                store.persist();
              }}
            />
            <span class="num text-[9px] text-ash shrink-0">{a.u.toFixed(0)};{a.v.toFixed(0)}</span>
            <button
              class="text-ash hover:text-ember transition-colors text-sm leading-none px-0.5"
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

  <section class="pt-3 hl-t space-y-2">
    <div class="flex gap-2">
      <button
        class="btn-primary flex-1 !py-2.5 text-sm"
        disabled={editor.busy !== null || !store.result}
        onclick={() => editor.exportPdf()}
      >
        {editor.busy === "pdf" ? "génération…" : "PDF 1:1"}
      </button>
      <button
        class="btn-ghost flex-1 justify-center !py-2.5 text-sm"
        disabled={editor.busy !== null || !store.result}
        onclick={() => editor.exportDxf()}
      >
        {editor.busy === "dxf" ? "génération…" : "DXF"}
      </button>
    </div>
    {#if editor.error}
      <p class="text-[12px] text-ember">{editor.error}</p>
    {/if}
    <p class="text-[10px] text-ash/70 leading-relaxed">
      Le PDF contient toutes les mises à plat ({editor.available.length}) avec cartouche et
      règle de contrôle. DXF R12 en mm — calques CUT / FRAME / AXIS / TEXT / ANNOT.
    </p>
  </section>
</div>
