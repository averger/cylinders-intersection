<script lang="ts">
  // Annotations — shown in the sidebar when the 2D workspace is active.
  // Layout options live in the header Options menu, exports in Exporter.
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
</script>

<div class="flex flex-col gap-5">
  {#if editor.annotations.length > 0}
    <section class="space-y-2">
      <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Annotations</span>
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
    </section>
  {/if}

  {#if editor.error}
    <p class="text-[12px] text-ember">{editor.error}</p>
  {/if}
</div>
