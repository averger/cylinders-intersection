<script lang="ts">
  // Study tab bar — click to activate, double-click to rename, × to close,
  // + for a new study (parameters copied from the current one).
  import { store } from "../lib/store.svelte";

  let renamingId = $state<string | null>(null);
  let renameValue = $state("");
  let renameInput = $state<HTMLInputElement | null>(null);

  function startRename(id: string, name: string) {
    renamingId = id;
    renameValue = name;
  }
  function commitRename() {
    if (renamingId) store.renameStudy(renamingId, renameValue);
    renamingId = null;
  }
  function onRenameKey(e: KeyboardEvent) {
    if (e.key === "Enter") {
      e.preventDefault();
      commitRename();
    } else if (e.key === "Escape") renamingId = null;
  }
  $effect(() => {
    if (renameInput) {
      renameInput.focus();
      renameInput.select();
    }
  });
</script>

<div
  class="no-print flex items-stretch gap-0.5 h-9 px-4 lg:px-6 border-b border-white/5 bg-black/30 overflow-x-auto"
  role="tablist"
  aria-label="Études"
>
  {#each store.studies as s (s.id)}
    {@const active = s.id === store.activeId}
    <div
      class="inline-flex items-center border-t-2 {active
        ? 'border-t-ember bg-white/[0.04]'
        : 'border-t-transparent'} border-r border-white/5"
    >
      {#if renamingId === s.id}
        <input
          class="w-32 h-6 mx-1.5 px-1.5 text-xs bg-black/60 border border-ember/70 rounded text-pearl focus:outline-none"
          bind:this={renameInput}
          bind:value={renameValue}
          onblur={commitRename}
          onkeydown={onRenameKey}
        />
      {:else}
        <button
          role="tab"
          aria-selected={active}
          class="max-w-[180px] h-full pl-3 pr-1.5 text-xs font-medium truncate transition-colors {active
            ? 'text-ember font-semibold'
            : 'text-ash hover:text-pearl'}"
          onclick={() => store.switchStudy(s.id)}
          ondblclick={() => startRename(s.id, s.name)}
          title="{s.name} — double-clic pour renommer"
        >
          {s.name}
        </button>
        {#if store.studies.length > 1}
          <button
            class="h-full px-1.5 text-sm leading-none text-ash/70 hover:text-ember transition-colors"
            aria-label="Fermer l'étude"
            onclick={(e) => {
              e.stopPropagation();
              store.closeStudy(s.id);
            }}
          >
            ×
          </button>
        {/if}
      {/if}
    </div>
  {/each}
  <button
    class="self-center ml-1.5 w-6 h-6 grid place-items-center rounded text-ash hover:text-ember hover:bg-white/5 border border-transparent hover:border-white/10 transition-colors text-base leading-none"
    aria-label="Nouvelle étude"
    title="Nouvelle étude (copie des paramètres courants)"
    onclick={() => store.addStudy()}
  >
    +
  </button>
</div>
