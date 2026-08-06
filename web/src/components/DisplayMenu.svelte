<script lang="ts">
  // Header dropdown for the 3D scene layers — keeps the viewport free of
  // any chrome.  Toggling a layer only flips visibility flags, no rebuild.
  import { store, type Show3D } from "../lib/store.svelte";

  let open = $state(false);

  const rows = $derived(
    [
      { key: "main", label: "tube principal" },
      {
        key: "cutters",
        label:
          store.params.mode === "multi"
            ? "piquages"
            : store.params.mode === "cyl_plane"
              ? "plan de coupe"
              : "tube incliné",
      },
      { key: "curves", label: "courbes de coupe" },
      { key: "angles", label: "angles φ · ψ" },
      { key: "labels", label: "étiquettes Ø" },
      { key: "axes", label: "axes · cotes Ø" },
      { key: "grid", label: "grille" },
    ] as { key: keyof Show3D; label: string }[],
  );

  function clickOutside(node: HTMLElement) {
    const onDown = (e: PointerEvent) => {
      if (!node.contains(e.target as Node)) open = false;
    };
    document.addEventListener("pointerdown", onDown, true);
    return { destroy: () => document.removeEventListener("pointerdown", onDown, true) };
  }
</script>

<div class="relative" use:clickOutside>
  <button
    class="btn-ghost !py-2 !px-3 text-[13px]"
    aria-haspopup="menu"
    aria-expanded={open}
    onclick={() => (open = !open)}
  >
    Affichage
    <svg
      viewBox="0 0 24 24"
      class="w-3.5 h-3.5 fill-none stroke-current transition-transform {open ? 'rotate-180' : ''}"
      stroke-width="2.5"
      stroke-linecap="round"
      stroke-linejoin="round"
    >
      <path d="m6 9 6 6 6-6" />
    </svg>
  </button>

  {#if open}
    <div
      class="absolute left-1/2 -translate-x-1/2 top-[calc(100%+10px)] w-[230px] glass-strong card-hairline rounded-2xl p-2 z-50"
      style="background: var(--bg-2); backdrop-filter: none; -webkit-backdrop-filter: none"
      role="menu"
    >
      {#each rows as row (row.key)}
        <button
          class="w-full flex items-center gap-2.5 px-2.5 py-1.5 rounded-lg hover:bg-carbon/60 transition-colors text-left"
          role="menuitemcheckbox"
          aria-checked={store.show3d[row.key]}
          onclick={() => store.toggle3d(row.key)}
        >
          <span
            class="w-3 h-3 rounded-[4px] border transition-colors shrink-0"
            style={store.show3d[row.key]
              ? "background: var(--ember); border-color: var(--ember)"
              : "border-color: var(--line-2)"}
          ></span>
          <span
            class="text-[10px] uppercase tracking-[0.14em] {store.show3d[row.key]
              ? 'text-silver'
              : 'text-ash/60'}">{row.label}</span
          >
        </button>
      {/each}
      <div class="mt-1 pt-1 border-t border-mist">
        <button
          class="w-full flex items-center gap-2.5 px-2.5 py-1.5 rounded-lg hover:bg-carbon/60 transition-colors text-left"
          role="menuitemcheckbox"
          aria-checked={store.autoRotate}
          onclick={() => (store.autoRotate = !store.autoRotate)}
        >
          <span
            class="w-3 h-3 rounded-[4px] border transition-colors shrink-0"
            style={store.autoRotate
              ? "background: var(--ember); border-color: var(--ember)"
              : "border-color: var(--line-2)"}
          ></span>
          <span
            class="text-[10px] uppercase tracking-[0.14em] {store.autoRotate
              ? 'text-silver'
              : 'text-ash/60'}">rotation auto</span
          >
        </button>
      </div>
    </div>
  {/if}
</div>
