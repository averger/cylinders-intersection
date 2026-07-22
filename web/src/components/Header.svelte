<script lang="ts">
  // Compact tool header: brand, view toggle, actions.
  import Segmented from "./Segmented.svelte";
  import { store, type View } from "../lib/store.svelte";

  const viewOptions = [
    { label: "3D", value: "3d" as View },
    { label: "2D", value: "2d" as View },
  ];
</script>

<header class="no-print shrink-0 hl-b bg-ink/70 backdrop-blur-xl z-30">
  <div class="px-4 lg:px-6 h-14 flex items-center gap-5">
    <a href="/" class="flex items-center gap-3 group shrink-0">
      <span
        class="relative grid place-items-center w-9 h-9 rounded-xl bg-carbon border border-mist"
      >
        <svg viewBox="0 0 32 32" class="w-5 h-5">
          <path
            d="M8 11h12a4 4 0 0 1 4 4 4 4 0 0 1-4 4H8z"
            fill="none"
            style="stroke: var(--ember)"
            stroke-width="2"
            stroke-linejoin="round"
          />
          <line x1="14" y1="6" x2="14" y2="26" style="stroke: var(--text-2)" stroke-width="2" stroke-linecap="round" />
        </svg>
        <span class="absolute inset-0 rounded-xl ring-1 ring-ember/0 group-hover:ring-ember/40 transition"></span>
      </span>
      <div class="leading-tight">
        <div class="text-pearl font-semibold tracking-tight text-[16px]">Cylix</div>
        <div class="text-[9px] uppercase tracking-[0.22em] text-ash">
          gabarits de découpe · échelle 1:1
        </div>
      </div>
    </a>

    <div class="mx-auto">
      <Segmented options={viewOptions} value={store.view} onchange={(v) => store.setView(v)} />
    </div>

    <div class="flex items-center gap-2 shrink-0">
      <button
        class="btn-ghost !py-2 !px-3"
        onclick={() => store.toggleTheme()}
        aria-label={store.theme === "light" ? "Passer en thème sombre" : "Passer en thème clair"}
        title={store.theme === "light" ? "Thème sombre" : "Thème clair"}
      >
        {#if store.theme === "light"}
          <svg viewBox="0 0 24 24" class="w-4 h-4 fill-none stroke-current" stroke-width="2" stroke-linecap="round">
            <path d="M21 12.8A9 9 0 1 1 11.2 3 7 7 0 0 0 21 12.8z" stroke-linejoin="round" />
          </svg>
        {:else}
          <svg viewBox="0 0 24 24" class="w-4 h-4 fill-none stroke-current" stroke-width="2" stroke-linecap="round">
            <circle cx="12" cy="12" r="4" />
            <path d="M12 2v2m0 16v2M4.9 4.9l1.4 1.4m11.4 11.4 1.4 1.4M2 12h2m16 0h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4" />
          </svg>
        {/if}
      </button>
      <button class="btn-ghost !py-2 text-sm" onclick={() => window.print()} aria-label="Imprimer 1:1">
        <svg
          viewBox="0 0 24 24"
          class="w-4 h-4 fill-none stroke-current"
          stroke-width="2"
          stroke-linecap="round"
          stroke-linejoin="round"
        >
          <path d="M6 9V3h12v6" />
          <rect x="4" y="9" width="16" height="9" rx="2" />
          <path d="M8 14h8v6H8z" />
        </svg>
        Imprimer
      </button>
    </div>
  </div>
</header>
