<script lang="ts">
  // Compact tool header: brand, view toggle, options, actions.
  import Segmented from "./Segmented.svelte";
  import OptionsMenu from "./OptionsMenu.svelte";
  import { store, type View } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";

  let exportMenu = $state(false);

  function pick(kind: "pdf" | "dxf") {
    exportMenu = false;
    if (kind === "pdf") editor.exportPdf();
    else editor.exportDxf();
  }

  function clickOutside(node: HTMLElement) {
    const onDown = (e: PointerEvent) => {
      if (!node.contains(e.target as Node)) exportMenu = false;
    };
    document.addEventListener("pointerdown", onDown, true);
    return { destroy: () => document.removeEventListener("pointerdown", onDown, true) };
  }

  const viewOptions = [
    { label: "3D", value: "3d" as View },
    { label: "2D", value: "2d" as View },
  ];
</script>

<header class="no-print shrink-0 hl-b bg-ink/70 backdrop-blur-xl z-30">
  <div class="px-4 lg:px-6 h-14 flex items-center gap-5">
    <a href="/" class="flex items-center gap-3 group shrink-0">
      <span class="relative grid place-items-center w-9 h-9">
        <!-- Cylix mark: the development sinusoid on its tile -->
        <svg viewBox="0 0 96 96" class="w-9 h-9" fill="none" aria-hidden="true">
          <defs>
            <linearGradient id="cylix-mark-g" x1="0" y1="0" x2="1" y2="1">
              <stop offset="0" stop-color="#FF7A33" />
              <stop offset="1" stop-color="#E8500F" />
            </linearGradient>
          </defs>
          <rect width="96" height="96" rx="24" fill="#16161A" />
          <path
            d="M16.0 48.0 L17.8 45.0 L19.6 42.2 L21.3 39.5 L23.1 37.1 L24.9 35.0 L26.7 33.3 L28.4 32.0 L30.2 31.3 L32.0 31.0 L33.8 31.3 L35.6 32.0 L37.3 33.3 L39.1 35.0 L40.9 37.1 L42.7 39.5 L44.4 42.2 L46.2 45.0 L48.0 48.0 L49.8 51.0 L51.6 53.8 L53.3 56.5 L55.1 58.9 L56.9 61.0 L58.7 62.7 L60.4 64.0 L62.2 64.7 L64.0 65.0 L65.8 64.7 L67.6 64.0 L69.3 62.7 L71.1 61.0 L72.9 58.9 L74.7 56.5 L76.4 53.8 L78.2 51.0 L80.0 48.0"
            stroke="url(#cylix-mark-g)"
            stroke-width="9"
            stroke-linecap="round"
            stroke-linejoin="round"
          />
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

    <div class="mx-auto flex items-center gap-2">
      <Segmented options={viewOptions} value={store.view} onchange={(v) => store.setView(v)} />
      <OptionsMenu />
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
      <div class="relative" use:clickOutside>
        <button
          class="btn-primary !py-2 text-sm"
          disabled={editor.busy !== null || !store.result}
          onclick={() => (exportMenu = !exportMenu)}
          aria-haspopup="menu"
          aria-expanded={exportMenu}
          aria-label="Exporter"
          title="Exporter — le cliché exact des vues 2D, en PDF 1:1 ou DXF"
        >
          <svg
            viewBox="0 0 24 24"
            class="w-4 h-4 fill-none stroke-current"
            stroke-width="2"
            stroke-linecap="round"
            stroke-linejoin="round"
          >
            <path d="M12 3v12" />
            <path d="m7 10 5 5 5-5" />
            <path d="M4 19h16" />
          </svg>
          {editor.busy !== null ? "export…" : "Exporter"}
          <svg
            viewBox="0 0 24 24"
            class="w-3.5 h-3.5 fill-none stroke-current transition-transform {exportMenu
              ? 'rotate-180'
              : ''}"
            stroke-width="2.5"
            stroke-linecap="round"
            stroke-linejoin="round"
          >
            <path d="m6 9 6 6 6-6" />
          </svg>
        </button>

        {#if exportMenu}
          <div
            class="absolute right-0 top-[calc(100%+8px)] w-64 glass-strong card-hairline rounded-2xl p-1.5 z-50"
            role="menu"
          >
            <button
              class="w-full text-left px-3.5 py-2.5 rounded-xl hover:bg-carbon transition-colors"
              role="menuitem"
              onclick={() => pick("pdf")}
            >
              <div class="text-[13.5px] font-medium text-pearl">PDF · échelle 1:1</div>
              <div class="text-[11px] text-ash leading-snug mt-0.5">
                Tuilé avec repères de collage — imprimer, enrouler, couper.
              </div>
            </button>
            <button
              class="w-full text-left px-3.5 py-2.5 rounded-xl hover:bg-carbon transition-colors"
              role="menuitem"
              onclick={() => pick("dxf")}
            >
              <div class="text-[13.5px] font-medium text-pearl">DXF · CAO / CNC</div>
              <div class="text-[11px] text-ash leading-snug mt-0.5">
                R12 en mm, calques CUT / FRAME / AXIS — laser, plasma, CAO.
              </div>
            </button>
          </div>
        {/if}
      </div>
    </div>
  </div>
</header>
