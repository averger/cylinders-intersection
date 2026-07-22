<script lang="ts">
  import { onMount } from "svelte";
  import Header from "./components/Header.svelte";
  import StudyTabs from "./components/StudyTabs.svelte";
  import ControlPanel from "./components/ControlPanel.svelte";
  import Viewer3D from "./components/Viewer3D.svelte";
  import Studio2D from "./components/Studio2D.svelte";
  import StatusBar from "./components/StatusBar.svelte";
  import PrintLayout from "./components/PrintLayout.svelte";
  import { store } from "./lib/store.svelte";

  onMount(() => {
    store.compute();
  });

  $effect(() => {
    document.documentElement.dataset.theme = store.theme;
  });
</script>

<div class="h-screen flex flex-col overflow-hidden">
  <Header />
  <StudyTabs />

  <div class="flex-1 flex min-h-0">
    <ControlPanel />

    <main class="flex-1 relative min-w-0 no-print">
      {#if store.view === "3d"}
        <Viewer3D />
      {:else}
        <Studio2D />
      {/if}

      {#if store.error}
        <div
          class="absolute left-4 right-4 bottom-4 ember-glow rounded-xl glass-strong px-4 py-3 text-sm text-pearl z-20"
        >
          <strong class="text-ember">Erreur de calcul.</strong>
          {store.error}
        </div>
      {/if}
    </main>
  </div>

  <StatusBar />
</div>

<!-- Print host: hidden on screen, becomes the whole page in @media print. -->
<div class="print-host" aria-hidden="true">
  <PrintLayout />
</div>

<style>
  .print-host {
    position: fixed;
    left: -200vw;
    top: 0;
    width: 100vw;
    pointer-events: none;
  }
  @media print {
    .print-host {
      position: static;
      left: 0;
      pointer-events: auto;
    }
  }
</style>
