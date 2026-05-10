<script lang="ts">
  import { onMount } from "svelte";
  import Header from "./components/Header.svelte";
  import Hero from "./components/Hero.svelte";
  import ControlPanel from "./components/ControlPanel.svelte";
  import Viewer3D from "./components/Viewer3D.svelte";
  import DevelopedView from "./components/DevelopedView.svelte";
  import PrintLayout from "./components/PrintLayout.svelte";
  import { store } from "./lib/store.svelte";

  onMount(() => {
    store.compute();
  });
</script>

<Header />

<main class="no-print">
  <Hero />

  <!-- 3D + controls -->
  <section id="scene" class="max-w-[1600px] mx-auto px-6 lg:px-10 pb-12">
    <div class="grid grid-cols-1 lg:grid-cols-[380px_1fr] gap-6">
      <ControlPanel />

      <div class="glass-strong card-hairline rounded-[var(--radius-card)] overflow-hidden relative min-h-[520px] lg:min-h-[640px]">
        <Viewer3D />

        {#if store.error}
          <div class="absolute left-4 right-4 bottom-4 ember-glow rounded-xl bg-black/70 px-4 py-3 text-sm text-pearl">
            <strong class="text-ember">Erreur de calcul.</strong> {store.error}
          </div>
        {/if}
      </div>
    </div>
  </section>

  <!-- Developed views -->
  <section id="patterns" class="max-w-[1600px] mx-auto px-6 lg:px-10 pb-16">
    <div class="flex flex-col lg:flex-row lg:items-end lg:justify-between gap-3 mb-6">
      <div>
        <span class="pill">développés à plat</span>
        <h2 class="text-pearl text-2xl lg:text-3xl font-semibold tracking-tight mt-2">
          Mises à plat &mdash; échelle 1:1
        </h2>
        <p class="text-ash text-sm max-w-2xl mt-1">
          Le développé est l’unwrap métrique du cylindre.  L’abscisse est la longueur
          d’arc, l’ordonnée est la coordonnée axiale.  Tracez ces courbes au feutre
          fin, ou exportez-les en SVG vectoriel pour DXF/CNC/laser.
        </p>
      </div>
      <div class="flex items-center gap-3 text-[11px] text-ash">
        {#if store.result}
          <span class="num">
            {store.result.dev_branch.length} pts · branche
          </span>
          {#if store.result.dev_main}
            <span class="num">{store.result.dev_main.length} pts · principal</span>
          {/if}
        {/if}
      </div>
    </div>

    {#if store.result}
      <div class="grid grid-cols-1 {store.result.dev_main ? 'xl:grid-cols-2' : ''} gap-6">
        <DevelopedView
          title="Tube incliné"
          subtitle={store.result.mode === "cyl_cyl"
            ? "Gabarit de découpe Ø₂"
            : "Coupe selon plan incliné Ø₁"}
          points={store.result.dev_branch}
          diameter={store.result.mode === "cyl_cyl" && store.result.r2
            ? store.result.r2 * 2
            : store.result.r1 * 2}
          circumference={store.result.circumference_branch ?? store.result.circumference_main}
          closed={true}
          accent="#ff5b1a"
          filename="developpé-tube"
        />

        {#if store.result.dev_main && store.result.bbox_main}
          <DevelopedView
            title="Gueule de loup"
            subtitle="Cylindre principal Ø₁"
            points={store.result.dev_main}
            diameter={store.result.r1 * 2}
            circumference={store.result.circumference_main}
            closed={false}
            accent="#29c2ff"
            filename="gueule-de-loup"
          />
        {/if}
      </div>
    {/if}
  </section>

  <!-- Print preview -->
  <section class="max-w-[1600px] mx-auto px-6 lg:px-10 pb-24">
    <div class="mb-6">
      <span class="pill">impression CNC / fabrication</span>
      <h2 class="text-pearl text-2xl lg:text-3xl font-semibold tracking-tight mt-2">
        Aperçu d’impression
      </h2>
      <p class="text-ash text-sm max-w-2xl mt-1">
        Lancez l’impression depuis la barre du haut.  Le rendu masque
        automatiquement l’interface et présente les gabarits sur fond blanc, à
        l’échelle réelle.  Une règle de référence&nbsp;100&nbsp;mm permet de vérifier
        l’étalonnage.
      </p>
    </div>
    <PrintLayout />
  </section>

  <footer class="border-t border-white/5 mt-12">
    <div class="max-w-[1600px] mx-auto px-6 lg:px-10 py-8 flex flex-col sm:flex-row items-start sm:items-center gap-3 justify-between text-[11px] text-ash">
      <div>
        Moteur Rust + Tauri &middot; UI Svelte 5 + Tailwind 4 &middot;
        rendu Three.js
      </div>
      <div class="flex items-center gap-4">
        <span>tout est calculé en local sur votre machine.</span>
      </div>
    </div>
  </footer>
</main>
