<script lang="ts">
  import { store } from "../lib/store.svelte";

  let mode = $derived(store.params.mode);
</script>

<section class="relative pt-12 pb-10 lg:pt-20 lg:pb-14 overflow-hidden">
  <!-- ambient glow -->
  <div
    class="absolute -top-32 right-1/4 w-[60vw] h-[60vw] rounded-full blur-3xl opacity-30 pointer-events-none"
    style="background: radial-gradient(circle at center, rgba(255,91,26,0.55), transparent 60%);"
  ></div>

  <div class="relative max-w-[1600px] mx-auto px-6 lg:px-10">
    <div class="flex flex-col lg:flex-row lg:items-end gap-8 lg:gap-12">
      <div class="flex-1 space-y-5 max-w-3xl">
        <span class="pill">
          <span class="inline-block w-1.5 h-1.5 rounded-full bg-ember"></span>
          calculs matriciels en Rust · UI Svelte
        </span>
        <h1 class="font-display font-semibold tracking-[-0.02em] text-pearl text-4xl lg:text-6xl leading-[1.02]">
          Découpez deux cylindres
          <span class="block bg-gradient-to-r from-ember via-ember-soft to-cyan-glow bg-clip-text text-transparent">
            avec la précision d’une CNC.
          </span>
        </h1>
        <p class="text-silver text-base lg:text-lg max-w-2xl leading-relaxed">
          Définissez les diamètres et l’angle. Le moteur calcule l’intersection, déroule
          la <em class="text-pearl not-italic font-medium">gueule de loup</em> sur le tube
          principal et le <em class="text-pearl not-italic font-medium">développé</em>
          sur le tube incliné. Tout est exporté à l’échelle&nbsp;1&thinsp;:&thinsp;1 pour
          impression, traçage ou usinage.
        </p>
        <div class="flex flex-wrap items-center gap-3 text-[12px] text-ash">
          <span class="pill" style="border-color: rgba(255,91,26,0.4); color:#ffb597">
            mode actif · {mode === "cyl_cyl" ? "cylindre × cylindre" : "cylindre × plan incliné"}
          </span>
          <span class="num">tolérance numérique &lt; 1e-9</span>
          <span class="hidden md:inline">·</span>
          <span class="hidden md:inline">SVG vectoriel ·  CSV brut · table d’usinage</span>
        </div>
      </div>

      <div class="grid grid-cols-3 gap-3 lg:gap-4 shrink-0 self-end">
        <div class="glass card-hairline px-4 py-4 min-w-[110px]">
          <div class="text-[10px] uppercase tracking-[0.22em] text-ash">Ø₁</div>
          <div class="num text-pearl text-2xl mt-1">{store.params.d1.toFixed(1)}</div>
          <div class="text-[10px] text-ash">millimètres</div>
        </div>
        {#if mode === "cyl_cyl"}
          <div class="glass card-hairline px-4 py-4">
            <div class="text-[10px] uppercase tracking-[0.22em] text-ash">Ø₂</div>
            <div class="num text-pearl text-2xl mt-1">{store.params.d2.toFixed(1)}</div>
            <div class="text-[10px] text-ash">millimètres</div>
          </div>
        {:else}
          <div class="glass card-hairline px-4 py-4">
            <div class="text-[10px] uppercase tracking-[0.22em] text-ash">z₀</div>
            <div class="num text-pearl text-2xl mt-1">{store.params.z0.toFixed(1)}</div>
            <div class="text-[10px] text-ash">millimètres</div>
          </div>
        {/if}
        <div class="glass card-hairline px-4 py-4">
          <div class="text-[10px] uppercase tracking-[0.22em] text-ash">φ</div>
          <div class="num text-pearl text-2xl mt-1">{store.params.angleDeg.toFixed(1)}°</div>
          <div class="text-[10px] text-ash">angle d’axe</div>
        </div>
      </div>
    </div>
  </div>
</section>
