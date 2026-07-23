<script lang="ts">
  // Slim bottom status bar — live geometry summary + engine heartbeat.
  import { store } from "../lib/store.svelte";

  let computedAgo = $state("—");
  $effect(() => {
    const t = store.lastComputedAt;
    if (!t) return;
    const update = () => {
      const dt = (Date.now() - t) / 1000;
      computedAgo = dt < 1 ? "à l'instant" : dt < 60 ? `il y a ${dt.toFixed(0)} s` : `il y a ${Math.floor(dt / 60)} min`;
    };
    update();
    const id = setInterval(update, 1000);
    return () => clearInterval(id);
  });

  let p = $derived(store.params);
  let r = $derived(store.result);
</script>

<footer
  class="no-print shrink-0 h-7 px-4 lg:px-6 flex items-center justify-between gap-4 hl-t bg-carbon/50 text-[10px] text-ash num overflow-hidden"
>
  <div class="flex items-center gap-3 truncate">
    <span class="text-silver">
      {p.mode === "cyl_cyl"
        ? "cylindre × cylindre"
        : p.mode === "cyl_plane"
          ? "cylindre × plan"
          : "châssis · nœud multi-piquages"}
    </span>
    <span>Ø₁ {p.d1.toFixed(1)}</span>
    {#if p.mode === "cyl_cyl"}<span>Ø₂ {p.d2.toFixed(1)}</span>{/if}
    {#if p.mode === "multi"}
      <span>{p.branches.length} piquage{p.branches.length > 1 ? "s" : ""}</span>
      {#if store.multiResult}
        <span class="hidden sm:inline"
          >périm. {store.multiResult.circumference_main.toFixed(1)} mm</span
        >
      {/if}
    {:else}
      <span>{p.mode === "cyl_plane" ? "φx" : "φ"} {p.angleDeg.toFixed(1)}°</span>
      {#if p.mode === "cyl_plane"}<span>φy {p.angleYDeg.toFixed(1)}°</span><span>z₀ {p.z0.toFixed(1)}</span>{/if}
      {#if r}
        <span class="hidden sm:inline">périm. {r.circumference_main.toFixed(1)} mm</span>
        <span class="hidden md:inline">{r.curve3d.length} pts</span>
      {/if}
    {/if}
  </div>
  <div class="flex items-center gap-2 shrink-0">
    {#if store.error}
      <span class="text-ember truncate max-w-[40vw]">{store.error}</span>
    {:else}
      <span class="inline-block w-1.5 h-1.5 rounded-full bg-ember {store.loading ? 'animate-pulse' : ''}"></span>
      <span>moteur Rust · {store.loading ? "calcul…" : computedAgo}</span>
    {/if}
  </div>
</footer>
