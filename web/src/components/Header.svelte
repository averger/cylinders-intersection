<script lang="ts">
  import { store } from "../lib/store.svelte";

  let computedAgo = $state("—");

  $effect(() => {
    const t = store.lastComputedAt;
    if (!t) return;
    const update = () => {
      const dt = (Date.now() - t) / 1000;
      computedAgo = dt < 1 ? "à l’instant" : dt < 60 ? `il y a ${dt.toFixed(0)} s` : "—";
    };
    update();
    const id = setInterval(update, 1000);
    return () => clearInterval(id);
  });
</script>

<header class="no-print sticky top-0 z-30 backdrop-blur-xl bg-ink/60 border-b border-white/5">
  <div class="max-w-[1600px] mx-auto px-6 lg:px-10 py-4 flex items-center gap-6">
    <a href="/" class="flex items-center gap-3 group">
      <span class="relative grid place-items-center w-9 h-9 rounded-xl bg-gradient-to-b from-graphite to-black border border-white/10">
        <svg viewBox="0 0 32 32" class="w-5 h-5">
          <rect width="32" height="32" rx="6" fill="none" />
          <path d="M8 11h12a4 4 0 0 1 4 4 4 4 0 0 1-4 4H8z" fill="none" stroke="#ff5b1a" stroke-width="2" stroke-linejoin="round" />
          <line x1="14" y1="6" x2="14" y2="26" stroke="#e2e2e2" stroke-width="2" stroke-linecap="round" />
        </svg>
        <span class="absolute inset-0 rounded-xl ring-1 ring-ember/0 group-hover:ring-ember/40 transition"></span>
      </span>
      <div class="leading-tight">
        <div class="text-pearl font-semibold tracking-tight text-[15px]">Intersection</div>
        <div class="text-[10px] uppercase tracking-[0.22em] text-ash">cylindres · gabarits 1:1</div>
      </div>
    </a>

    <nav class="hidden md:flex items-center gap-1 ml-6 text-sm text-silver">
      <a href="#scene" class="px-3 py-1.5 rounded-full hover:text-pearl hover:bg-white/5 transition">Scène 3D</a>
      <a href="#patterns" class="px-3 py-1.5 rounded-full hover:text-pearl hover:bg-white/5 transition">Développés</a>
      <a href="#studio" class="px-3 py-1.5 rounded-full hover:text-pearl hover:bg-white/5 transition">Atelier</a>
      <a href="#theorie" class="px-3 py-1.5 rounded-full hover:text-pearl hover:bg-white/5 transition">Théorie</a>
    </nav>

    <div class="ml-auto flex items-center gap-3">
      <div class="hidden sm:flex items-center gap-2 text-[11px] text-ash">
        <span class="inline-block w-1.5 h-1.5 rounded-full bg-ember animate-pulse"></span>
        <span class="num">moteur Rust · {computedAgo}</span>
      </div>
      <button
        class="btn-primary"
        onclick={() => window.print()}
        aria-label="Imprimer le gabarit"
      >
        <svg viewBox="0 0 24 24" class="w-4 h-4 fill-none stroke-current" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <path d="M6 9V3h12v6" />
          <rect x="4" y="9" width="16" height="9" rx="2" />
          <path d="M8 14h8v6H8z" />
        </svg>
        Imprimer 1:1
      </button>
    </div>
  </div>
</header>
