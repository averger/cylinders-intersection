<script lang="ts">
  import { store } from "../lib/store.svelte";
  import { buildPatternSVG } from "../lib/svg";

  // We render patterns at print scale — 1 mm = 1 SVG user unit, with the
  // <svg> dimensioned in millimetres so that the printer reproduces 1:1.
  let mainSVG = $derived.by(() => {
    const r = store.result;
    if (!r || !r.dev_main) return null;
    return buildPatternSVG(r.dev_main, {
      title: "Gueule de loup — cylindre principal",
      subtitle: `Ø ${(r.r1 * 2).toFixed(2)} mm — angle ${((r.phi * 180) / Math.PI).toFixed(2)}°`,
      diameter: r.r1 * 2,
      circumference: r.circumference_main,
      closed: false,
    });
  });

  let branchSVG = $derived.by(() => {
    const r = store.result;
    if (!r) return null;
    const isCyl = r.mode === "cyl_cyl";
    const d = isCyl && r.r2 != null ? r.r2 * 2 : r.r1 * 2;
    return buildPatternSVG(r.dev_branch, {
      title: isCyl ? "Tube incliné — gabarit de coupe" : "Tube — coupe selon plan incliné",
      subtitle: `Ø ${d.toFixed(2)} mm — angle ${((r.phi * 180) / Math.PI).toFixed(2)}°`,
      diameter: d,
      circumference: r.circumference_branch ?? r.circumference_main,
      closed: true,
    });
  });
</script>

<section id="print" class="print-area space-y-8">
  <div class="no-print pill">prévisualisation impression · 1:1 garanti</div>

  {#if branchSVG}
    <article class="bg-white text-black rounded-[14px] p-4 lg:p-8 border border-mist/30 shadow-2xl">
      <div class="flex flex-wrap items-baseline justify-between gap-2 no-print">
        <h3 class="text-lg font-semibold tracking-tight text-black">
          Gabarit du tube {store.result?.mode === "cyl_cyl" ? "incliné" : "coupé"}
        </h3>
        <span class="num text-xs text-mist">
          Ø {(store.result?.mode === "cyl_cyl" && store.result.r2
            ? store.result.r2 * 2
            : (store.result?.r1 ?? 0) * 2).toFixed(2)} mm
        </span>
      </div>
      <div class="overflow-auto scroll-fade max-w-full">
        <!-- bind raw SVG: no escaping needed because we generate it ourselves -->
        {@html branchSVG}
      </div>
    </article>
  {/if}

  {#if mainSVG}
    <article class="bg-white text-black rounded-[14px] p-4 lg:p-8 border border-mist/30 shadow-2xl">
      <div class="flex flex-wrap items-baseline justify-between gap-2 no-print">
        <h3 class="text-lg font-semibold tracking-tight text-black">
          Gueule de loup — cylindre principal
        </h3>
        <span class="num text-xs text-mist">
          Ø {((store.result?.r1 ?? 0) * 2).toFixed(2)} mm
        </span>
      </div>
      <div class="overflow-auto scroll-fade max-w-full">
        {@html mainSVG}
      </div>
    </article>
  {/if}

  <p class="no-print text-[11px] text-ash">
    Astuce impression — sélectionnez « Échelle réelle » / « 100 % » / « Actual size »
    et désactivez « ajuster à la page ». La règle de 100 mm imprimée doit mesurer
    exactement 100 mm une fois sur papier ; sinon, ajustez le réglage de votre pilote.
  </p>
</section>
