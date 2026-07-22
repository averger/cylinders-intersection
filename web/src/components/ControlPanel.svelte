<script lang="ts">
  import Slider from "./Slider.svelte";
  import Segmented from "./Segmented.svelte";
  import ExportPanel from "./ExportPanel.svelte";
  import { store, type Mode } from "../lib/store.svelte";
  import type { Branch } from "../lib/api";

  let p = $derived(store.params);

  const modeOptions = [
    { label: "Cyl. ↺ Cyl.", value: "cyl_cyl" as Mode },
    { label: "Cyl. ↺ Plan", value: "cyl_plane" as Mode },
  ];

  const branchOptions = [
    { label: "Outer", value: "outer" as Branch },
    { label: "Inner", value: "inner" as Branch },
  ];
</script>

<aside
  class="no-print w-[320px] xl:w-[350px] shrink-0 hl-r bg-graphite/40 overflow-y-auto"
>
  <div class="p-5 flex flex-col gap-6">
    <section class="space-y-3">
      <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Mode</span>
      <Segmented
        options={modeOptions}
        value={p.mode}
        onchange={(v) => {
          store.params = { ...store.params, mode: v };
          store.compute();
        }}
      />
    </section>

    <section class="space-y-5">
      <Slider
        label="Diamètre Ø₁ (cylindre principal)"
        unit="mm"
        min={5}
        max={500}
        step={0.5}
        value={p.d1}
        decimals={2}
        onchange={(v) => {
          store.params = { ...store.params, d1: v };
          store.compute();
        }}
      />
      {#if p.mode === "cyl_cyl"}
        <Slider
          label="Diamètre Ø₂ (cylindre incliné)"
          unit="mm"
          min={5}
          max={500}
          step={0.5}
          value={p.d2}
          decimals={2}
          onchange={(v) => {
            store.params = { ...store.params, d2: v };
            store.compute();
          }}
        />
      {/if}
      <Slider
        label={p.mode === "cyl_cyl" ? "Angle entre les axes φ" : "Inclinaison du plan φ"}
        unit="°"
        min={p.mode === "cyl_cyl" ? 1 : -85}
        max={p.mode === "cyl_cyl" ? 90 : 85}
        step={0.5}
        value={p.angleDeg}
        decimals={1}
        onchange={(v) => {
          store.params = { ...store.params, angleDeg: v };
          store.compute();
        }}
      />

      {#if p.mode === "cyl_plane"}
        <Slider
          label="Décalage z₀ du plan"
          unit="mm"
          min={-200}
          max={200}
          step={0.5}
          value={p.z0}
          decimals={1}
          hint="Hauteur de coupe sur l'axe du tube"
          onchange={(v) => {
            store.params = { ...store.params, z0: v };
            store.compute();
          }}
        />
      {/if}
    </section>

    {#if p.mode === "cyl_cyl"}
      <section class="space-y-3">
        <div class="flex items-center justify-between">
          <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Branche</span>
          <span class="text-[10px] text-ash/70">lèvre conservée</span>
        </div>
        <Segmented
          options={branchOptions}
          value={p.branch}
          onchange={(v) => {
            store.params = { ...store.params, branch: v };
            store.compute();
          }}
        />
      </section>
    {/if}

    <Slider
      label="Échantillonnage"
      unit="pts"
      min={240}
      max={4800}
      step={60}
      value={p.samples}
      decimals={0}
      hint="Plus de points → contour plus lisse"
      onchange={(v) => {
        store.params = { ...store.params, samples: Math.round(v) };
        store.compute();
      }}
    />

    {#if store.view === "2d"}
      <div class="pt-5 border-t border-mist">
        <ExportPanel />
      </div>
    {:else}
      <footer class="pt-4 hl-t space-y-2 text-[11px] text-ash">
        <div class="flex justify-between">
          <span>Périmètre cyl. principal</span>
          <span class="num text-silver"
            >{store.result ? store.result.circumference_main.toFixed(2) : "—"} mm</span
          >
        </div>
        {#if p.mode === "cyl_cyl"}
          <div class="flex justify-between">
            <span>Périmètre cyl. branche</span>
            <span class="num text-silver"
              >{store.result?.circumference_branch
                ? store.result.circumference_branch.toFixed(2)
                : "—"} mm</span
            >
          </div>
        {/if}
        <div class="flex justify-between">
          <span>Points d'intersection</span>
          <span class="num text-silver">{store.result?.curve3d.length ?? 0}</span>
        </div>
      </footer>
    {/if}
  </div>
</aside>
