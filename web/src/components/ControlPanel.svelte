<script lang="ts">
  import Slider from "./Slider.svelte";
  import Segmented from "./Segmented.svelte";
  import Switch from "./Switch.svelte";
  import ExportPanel from "./ExportPanel.svelte";
  import { store, type Mode } from "../lib/store.svelte";
  import type { Branch } from "../lib/api";

  let p = $derived(store.params);

  const modeOptions = [
    { label: "Cyl / Cyl", value: "cyl_cyl" as Mode },
    { label: "Cyl / Plan", value: "cyl_plane" as Mode },
    { label: "Multi tube", value: "multi" as Mode },
  ];

  function patchBranch(index: number, patch: Partial<(typeof store.params.branches)[number]>) {
    const branches = store.params.branches.map((b, i) => (i === index ? { ...b, ...patch } : b));
    store.params = { ...store.params, branches };
    store.compute();
  }

  function addBranch() {
    if (store.params.branches.length >= 8) return;
    const last = store.params.branches[store.params.branches.length - 1];
    // Axes are concurrent: an exact clone would die entirely on its twin —
    // tilt the newcomer so it lands somewhere of its own.
    const angleDeg = last.angleDeg + 30 <= 175 ? last.angleDeg + 30 : last.angleDeg - 30;
    const branches = [...store.params.branches, { ...last, angleDeg, z: last.z }];
    store.params = { ...store.params, branches };
    store.compute();
  }

  function removeBranch(index: number) {
    if (store.params.branches.length <= 1) return;
    const branches = store.params.branches.filter((_, i) => i !== index);
    store.params = { ...store.params, branches };
    store.compute();
  }

  // The second plane tilt is rare — revealed on demand.
  let showPhiY = $state(false);
  $effect(() => {
    if (store.params.angleYDeg !== 0) showPhiY = true;
  });

  const branchOptions = [
    { label: "Extérieure", value: "outer" as Branch },
    { label: "Intérieure", value: "inner" as Branch },
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
          // Un piquage plus gros que le tube principal traverse de part en
          // part : hors du domaine gueule de loup — on borne Ø₂ à Ø₁.
          store.params = {
            ...store.params,
            d1: v,
            d2: Math.min(store.params.d2, v),
            branches: store.params.branches.map((b) => ({ ...b, d: Math.min(b.d, v) })),
          };
          store.compute();
        }}
      />
      {#if p.mode === "cyl_cyl"}
        <Slider
          label="Diamètre Ø₂ (cylindre incliné)"
          unit="mm"
          min={5}
          max={p.d1}
          step={0.5}
          value={p.d2}
          hint="au plus égal à Ø₁"
          decimals={2}
          onchange={(v) => {
            store.params = { ...store.params, d2: v };
            store.compute();
          }}
        />
      {/if}
      {#if p.mode !== "multi"}
        <Slider
          label={p.mode === "cyl_cyl" ? "Angle entre les axes φ" : "Inclinaison selon X · φx"}
          unit="°"
          min={p.mode === "cyl_cyl" ? 1 : -85}
          max={p.mode === "cyl_cyl" ? 90 : 85}
          step={0.5}
          value={p.angleDeg}
          decimals={1}
          hint={p.mode === "cyl_plane" ? "bascule du plan autour de l'axe X" : undefined}
          onchange={(v) => {
            store.params = { ...store.params, angleDeg: v };
            store.compute();
          }}
        />
      {/if}

      {#if p.mode === "cyl_plane"}
        <div class="flex items-center justify-between -mt-1">
          <div class="flex flex-col">
            <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Plan orienté</span>
            <span class="text-[10px] text-ash/70">second angle φy</span>
          </div>
          <Switch
            checked={showPhiY}
            label="Plan orienté — second angle φy"
            onchange={(v) => {
              showPhiY = v;
              if (!v && store.params.angleYDeg !== 0) {
                store.params = { ...store.params, angleYDeg: 0 };
                store.compute();
              }
            }}
          />
        </div>
        {#if showPhiY}
          <Slider
            label="Inclinaison selon Y · φy"
            unit="°"
            min={-85}
            max={85}
            step={0.5}
            value={p.angleYDeg}
            decimals={1}
            hint="déphase le gabarit autour du tube — caler sur la génératrice 0°"
            onchange={(v) => {
              store.params = { ...store.params, angleYDeg: v };
              store.compute();
            }}
          />
        {/if}
      {/if}

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

    {#if p.mode === "multi"}
      <section class="space-y-3">
        <div class="flex items-center justify-between">
          <div class="flex flex-col">
            <span class="text-[10px] uppercase tracking-[0.18em] text-ash">
              Piquages ({p.branches.length})
            </span>
            <span class="text-[10px] text-ash/70">
              l'ordre donne la priorité : un piquage meurt sur ceux au-dessus
            </span>
          </div>
          <button
            class="text-[11px] text-ember hover:text-ember-soft transition-colors font-medium disabled:opacity-40"
            disabled={p.branches.length >= 8}
            onclick={addBranch}
          >
            + Ajouter
          </button>
        </div>

        {#each p.branches as b, i (i)}
          <div class="rounded-xl border border-mist bg-carbon/50 p-3 space-y-4">
            <div class="flex items-center justify-between">
              <span class="text-[11px] font-medium text-silver">Piquage {i + 1}</span>
              {#if p.branches.length > 1}
                <button
                  class="text-ash hover:text-ember transition-colors text-sm leading-none"
                  aria-label={`Supprimer le piquage ${i + 1}`}
                  onclick={() => removeBranch(i)}
                >
                  ×
                </button>
              {/if}
            </div>
            <Slider
              label="Diamètre Ø"
              unit="mm"
              min={5}
              max={p.d1}
              step={0.5}
              value={b.d}
              decimals={1}
              onchange={(v) => patchBranch(i, { d: v })}
            />
            <Slider
              label="Excentrement z"
              unit="mm"
              min={-200}
              max={200}
              step={1}
              value={b.z}
              decimals={0}
              hint="0 = axe au centre du tube · décale le point de croisement"
              onchange={(v) => patchBranch(i, { z: v })}
            />
            <Slider
              label="Inclinaison φ"
              unit="°"
              min={5}
              max={175}
              step={0.5}
              value={b.angleDeg}
              decimals={1}
              hint="depuis l'axe principal — > 90° : penche vers le bas"
              onchange={(v) => patchBranch(i, { angleDeg: v })}
            />
            <Slider
              label="Azimut ψ"
              unit="°"
              min={-180}
              max={180}
              step={1}
              value={b.azimutDeg}
              decimals={0}
              hint="rotation autour du tube principal"
              onchange={(v) => patchBranch(i, { azimutDeg: v })}
            />
          </div>
        {/each}

        {#if store.multiResult && store.multiResult.pairs.length > 0}
          <div class="rounded-xl border border-mist bg-carbon/30 p-3 space-y-2">
            <div class="flex flex-col">
              <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Nœud</span>
              <span class="text-[10px] text-ash/70">
                excentrement, jeu ou recouvrement — les cotes du bureau d'études
              </span>
            </div>
            {#each store.multiResult.pairs as pr (`${pr.i}-${pr.j}`)}
              <div class="space-y-1 text-[11px]">
                <div class="flex justify-between">
                  <span class="text-ash"
                    >P{pr.i + 1} · P{pr.j + 1}
                    <span class="text-ash/60">{pr.same_side ? "même côté" : "opposés"}</span></span
                  >
                  <span class="num text-silver">
                    e = {pr.eccentricity === null ? "—" : pr.eccentricity.toFixed(1)} mm
                  </span>
                </div>
                {#if pr.overlap_pct !== null}
                  <div class="flex justify-between">
                    <span class="text-ash">recouvrement λov</span>
                    <span
                      class="num {pr.overlap_pct < 25 ? 'text-ember' : 'text-silver'}"
                      >{pr.overlap_pct.toFixed(0)} %</span
                    >
                  </div>
                {:else if pr.gap !== null}
                  <div class="flex justify-between">
                    <span class="text-ash">jeu g</span>
                    <span class="num text-silver">{pr.gap.toFixed(1)} mm</span>
                  </div>
                {/if}
              </div>
            {/each}
          </div>
        {/if}

        {#if store.multiResult && store.multiResult.warnings.length > 0}
          <div class="space-y-1.5">
            {#each store.multiResult.warnings as w, wi (wi)}
              <p class="text-[10px] text-ember/90 leading-relaxed">⚠ {w}</p>
            {/each}
          </div>
        {/if}
      </section>
    {/if}

    {#if p.mode === "cyl_cyl"}
      <section class="space-y-3">
        <div class="flex items-center justify-between">
          <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Lèvre de coupe</span>
        </div>
        <Segmented
          options={branchOptions}
          value={p.branch}
          onchange={(v) => {
            store.params = { ...store.params, branch: v };
            store.compute();
          }}
        />
        <p class="text-[10px] text-ash/70 leading-relaxed">
          Côté d'où arrive le tube incliné : <b>extérieure</b> = il s'appuie sur le
          gros tube et s'arrête au premier contact (piquage en selle, le cas
          courant) ; <b>intérieure</b> = il arrive du côté opposé — gabarit miroir.
        </p>
      </section>
    {/if}

    {#if store.view === "2d"}
      <div class="pt-5 border-t border-mist">
        <ExportPanel />
      </div>
    {:else}
      <footer class="pt-4 hl-t space-y-2 text-[11px] text-ash">
        {#if p.mode === "multi"}
          <div class="flex justify-between">
            <span>Périmètre cyl. principal</span>
            <span class="num text-silver"
              >{store.multiResult ? store.multiResult.circumference_main.toFixed(2) : "—"} mm</span
            >
          </div>
          <div class="flex justify-between">
            <span>Piquages</span>
            <span class="num text-silver">{store.multiResult?.branches.length ?? 0}</span>
          </div>
          <div class="flex justify-between">
            <span>Coutures mutuelles</span>
            <span class="num text-silver"
              >{store.multiResult?.branches.filter((b) => b.cut_by_neighbor).length ?? 0}</span
            >
          </div>
        {:else}
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
        {/if}
      </footer>
    {/if}
  </div>
</aside>
