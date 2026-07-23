<script lang="ts">
  // Header dropdown for layout/export options — keeps the left panel clean.
  import Segmented from "./Segmented.svelte";
  import Slider from "./Slider.svelte";
  import { editor } from "../lib/editor.svelte";
  import type { PageFormat, PageOrientation, ScaleMode } from "../lib/export";

  let open = $state(false);

  const formatOptions = [
    { label: "A4", value: "a4" as PageFormat },
    { label: "A3", value: "a3" as PageFormat },
    { label: "A2", value: "a2" as PageFormat },
  ];
  const orientationOptions = [
    { label: "Paysage", value: "landscape" as PageOrientation },
    { label: "Portrait", value: "portrait" as PageOrientation },
  ];
  const scaleOptions = [
    { label: "1:1 · tuiles", value: "one_to_one" as ScaleMode },
    { label: "Ajusté", value: "fit" as ScaleMode },
  ];

  const layerDefs: [keyof typeof editor.layers & string, string][] = [
    ["grid", "Grille 10 mm"],
    ["frame", "Emprise du tube"],
    ["axis", "Génératrices 90°"],
    ["labels", "Étiquettes"],
    ["scale_bar", "Règle 100 mm"],
  ];

  function clickOutside(node: HTMLElement) {
    const onDown = (e: PointerEvent) => {
      if (!node.contains(e.target as Node)) open = false;
    };
    document.addEventListener("pointerdown", onDown, true);
    return { destroy: () => document.removeEventListener("pointerdown", onDown, true) };
  }
</script>

<div class="relative" use:clickOutside>
  <button
    class="btn-ghost !py-2 !px-3 text-[13px]"
    aria-haspopup="menu"
    aria-expanded={open}
    onclick={() => (open = !open)}
  >
    Options
    <svg
      viewBox="0 0 24 24"
      class="w-3.5 h-3.5 fill-none stroke-current transition-transform {open ? 'rotate-180' : ''}"
      stroke-width="2.5"
      stroke-linecap="round"
      stroke-linejoin="round"
    >
      <path d="m6 9 6 6 6-6" />
    </svg>
  </button>

  {#if open}
    <div
      class="absolute left-1/2 -translate-x-1/2 top-[calc(100%+10px)] w-[360px] glass-strong card-hairline rounded-2xl p-5 z-50 flex flex-col gap-6 max-h-[calc(100vh-96px)] overflow-y-auto"
      style="background: var(--bg-2); backdrop-filter: none; -webkit-backdrop-filter: none"
    >
      <section class="space-y-3">
        <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Mise en page</span>
        <div class="flex flex-wrap gap-2">
          <Segmented
            options={formatOptions}
            value={editor.page.format}
            onchange={(v) => (editor.page = { ...editor.page, format: v })}
          />
          <Segmented
            options={orientationOptions}
            value={editor.page.orientation}
            onchange={(v) => (editor.page = { ...editor.page, orientation: v })}
          />
        </div>
        <Segmented
          options={scaleOptions}
          value={editor.scale}
          onchange={(v) => (editor.scale = v)}
        />
        <p class="text-[10px] text-ash/70 leading-relaxed">
          1:1 = tracé direct sur tôle/tube après collage des tuiles (repères inclus).
          « Ajusté » n'est qu'un aperçu, jamais pour le traçage.
        </p>
      </section>

      <section class="space-y-3">
        <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Calques</span>
        <div class="grid grid-cols-2 gap-x-3 gap-y-2">
          {#each layerDefs as [key, label] (key)}
            <label class="flex items-center gap-2 cursor-pointer text-[13px] text-silver">
              <input
                type="checkbox"
                class="accent-[#ff5b1a] w-3.5 h-3.5"
                checked={editor.layers[key]}
                onchange={(e) =>
                  (editor.layers = {
                    ...editor.layers,
                    [key]: (e.currentTarget as HTMLInputElement).checked,
                  })}
              />
              {label}
            </label>
          {/each}
        </div>
      </section>

      <Slider
        label="Trait de coupe"
        unit="mm"
        min={0.1}
        max={1}
        step={0.05}
        value={editor.cutWidth}
        decimals={2}
        onchange={(v) => (editor.cutWidth = v)}
      />

      <section class="space-y-2">
        <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Cartouche</span>
        <div class="grid grid-cols-2 gap-2">
          {#each [["title", "Titre"], ["project", "Projet"], ["author", "Auteur"], ["date", "Date"]] as [key, ph] (key)}
            <input
              type="text"
              placeholder={ph}
              class="panel-input px-2.5 py-1.5 text-[13px] min-w-0"
              value={editor.titleBlock[key as keyof typeof editor.titleBlock]}
              oninput={(e) =>
                (editor.titleBlock = {
                  ...editor.titleBlock,
                  [key]: (e.currentTarget as HTMLInputElement).value,
                })}
            />
          {/each}
        </div>
        <input
          type="text"
          placeholder="Notes (procédé, matière, consignes…)"
          class="w-full panel-input px-2.5 py-1.5 text-[13px]"
          value={editor.titleBlock.notes}
          oninput={(e) =>
            (editor.titleBlock = {
              ...editor.titleBlock,
              notes: (e.currentTarget as HTMLInputElement).value,
            })}
        />
      </section>
    </div>
  {/if}
</div>
