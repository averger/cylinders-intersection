<script lang="ts">
  // Export settings — shown in the sidebar when the 2D workspace is active.
  import Segmented from "./Segmented.svelte";
  import Slider from "./Slider.svelte";
  import { store } from "../lib/store.svelte";
  import { editor } from "../lib/editor.svelte";
  import type { PageFormat, PageOrientation, ScaleMode } from "../lib/export";

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
</script>

<div class="flex flex-col gap-6">
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
    <Segmented options={scaleOptions} value={editor.scale} onchange={(v) => (editor.scale = v)} />
    <p class="text-[10px] text-ash/70 leading-relaxed">
      1:1 = tracé direct sur tôle/tube après collage des tuiles (repères inclus). « Ajusté »
      n'est qu'un aperçu — jamais pour le traçage.
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

  <section class="space-y-2">
    <span class="text-[10px] uppercase tracking-[0.18em] text-ash">Annotations</span>
    {#if editor.annotations.length === 0}
      <p class="text-[11px] text-ash/70">
        Double-cliquez sur une mise à plat pour placer un repère texte.
      </p>
    {:else}
      <ul class="space-y-1.5">
        {#each editor.annotations as a, index (index)}
          <li class="flex items-center gap-1.5">
            <span
              class="shrink-0 w-1.5 h-1.5 rounded-full"
              style="background: {a.pattern === 'branch' ? 'var(--ember)' : 'var(--cyan)'}"
              title={a.pattern === "branch" ? "tube incliné" : "gueule de loup"}
            ></span>
            <input
              type="text"
              class="flex-1 min-w-0 panel-input px-2 py-1 text-[12px] {editor.selected ===
              index
                ? 'border-ember/80'
                : 'border-mist/50 focus:border-ember/70'}"
              value={a.text}
              onfocus={() => (editor.selected = index)}
              oninput={(e) => {
                editor.annotations[index].text = (e.currentTarget as HTMLInputElement).value;
                store.persist();
              }}
            />
            <span class="num text-[9px] text-ash shrink-0">{a.u.toFixed(0)};{a.v.toFixed(0)}</span>
            <button
              class="text-ash hover:text-ember transition-colors text-sm leading-none px-0.5"
              aria-label="Supprimer l'annotation"
              onclick={() => editor.removeAnnotation(index)}
            >
              ×
            </button>
          </li>
        {/each}
      </ul>
    {/if}
  </section>

  <section class="pt-3 hl-t space-y-2">
    <div class="flex gap-2">
      <button
        class="btn-primary flex-1 !py-2.5 text-sm"
        disabled={editor.busy !== null || !store.result}
        onclick={() => editor.exportPdf()}
      >
        {editor.busy === "pdf" ? "génération…" : "PDF 1:1"}
      </button>
      <button
        class="btn-ghost flex-1 justify-center !py-2.5 text-sm"
        disabled={editor.busy !== null || !store.result}
        onclick={() => editor.exportDxf()}
      >
        {editor.busy === "dxf" ? "génération…" : "DXF"}
      </button>
    </div>
    {#if editor.error}
      <p class="text-[12px] text-ember">{editor.error}</p>
    {/if}
    <p class="text-[10px] text-ash/70 leading-relaxed">
      Le PDF contient toutes les mises à plat ({editor.available.length}) avec cartouche et
      règle de contrôle. DXF R12 en mm — calques CUT / FRAME / AXIS / TEXT / ANNOT.
    </p>
  </section>
</div>
