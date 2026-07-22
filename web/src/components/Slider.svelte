<script lang="ts">
  interface Props {
    label: string;
    unit?: string;
    min: number;
    max: number;
    step?: number;
    value: number;
    onchange: (v: number) => void;
    decimals?: number;
    hint?: string;
  }

  let {
    label,
    unit = "",
    min,
    max,
    step = 1,
    value,
    onchange,
    decimals = 1,
    hint,
  }: Props = $props();

  let pct = $derived(((value - min) / (max - min)) * 100);
  let editing = $state(false);
  let buffer = $state(value.toFixed(decimals));

  function commit(v: number) {
    if (Number.isFinite(v)) {
      const clamped = Math.min(max, Math.max(min, v));
      onchange(clamped);
    }
    editing = false;
  }

  function focusInput(node: HTMLInputElement) {
    node.focus();
    node.select();
  }
</script>

<div class="space-y-2 select-none">
  <div class="flex items-baseline justify-between gap-3">
    <div class="flex flex-col">
      <span class="text-[10px] uppercase tracking-[0.18em] text-ash">{label}</span>
      {#if hint}<span class="text-[10px] text-ash/70">{hint}</span>{/if}
    </div>
    {#if editing}
      <input
        type="number"
        class="num panel-input px-2 py-1 w-28 text-right text-sm"
        bind:value={buffer}
        onblur={() => commit(parseFloat(buffer))}
        onkeydown={(e) => {
          if (e.key === "Enter") (e.currentTarget as HTMLInputElement).blur();
          if (e.key === "Escape") {
            buffer = value.toFixed(decimals);
            editing = false;
          }
        }}
        oninput={(e) => (buffer = (e.currentTarget as HTMLInputElement).value)}
        use:focusInput
      />
    {:else}
      <button
        type="button"
        class="num text-pearl text-base hover:text-ember-soft transition-colors"
        onclick={() => {
          buffer = value.toFixed(decimals);
          editing = true;
        }}
      >
        {value.toFixed(decimals)}<span class="text-ash text-xs ml-1">{unit}</span>
      </button>
    {/if}
  </div>
  <input
    type="range"
    {min}
    {max}
    {step}
    {value}
    oninput={(e) => onchange(parseFloat((e.currentTarget as HTMLInputElement).value))}
    style="--fill: {pct}%"
    class="w-full"
  />
</div>
