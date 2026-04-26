<script lang="ts" generics="T extends string">
  interface Option {
    label: string;
    value: T;
  }
  interface Props {
    options: Option[];
    value: T;
    onchange: (v: T) => void;
  }
  let { options, value, onchange }: Props = $props();

  let container = $state<HTMLDivElement | undefined>(undefined);
  let thumbStyle = $state("");

  function updateThumb() {
    if (!container) return;
    const idx = options.findIndex((o) => o.value === value);
    const buttons = container.querySelectorAll<HTMLButtonElement>("button[data-seg]");
    const btn = buttons[idx];
    if (!btn) return;
    const containerRect = container.getBoundingClientRect();
    const r = btn.getBoundingClientRect();
    thumbStyle = `left: ${r.left - containerRect.left}px; width: ${r.width}px;`;
  }

  $effect(() => {
    // depend on value & options to re-measure
    void value;
    void options;
    requestAnimationFrame(updateThumb);
  });

  $effect(() => {
    const obs = new ResizeObserver(updateThumb);
    if (container) obs.observe(container);
    return () => obs.disconnect();
  });
</script>

<div class="seg" bind:this={container}>
  <span class="seg-thumb" style={thumbStyle}></span>
  {#each options as opt (opt.value)}
    <button
      data-seg
      type="button"
      aria-pressed={value === opt.value}
      onclick={() => onchange(opt.value)}
    >
      {opt.label}
    </button>
  {/each}
</div>
