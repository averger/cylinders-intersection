<script lang="ts">
  // Apple-style on/off switch.
  interface Props {
    checked: boolean;
    onchange: (v: boolean) => void;
    label?: string;
  }
  let { checked, onchange, label }: Props = $props();
</script>

<button
  type="button"
  role="switch"
  aria-checked={checked}
  aria-label={label}
  class="sw {checked ? 'on' : ''}"
  onclick={() => onchange(!checked)}
>
  <span class="knob"></span>
</button>

<style>
  .sw {
    position: relative;
    width: 42px;
    height: 26px;
    flex: 0 0 auto;
    border-radius: 999px;
    background: var(--seg-bg);
    border: 1px solid var(--line-2);
    transition: background 220ms cubic-bezier(0.22, 1, 0.36, 1), border-color 220ms ease;
    cursor: pointer;
  }
  .sw.on {
    background: linear-gradient(180deg, var(--ember-soft), var(--ember));
    border-color: transparent;
  }
  .knob {
    position: absolute;
    top: 2px;
    left: 2px;
    width: 20px;
    height: 20px;
    border-radius: 999px;
    background: #fff;
    box-shadow:
      0 1px 3px rgba(20, 20, 25, 0.35),
      0 0 0 0.5px rgba(20, 20, 25, 0.08);
    transition: transform 220ms cubic-bezier(0.22, 1, 0.36, 1);
  }
  .sw.on .knob {
    transform: translateX(16px);
  }
  .sw:active .knob {
    transform: scaleX(1.1) translateX(0);
    transform-origin: left;
  }
  .sw.on:active .knob {
    transform: scaleX(1.1) translateX(14px);
    transform-origin: right;
  }
</style>
