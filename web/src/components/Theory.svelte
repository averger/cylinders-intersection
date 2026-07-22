<script lang="ts">
  import katex from "katex";
  import "katex/dist/katex.min.css";

  // Render KaTeX at mount time via an action — content is static.
  function math(node: HTMLElement, tex: string) {
    katex.render(tex, node, { displayMode: true, throwOnError: false });
  }
  function imath(node: HTMLElement, tex: string) {
    katex.render(tex, node, { displayMode: false, throwOnError: false });
  }

  interface Block {
    title: string;
    tag: string;
    body: string;
    tex: string;
    note?: string;
  }

  const blocks: Block[] = [
    {
      title: "Paramétrisation",
      tag: "01",
      body:
        "Le cylindre principal (rayon R₁) a son axe sur Oz. Le cylindre secondaire " +
        "(rayon R₂) est un cylindre vertical basculé d'un angle φ autour de Ox — " +
        "sa génératrice d'angle θ reste dans le plan x = R₂ cos θ.",
      tex: String.raw`\begin{aligned}
x &= R_2\cos\theta\\
y &= R_2\sin\theta\,\cos\varphi - t\,\sin\varphi\\
z &= R_2\sin\theta\,\sin\varphi + t\,\cos\varphi
\end{aligned}`,
    },
    {
      title: "L'intersection est un second degré",
      tag: "02",
      body:
        "Imposer x² + y² = R₁² donne, à θ fixé, une équation quadratique en t " +
        "dont le discriminant se réduit à une forme remarquablement simple :",
      tex: String.raw`t_\pm(\theta)=\frac{R_2\sin\theta\cos\varphi\;\pm\;\sqrt{R_1^{2}-R_2^{2}\cos^{2}\theta}}{\sin\varphi}`,
      note:
        "Condition d'existence : R₂|cos θ| ≤ R₁ — la génératrice doit atteindre le " +
        "gros tube. Les racines ± sont l'entrée et la sortie de la paroi (branches inner/outer).",
    },
    {
      title: "Développer = isométrie",
      tag: "03",
      body:
        "Un cylindre a une courbure de Gauss nulle : il se déroule sur le plan sans " +
        "distorsion. Le déroulage préserve exactement longueurs et angles — c'est ce " +
        "qui autorise le gabarit papier à l'échelle 1:1.",
      tex: String.raw`\mathrm{d}s^{2}=R^{2}\,\mathrm{d}\theta^{2}+\mathrm{d}t^{2}
\;\;\xrightarrow{\;(u,v)=(R\theta,\;t)\;}\;\;
\mathrm{d}u^{2}+\mathrm{d}v^{2}`,
    },
    {
      title: "Gueule de loup",
      tag: "04",
      body:
        "Vue du cylindre principal, la même courbe se déroule avec (u₁, v₁) = (R₁α, z). " +
        "Les deux angles sont liés par une relation de transfert exacte ; à 90°, on " +
        "retrouve la formule traditionnelle de traçage :",
      tex: String.raw`\cos\alpha=\frac{R_2}{R_1}\cos\theta
\qquad\Longrightarrow\qquad
v_1(\alpha)\Big|_{\varphi=\pi/2}=\pm\sqrt{R_2^{2}-R_1^{2}\cos^{2}\alpha}`,
    },
    {
      title: "Coupe en sifflet",
      tag: "05",
      body:
        "Un tube coupé par un plan incliné donne une ellipse dans l'espace… et une " +
        "sinusoïde parfaite une fois déroulé — période = périmètre, amplitude = R tan φ.",
      tex: String.raw`v(u)=z_0-R_1\tan\varphi\,\sin\!\frac{u}{R_1}`,
    },
  ];
</script>

<section id="theorie" class="max-w-[1600px] mx-auto px-6 lg:px-10 pb-20 scroll-mt-24">
  <div class="flex flex-col lg:flex-row lg:items-end lg:justify-between gap-3 mb-8">
    <div>
      <span class="pill">théorie</span>
      <h2 class="text-pearl text-2xl lg:text-3xl font-semibold tracking-tight mt-2">
        La géométrie exacte, sans approximation
      </h2>
      <p class="text-ash text-sm max-w-2xl mt-1">
        Le moteur ne discrétise jamais les équations&nbsp;: chaque point provient d'une
        formule fermée. Dérivation complète, cas limites et bornes d'erreur dans
        <a
          href="https://github.com/averger/cylinders-intersection/blob/main/docs/THEORY.md"
          target="_blank"
          rel="noopener"
          class="text-cyan-glow hover:text-cyan-soft transition-colors underline underline-offset-2"
          >docs/THEORY.md</a
        >.
      </p>
    </div>
    <div class="text-[11px] text-ash num">
      5 résultats clés · démonstrations vérifiées par tests
    </div>
  </div>

  <div class="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-5">
    {#each blocks as b (b.tag)}
      <article class="glass card-hairline p-6 flex flex-col gap-4 min-w-0">
        <header class="flex items-start justify-between gap-3">
          <h3 class="text-pearl text-base font-semibold tracking-tight">{b.title}</h3>
          <span class="num text-[11px] text-ash/70">{b.tag}</span>
        </header>
        <p class="text-[13px] leading-relaxed text-silver">{b.body}</p>
        <div class="theory-math rounded-xl bg-black/35 border border-white/5 px-4 py-3 overflow-x-auto">
          <div use:math={b.tex}></div>
        </div>
        {#if b.note}
          <p class="text-[11px] leading-relaxed text-ash">{b.note}</p>
        {/if}
      </article>
    {/each}

    <article class="glass card-hairline p-6 flex flex-col gap-4 justify-between min-w-0">
      <div class="space-y-4">
        <header class="flex items-start justify-between gap-3">
          <h3 class="text-pearl text-base font-semibold tracking-tight">
            Pourquoi c'est irréprochable
          </h3>
          <span class="num text-[11px] text-ash/70">∎</span>
        </header>
        <ul class="text-[13px] leading-relaxed text-silver space-y-2 list-none">
          <li>
            <span class="text-ember mr-1.5">→</span>Discriminant en forme fermée&nbsp;:
            <span use:imath={String.raw`\Delta=4\sin^{2}\varphi\,(R_1^{2}-R_2^{2}\cos^{2}\theta)`}
            ></span>
          </li>
          <li>
            <span class="text-ember mr-1.5">→</span>Seule approximation&nbsp;: la polyligne
            d'échantillonnage, erreur bornée par la flèche
            <span use:imath={String.raw`e\le\tfrac{\kappa}{8}\ell^{2}`}></span>
            — &lt;&nbsp;0,01&nbsp;mm au réglage par défaut.
          </li>
          <li>
            <span class="text-ember mr-1.5">→</span>Cas dégénérés traités&nbsp;: Steinmetz
            (R₁&nbsp;=&nbsp;R₂), tangences Δ&nbsp;=&nbsp;0, cylindres coaxiaux rejetés
            proprement.
          </li>
          <li>
            <span class="text-ember mr-1.5">→</span>Vérifié par tests Rust et
            non-régression contre l'implémentation NumPy d'origine.
          </li>
        </ul>
      </div>
      <a
        href="https://github.com/averger/cylinders-intersection/blob/main/docs/THEORY.md"
        target="_blank"
        rel="noopener"
        class="btn-ghost justify-center"
      >
        Lire la dérivation complète
      </a>
    </article>
  </div>
</section>

<style>
  .theory-math :global(.katex-display) {
    margin: 0;
  }
  .theory-math :global(.katex) {
    color: var(--color-pearl);
    font-size: 1.02rem;
  }
</style>
