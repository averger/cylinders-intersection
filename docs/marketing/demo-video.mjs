/**
 * Démo vidéo Cylix pour LinkedIn — cas cylindre / cylindre, avant / après.
 *
 * Scénario (~15 s) :
 *   1. les deux gabarits 2D AVANT modification (Ø₂ = Ø₁, 45°)
 *   2. passage en 3D, on tourne la pièce
 *   3. Ø₂ : 100 → 80 mm (0,8·Ø₁)
 *   4. angle des axes : 45° → 30°
 *   5. retour en 2D : les deux gabarits APRÈS, ils ont suivi
 * Thème clair.
 *
 * Aucun sous-titre : la démo se lit sur les valeurs du panneau et la forme
 * des gabarits.
 *
 * Notes d'enregistrement : le rendu WebGL logiciel sature le thread
 * principal et chaque événement d'entrée déclenche un recalcul, donc on
 * minimise les allers-retours CDP (mouse.move(…, {steps}) en un appel,
 * curseur factice suivi EN PAGE).  On filme en 1920×1080, la taille pour
 * laquelle l'app est dessinée, et on réduit à l'encodage.  Les marques de
 * phase sont écrites dans marks.json pour le montage.
 */
import { chromium } from "playwright-core";
import { writeFileSync } from "node:fs";

const W = 1920;
const H = 1080;

const b = await chromium.launch({
  executablePath: "/opt/pw-browsers/chromium",
  args: ["--no-sandbox", "--use-gl=swiftshader", "--enable-unsafe-swiftshader", "--hide-scrollbars"],
});
const ctx = await b.newContext({
  viewport: { width: W, height: H },
  recordVideo: { dir: "video-raw", size: { width: W, height: H } },
});
const t0 = Date.now();
const pg = await ctx.newPage();
pg.setDefaultTimeout(60000);
pg.on("pageerror", (e) => console.log("[STACK]", (e.stack || e.message).split("\n").slice(0, 5).join(" | ")));
pg.on("crash", () => console.log("[CRASH]"));

// État de départ : Cyl/Cyl, Ø₂ = Ø₁ = 100, 45°, thème clair, VUE 2D.
await pg.addInitScript(() => {
  const study = {
    id: "s1",
    name: "Étude 1",
    params: {
      mode: "cyl_cyl",
      d1: 100,
      d2: 100,
      angleDeg: 45,
      angleYDeg: 0,
      branch: "outer",
      z0: 0,
      samples: 1440,
      branches: [{ d: 60, z: 0, angleDeg: 45, azimutDeg: 0 }],
    },
    editor: {
      page: { format: "a4", orientation: "landscape", margin_mm: 10 },
      scale: "one_to_one",
      layers: { grid: true, frame: true, axis: true, labels: true, scale_bar: true },
      titleBlock: { title: "", project: "", author: "", date: "2026-08-06", notes: "" },
      annotations: [],
      cutWidth: 0.35,
    },
  };
  localStorage.setItem(
    "cylix.studies.v1",
    JSON.stringify({ studies: [study], activeId: "s1", view: "2d", theme: "light" }),
  );
});

await pg.goto("http://127.0.0.1:8787/", { waitUntil: "networkidle" });
await pg.waitForTimeout(2200); // gabarits 2D tracés

// --- curseur factice, suivi en page (Playwright ne filme pas le pointeur)
await pg.evaluate(() => {
  const cur = document.createElement("div");
  cur.style.cssText =
    "position:fixed;left:-50px;top:-50px;width:20px;height:20px;border-radius:50%;" +
    "background:rgba(224,76,13,0.30);border:2px solid #e04c0d;pointer-events:none;" +
    "z-index:99999;transform:translate(-50%,-50%)";
  document.body.appendChild(cur);
  addEventListener(
    "mousemove",
    (e) => {
      cur.style.left = e.clientX + "px";
      cur.style.top = e.clientY + "px";
    },
    true,
  );
  // La pastille « calcul… » est un état réel de l'app, mais ici le rendu
  // logiciel la fait clignoter en permanence alors qu'avec un GPU le calcul
  // est imperceptible : on la masque pour que la démo reste représentative.
  const st = document.createElement("style");
  st.textContent = ".pill.glass-strong { display: none !important }";
  document.head.appendChild(st);

});

const marks = {};
const mark = (name) => {
  marks[name] = (Date.now() - t0) / 1000;
};

/** Glisse un curseur à la souris jusqu'à `target`, puis cale la valeur exacte. */
async function dragSlider(index, target, chunks, pauseMs) {
  const el = pg.locator('input[type="range"]').nth(index);
  const box = await el.boundingBox();
  const min = parseFloat(await el.getAttribute("min"));
  const max = parseFloat(await el.getAttribute("max"));
  const value = parseFloat(await el.inputValue());
  const thumb = 15;
  const px = (v) => box.x + thumb / 2 + ((v - min) / (max - min)) * (box.width - thumb);
  const y = box.y + box.height / 2;
  await pg.mouse.move(px(value), y);
  await pg.mouse.down();
  for (let k = 1; k <= chunks; k++) {
    await pg.mouse.move(px(value + ((target - value) * k) / chunks), y, { steps: 5 });
    await pg.waitForTimeout(pauseMs);
  }
  await pg.mouse.up();
  await el.evaluate((node, v) => {
    const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, "value").set;
    setter.call(node, String(v));
    node.dispatchEvent(new Event("input", { bubbles: true }));
  }, target);
}

const clickView = async (name) => {
  const btn = pg.getByRole("button", { name, exact: true });
  const bb = await btn.boundingBox();
  await pg.mouse.move(bb.x + bb.width / 2, bb.y + bb.height / 2, { steps: 6 });
  await btn.click();
};

// --- 1. AVANT : les deux gabarits, Ø₂ = Ø₁, 45° ------------------------
mark("avant");
await pg.waitForTimeout(1600);

// --- 2. passage en 3D, on tourne la pièce -----------------------------
mark("vers3d");
await clickView("3D");
await pg.waitForTimeout(1200);
mark("orbite");
const cx = 1180;
const cy = 560;
await pg.mouse.move(cx, cy);
await pg.mouse.down();
await pg.mouse.move(cx + 140, cy - 40, { steps: 12 });
await pg.mouse.move(cx + 200, cy + 20, { steps: 10 });
await pg.mouse.up();

// --- 3. Ø₂ : 100 → 80 mm (0,8·Ø₁) -------------------------------------
mark("diametre");
await dragSlider(1, 80, 6, 110);

// --- 4. angle : 45° → 30° ---------------------------------------------
mark("angle");
await dragSlider(2, 30, 5, 110);
await pg.waitForTimeout(200);

// --- 5. APRÈS : retour 2D, les gabarits ont suivi ----------------------
mark("clic2d");
await clickView("2D");
await pg.waitForTimeout(1500);
mark("apres");
await pg.waitForTimeout(1800);
await pg.waitForTimeout(300);
mark("fin");

writeFileSync("marks.json", JSON.stringify(marks, null, 1));
console.log(JSON.stringify(marks));
await ctx.close();
await b.close();
console.log("done");
