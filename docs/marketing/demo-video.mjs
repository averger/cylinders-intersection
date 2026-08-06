/**
 * Démo vidéo Cylix pour LinkedIn — cas cylindre / cylindre.
 *
 * Scénario (~15 s) : on tourne la 3D, on passe Ø₂ de 100 à 80 mm
 * (0,8·Ø₁), on passe l'angle de 45° à 30°, puis on bascule en 2D pour
 * montrer que les deux gabarits ont suivi. Thème clair.
 *
 * Notes d'enregistrement : le rendu WebGL logiciel sature le thread
 * principal, donc chaque appel CDP coûte cher — on minimise les
 * allers-retours (mouse.move(…, {steps}) en un appel, curseur factice
 * suivi EN PAGE par un écouteur mousemove) et on enregistre en 720p.
 */
import { chromium } from "playwright-core";

const W = 1280;
const H = 720;
const OUT = "video-raw";

const b = await chromium.launch({
  executablePath: "/opt/pw-browsers/chromium",
  args: ["--no-sandbox", "--use-gl=swiftshader", "--enable-unsafe-swiftshader", "--hide-scrollbars"],
});
const ctx = await b.newContext({
  viewport: { width: W, height: H },
  recordVideo: { dir: OUT, size: { width: W, height: H } },
});
const t0 = Date.now();
const pg = await ctx.newPage();
pg.setDefaultTimeout(60000);
pg.on("pageerror", (e) => console.log("[STACK]", (e.stack || e.message).split("\n").slice(0, 5).join(" | ")));
pg.on("crash", () => console.log("[CRASH]"));

// État de départ : Cyl/Cyl, Ø₂ = Ø₁ = 100, 45°, thème clair, vue 3D.
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
    JSON.stringify({ studies: [study], activeId: "s1", view: "3d", theme: "light" }),
  );
});

await pg.goto("http://127.0.0.1:8787/", { waitUntil: "networkidle" });
await pg.waitForTimeout(3000); // première scène 3D construite

// --- curseur factice (suivi en page) + bandeau de légende --------------
await pg.evaluate(() => {
  const cur = document.createElement("div");
  cur.style.cssText =
    "position:fixed;left:-50px;top:-50px;width:15px;height:15px;border-radius:50%;" +
    "background:rgba(224,76,13,0.30);border:1.5px solid #e04c0d;pointer-events:none;" +
    "z-index:99999;transform:translate(-50%,-50%)";
  document.body.appendChild(cur);
  // Suivi en page : zéro aller-retour pendant la chorégraphie.
  addEventListener(
    "mousemove",
    (e) => {
      cur.style.left = e.clientX + "px";
      cur.style.top = e.clientY + "px";
    },
    true,
  );
  const cap = document.createElement("div");
  cap.id = "__cap";
  cap.style.cssText =
    "position:fixed;left:50%;bottom:74px;transform:translateX(-50%);z-index:99999;" +
    "font-family:'JetBrains Mono',monospace;font-size:14px;letter-spacing:0.04em;" +
    "color:#1b1b1f;background:rgba(255,255,255,0.94);border:1px solid rgba(0,0,0,0.10);" +
    "padding:7px 16px;border-radius:999px;box-shadow:0 4px 18px rgba(20,20,25,0.10);" +
    "opacity:0;transition:opacity .3s;pointer-events:none;white-space:nowrap";
  document.body.appendChild(cap);
  window.__cap = (t) => {
    const c = document.getElementById("__cap");
    if (!c) return;
    c.style.opacity = t ? "1" : "0";
    if (t) c.textContent = t;
  };
});

const cap = (t) => pg.evaluate((t2) => window.__cap(t2), t);

/** Glisse un curseur à la souris jusqu'à `target`, puis cale la valeur exacte. */
async function dragSlider(index, target, chunks, pauseMs) {
  const el = pg.locator('input[type="range"]').nth(index);
  const box = await el.boundingBox();
  const min = parseFloat(await el.getAttribute("min"));
  const max = parseFloat(await el.getAttribute("max"));
  const value = parseFloat(await el.inputValue());
  const thumb = 15; // le pouce est centré sur la valeur
  const px = (v) => box.x + thumb / 2 + ((v - min) / (max - min)) * (box.width - thumb);
  const y = box.y + box.height / 2;
  await pg.mouse.move(px(value), y);
  await pg.mouse.down();
  for (let k = 1; k <= chunks; k++) {
    const v = value + ((target - value) * k) / chunks;
    await pg.mouse.move(px(v), y, { steps: 5 });
    await pg.waitForTimeout(pauseMs);
  }
  await pg.mouse.up();
  await el.evaluate((node, v) => {
    const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, "value").set;
    setter.call(node, String(v));
    node.dispatchEvent(new Event("input", { bubbles: true }));
  }, target);
}

const tStart = Date.now();
const t = () => ((Date.now() - tStart) / 1000).toFixed(2);

// --- 1. on joue avec la 3D : rotation auto, puis orbite à la souris
await cap("Cyl / Cyl · Ø₂ = Ø₁ · 45°");
await pg.waitForTimeout(900);
const cx = 800;
const cy = 380;
await pg.mouse.move(cx, cy);
await pg.mouse.down();
await pg.mouse.move(cx + 130, cy - 45, { steps: 14 });
await pg.mouse.move(cx + 210, cy + 25, { steps: 12 });
await pg.mouse.move(cx + 90, cy - 10, { steps: 12 });
await pg.mouse.up();
console.log("orbite", t());

// --- 2. Ø₂ : 100 → 80 mm (0,8·Ø₁)
await cap("Ø₂ : 100 → 80 mm");
await dragSlider(1, 80, 6, 120);
console.log("diamètre", t());

// --- 3. angle : 45° → 30°
await cap("angle des axes : 45° → 30°");
await dragSlider(2, 30, 5, 120);
await pg.waitForTimeout(250);
console.log("angle", t());

// --- 4. bascule 2D : les gabarits ont suivi
await cap("les deux gabarits suivent · échelle 1:1");
const twoD = pg.getByRole("button", { name: "2D", exact: true });
const bb = await twoD.boundingBox();
await pg.mouse.move(bb.x + bb.width / 2, bb.y + bb.height / 2, { steps: 8 });
await twoD.click();
await pg.waitForTimeout(1400);
await pg.mouse.move(700, 460, { steps: 4 });
await pg.mouse.wheel(0, 260);
await pg.waitForTimeout(1200);
await cap("");
await pg.waitForTimeout(250);
console.log("2D", t());

console.log("OFFSET_MS", tStart - t0, "CHOREO_MS", Date.now() - tStart);
await ctx.close();
await b.close();
console.log("done");
