# `cylinders-intersection` — webapp Rust + Svelte

Interface graphique pour le calcul d’intersection de deux cylindres (ou d’un
cylindre et d’un plan incliné), avec génération des **développés** à plat
(« gueule de loup ») exportables à l’**échelle 1 : 1** pour traçage / découpe
laser / CNC.

> Le moteur géométrique vit en Rust ; l’interface est une SPA Svelte 5
> compilée par Vite et embarquée directement dans le binaire. Tout tourne
> en local — aucune donnée ne quitte la machine.

## Architecture

```
cylinders-intersection/
├── Cargo.toml              # crate `cylinders-intersection`
├── src/
│   ├── geometry.rs         # rotation X, paramétrisation cyl-2, BBox, unwrap
│   ├── intersection.rs     # cyl/cyl + cyl/plan + développés (gueule de loup)
│   ├── api.rs              # routes Axum + SPA fallback
│   ├── assets.rs           # rust-embed sur web/dist
│   ├── lib.rs              # re-exports
│   └── main.rs             # serveur Axum (CLI clap)
└── web/                    # SPA (Svelte 5 + Tailwind 4 + Three.js)
    ├── src/
    │   ├── main.ts
    │   ├── App.svelte
    │   ├── app.css
    │   ├── components/
    │   │   ├── Header.svelte
    │   │   ├── Hero.svelte
    │   │   ├── ControlPanel.svelte
    │   │   ├── Viewer3D.svelte         # rendu Three.js
    │   │   ├── DevelopedView.svelte    # SVG mm interactif
    │   │   ├── PrintLayout.svelte      # rendu impression 1:1
    │   │   ├── Slider.svelte
    │   │   └── Segmented.svelte
    │   └── lib/
    │       ├── api.ts          # client typé pour /api/intersect/*
    │       ├── store.svelte.ts # state Svelte 5 (runes)
    │       └── svg.ts          # générateur SVG mm + downloads
    ├── vite.config.ts
    ├── tsconfig.json
    └── index.html
```

## Build

```bash
# 1) compiler la SPA (produit web/dist)
cd web && npm install && npm run build && cd ..

# 2) compiler le binaire Rust qui embarque la SPA
cargo build --release

# 3) lancer
./target/release/cylinders-intersection --port 8787
# puis ouvrir http://127.0.0.1:8787
```

## Développement

```bash
# terminal 1 — backend Rust avec rebuild automatique
cargo run --release -- --port 8787

# terminal 2 — frontend en hot-reload, proxy /api -> :8787
cd web && npm run dev
# Vite démarre sur http://127.0.0.1:5173
```

## API

| Méthode | Route | Corps | Réponse |
|---------|-------|-------|---------|
| `GET`   | `/api/health` | — | `{name, version, uptime_ms}` |
| `POST`  | `/api/intersect/cyl-cyl` | `{r1, r2, phi, n_samples?, branch?}` | `IntersectionPayload` |
| `POST`  | `/api/intersect/cyl-plane` | `{r1, phi, z0?, n_samples?}` | `IntersectionPayload` |

`IntersectionPayload` :

```jsonc
{
  "mode": "cyl_cyl" | "cyl_plane",
  "r1": 50.0, "r2": 35.0, "phi": 0.7853, "branch": "outer",
  "curve3d":     [{ "x": ..., "y": ..., "z": ... }, ...],
  "dev_branch":  [{ "theta": ..., "u": ..., "v": ... }, ...],
  "dev_main":    [...] | null,
  "bbox_branch": { "u_min":..., "u_max":..., "v_min":..., "v_max":... },
  "bbox_main":   {...} | null,
  "circumference_branch": 219.91,
  "circumference_main":   314.16
}
```

Toutes les longueurs sont en **millimètres** ; les angles en **radians**.

## Modes de calcul

### Cylindre × cylindre

Le cylindre 1 est sur l’axe Oz, rayon `r1`. Le cylindre 2, rayon `r2`, est
obtenu par rotation `Rx(phi)` autour de Ox d’un cylindre coaxial à Oz. On
résout le quadratique en `t` :

> `sin²φ · t² − 2·r2·sinθ·cosφ·sinφ · t + (r2²(cos²θ + sin²θ·cos²φ) − r1²) = 0`

et on retient la racine selon `branch ∈ {outer, inner}`. Le développé sur le
cylindre 2 est `(u, v) = (r2·θ, t)` ; le développé sur le cylindre 1 est
`(u, v) = (r1·α_unwrap, z)`.

### Cylindre × plan incliné

Plan d’équation `z = z0 − y·tan(φ)`. Pour chaque `θ`, l’intersection est
analytique :
`x = r1·cosθ`, `y = r1·sinθ`, `z = z0 − r1·sinθ·tanφ`. Le développé est
directement `(u, v) = (r1·θ, z)`.

## Impression 1 : 1

Le bouton **Imprimer 1:1** déclenche `window.print()`. Le CSS d’impression
masque l’UI et ne laisse visibles que les `.print-area` rendus par
`PrintLayout.svelte`. Chaque SVG est dimensionné en `mm` avec `viewBox` =
mêmes valeurs en mm — le pilote d’impression doit être réglé sur **« Échelle
réelle »** / **« 100 % »** pour reproduire le gabarit au scale exact.

Une **règle de référence de 100 mm** est imprimée en bas de chaque feuille
pour vérifier l’étalonnage après impression.

Exports disponibles depuis chaque carte :

* `SVG 1:1` — vectoriel, prêt pour LightBurn / RDWorks / Illustrator / Inkscape.
* `CSV` — colonnes `theta_rad, u_mm, v_mm` pour intégration dans un post-processeur CNC.
