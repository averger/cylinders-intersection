# Intersection de deux cylindres inclinés — théorie & mode d’emploi

*Document de référence rédigé pour expliquer la géométrie et l’usage de l’application desktop Tauri (Rust + Svelte) de génération des gabarits (développés) d’intersection de deux cylindres. Le script Python original reste fourni en référence.*

— **A. Verger**

---

## 1) Théorie (version courte mais complète)

### Modèle géométrique

On travaille dans un repère orthonormé $(x,y,z)$.

**Cylindre 1** (*vertical*), axe $Oz$, rayon $R_1>0$ :

$$
x^2 + y^2 = R_1^2.
$$

**Cylindre 2** (*incliné*), rayon $R_2>0$, obtenu par rotation d’angle $\phi$ autour de l’axe $x$ d’un cylindre initialement coaxial à $Oz$.  
Paramétrisation **avant** rotation : $(R_2\cos\theta,\ R_2\sin\theta,\ t)$.  
Après la rotation $R_x(\phi)$ (autour de $x$), on obtient :

$$
\begin{aligned}
x(\theta,t) &= R_2\cos\theta,\\
y(\theta,t) &= R_2\sin\theta\cos\phi - t\sin\phi,\\
z(\theta,t) &= R_2\sin\theta\sin\phi + t\cos\phi.
\end{aligned}
$$

> *Remarque.* Pour $\phi=0$, le cylindre 2 redevient coaxial à $Oz$. Dans ce cas particulier, l’équation quadratique dégénère ($a=0$) et l’intersection est vide si $R_1\ne R_2$, ou bien un cercle de rayon $R_1$ si $R_1=R_2$.

### Cas particulier : cylindre × plan incliné

Quand le « cylindre 2 » dégénère vers un plan, on remplace l’équation $x^2+y^2=R_2^2$ par celle d’un plan d’inclinaison $\phi$ autour de $Ox$, passant par $(0,0,z_0)$ :

$$
y\sin\phi + z\cos\phi = z_0\cos\phi
\qquad\Longleftrightarrow\qquad
z(\theta) = z_0 - R_1\sin\theta\,\tan\phi.
$$

L’intersection est analytique : pour $\theta\in[0,2\pi)$, le point sur le cylindre 1 vaut $(R_1\cos\theta, R_1\sin\theta, z(\theta))$ et le développé sur le tube vertical est directement $(u,v) = (R_1\,\theta,\ z(\theta))$ — c’est la sinusoïde classique d’une **coupe d’onglet**. Aucun choix de branche n’est nécessaire.

### Condition d’intersection et équation en $t$

Un point du cylindre 2 appartient au cylindre 1 ssi $x^2(\theta,t)+y^2(\theta,t)=R_1^2$.  
En remplaçant $x$ et $y$ par les expressions ci-dessus, on obtient, pour chaque $\theta\in[0,2\pi)$, un **quadratique en $t$** :

$$
a\,t^2 + b(\theta)\,t + c_0(\theta)=0,
$$

avec

$$
\begin{aligned}
a &= \sin^2\phi,\\
b(\theta) &= -\,2\,R_2\,\sin\theta\,\cos\phi\,\sin\phi,\\
c_0(\theta) &= R_2^2\!\big(\cos^2\theta+\sin^2\theta\,\cos^2\phi\big)-R_1^2.
\end{aligned}
$$

Le **discriminant** $\Delta(\theta)=b^2(\theta)-4\,a\,c_0(\theta)$ décide de l’existence de solutions réelles. Dès qu’il existe des $\theta$ tels que $\Delta(\theta)\ge 0$, les deux surfaces se coupent. Les solutions sont

$$
t_\pm(\theta)=\frac{-\,b(\theta)\pm\sqrt{\Delta(\theta)}}{2a}\qquad(\phi\ne 0).
$$

### Deux branches (« outer » / « inner »)

Pour une génératrice donnée ($\theta$ fixé), il y a en général **deux points** d’intersection — « entrée »/« sortie ». En balayant $\theta$, cela engendre **deux courbes**.  
Dans le code, on choisit la branche via `branch="outer"` (racine $+$) ou `branch="inner"` (racine $-$). **En fabrication**, on n’utilise généralement **qu’une seule lèvre** (souvent `outer`).

### Développés (mises à plat)

Le développement isométrique d’un cylindre de rayon $R$ vers le plan se fait via la carte $(\theta,v)\mapsto(u,v)$ avec

$$
u = R\,\theta,\qquad v = \text{coordonnée axiale (inchangée)}.
$$

On en déduit :

**(i) Gabarit sur le cylindre incliné (tube à couper)**

$$
\big(u_2(\theta),\,v_2(\theta)\big)=\big(R_2\,\theta,\ t_\star(\theta)\big),
$$

où $t_\star$ désigne l’une des deux branches $t_\pm$.

**(ii) « Gueule de loup » sur le cylindre vertical**

On passe en cylindriques du cylindre 1 : $\alpha=\mathrm{atan2}(y,x)$ et $z=z$, puis

$$
(u_1,v_1)=\big(R_1\,\alpha^\uparrow,\ z\big),
$$

où $\alpha^\uparrow$ est l’angle **déroulé** (*unwrap*) pour supprimer le saut à $2\pi$.

> Les deux courbes 2D, **à l’échelle 1**, servent de gabarits : l’une pour **découper** le tube incliné, l’autre pour **présenter** le tube vertical (contact « fish-mouth »).

---

## 2) Application desktop Tauri + Rust + Svelte (gabarits 1 : 1)

L’implémentation de référence est désormais une **application desktop** packagée avec **Tauri 2** — moteur géométrique **Rust + nalgebra**, interface **Svelte 5 + Tailwind 4**, vue 3D **Three.js**. Plus de serveur HTTP : les calculs sont exposés directement sous forme de commandes `invoke` consommées par le frontend.

* Calcul matriciel côté Rust (commande Tauri, zéro réseau).
* Choix entre **cylindre × cylindre** et **cylindre × plan incliné**.
* Export **SVG vectoriel à l’échelle 1 : 1** et **CSV** brut pour CNC / laser / découpe plasma.
* Mode impression dédié, avec règle de référence 100 mm pour vérifier le scaling de l’imprimante.
* Bundles natifs Windows (`.msi`, `.exe`), macOS (`.dmg`, `.app`) et Linux (`.deb`, `.AppImage`, `.rpm`).

### 2.1 Pré-requis

| Outil       | Version testée |
|-------------|----------------|
| Rust        | 1.94 +         |
| Node        | 22 +           |
| npm         | 10 +           |
| Tauri CLI   | 2.x (installé via `npm i`) |

Dépendances système Linux (Ubuntu/Debian) : `libwebkit2gtk-4.1-dev`, `libgtk-3-dev`, `libayatana-appindicator3-dev`, `librsvg2-dev`, `libsoup-3.0-dev`, `pkg-config`. Sur macOS : Xcode Command Line Tools. Sur Windows : Microsoft C++ Build Tools + WebView2 (préinstallé sur Windows 11).

### 2.2 Build complet

```bash
# 1) installer les dépendances frontend (inclut le CLI Tauri)
cd web && npm install && cd ..

# 2) compiler l'app desktop (déclenche automatiquement npm run build sur la SPA)
cd web && npm run tauri:build
# bundles produits dans src-tauri/target/release/bundle/
```

Le binaire seul (sans installateur) reste disponible dans `src-tauri/target/release/`.

### 2.3 Développement

```bash
cd web && npm run tauri:dev
# Vite démarre sur http://127.0.0.1:5173 (hot-reload Svelte)
# Tauri lance simultanément la fenêtre native et y branche le dev server
```

> Note : `tauri:dev` orchestre déjà le frontend Vite via `beforeDevCommand` (cf. `src-tauri/tauri.conf.json`) — pas besoin de second terminal.

Pour tester uniquement le frontend dans le navigateur (sans Tauri), `npm run dev` reste utilisable, mais les appels `invoke` échoueront hors du runtime Tauri.

### 2.4 Commandes Tauri

L’ancienne API HTTP (Axum) a été remplacée par deux **commandes Tauri** invoquées depuis le frontend via `@tauri-apps/api`.

| Commande Rust            | Argument             | Réponse              |
|--------------------------|----------------------|----------------------|
| `intersect_cyl_cyl`      | `CylCylInput`        | `IntersectionPayload` |
| `intersect_cyl_plane`    | `CylPlaneInput`      | `IntersectionPayload` |

```ts
import { invoke } from "@tauri-apps/api/core";

const res = await invoke("intersect_cyl_cyl", {
  input: { r1: 50, r2: 35, phi: Math.PI / 4, n_samples: 1440, branch: "outer" },
});
```

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

* Toutes les longueurs sont en **millimètres**.
* Les angles sont en **radians**.
* `branch ∈ {"outer", "inner"}` choisit la racine $t_\pm$ (cf. §1).
* `dev_main` n’est renvoyé qu’en mode `cyl_cyl`.
* Les erreurs de validation reviennent sous forme de `Promise.reject(string)` — interceptées par `web/src/lib/api.ts`.

### 2.5 Impression 1 : 1

* Le bouton **« Imprimer 1 : 1 »** déclenche `window.print()`. Le CSS d’impression masque l’UI et ne laisse visibles que les gabarits sur fond blanc.
* Chaque SVG est dimensionné en `mm` (`width="…mm"`, `viewBox` identique). Le pilote d’impression doit être réglé sur **« Échelle réelle »** / **« 100 % »** / **« Actual size »** ; désactiver l’option « ajuster à la page ».
* Une **règle de 100 mm** est imprimée en bas de chaque feuille pour valider l’étalonnage après impression. Si elle ne mesure pas exactement 100 mm une fois sur papier, ajuster le pilote.

Exports disponibles depuis chaque carte :

| Format | Usage                                                            |
|--------|------------------------------------------------------------------|
| `SVG`  | LightBurn, RDWorks, Inkscape, Illustrator, Fusion 360, Onshape — vectoriel pur |
| `CSV`  | Colonnes `theta_rad, u_mm, v_mm` — intégration dans un post-processeur CNC |

### 2.6 Arborescence

```
cylinders-intersection/
├── src-tauri/                  # crate Tauri (anciennement Axum)
│   ├── Cargo.toml
│   ├── tauri.conf.json         # config bundle, fenêtre, beforeDev/Build
│   ├── build.rs
│   ├── capabilities/
│   │   └── default.json        # permissions de la fenêtre principale
│   ├── icons/                  # icônes (PNG/ICO/ICNS)
│   └── src/
│       ├── geometry.rs         # rotation X, paramétrisation cyl-2, BBox, unwrap
│       ├── intersection.rs     # cyl/cyl + cyl/plan + développés (gueule de loup)
│       ├── lib.rs              # commandes Tauri + run()
│       └── main.rs             # entrée binaire native
└── web/                        # SPA (Svelte 5 + Tailwind 4 + Three.js)
    ├── package.json            # scripts tauri:dev / tauri:build
    ├── vite.config.ts
    ├── tsconfig.json
    ├── index.html
    └── src/
        ├── App.svelte
        ├── main.ts
        ├── app.css
        ├── components/
        │   ├── Header.svelte
        │   ├── Hero.svelte
        │   ├── ControlPanel.svelte
        │   ├── Viewer3D.svelte         # rendu Three.js
        │   ├── DevelopedView.svelte    # SVG mm interactif
        │   ├── PrintLayout.svelte      # rendu impression 1:1
        │   ├── Slider.svelte
        │   └── Segmented.svelte
        └── lib/
            ├── api.ts          # wrapper invoke() Tauri
            ├── store.svelte.ts
            └── svg.ts
```

### 2.7 Validation numérique

Tests unitaires Rust :

```bash
cd src-tauri && cargo test --lib --release
```

Deux cas couverts par défaut :

* **Cylindres égaux à 90°** — la racine `outer` doit atteindre le maximum $v=R_1$ ($\pm 10^{-6}$).
* **Plan à 45° sur un tube Ø 60 mm** — l’abscisse couvre la totalité de $2\pi R_1$ (à un pas d’échantillonnage près).

### 2.8 Script Python historique (legacy)

Le script de référence `tubes_intersection_and_unwrap.py` reste fourni pour comparaison ou usage en notebook. Il implémente la même algèbre que le backend Rust (cas cyl/cyl uniquement). Dépendances : `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib
python tubes_intersection_and_unwrap.py     # exemple d'exécution intégré
```

---

## 3) Notes de validité et cas limites

* Pour $\phi=0$, l’équation quadratique **dégénère** ($a=0$) : on retombe sur deux cylindres coaxiaux. L’intersection est un cercle seulement si $R_1=R_2$. Le backend renvoie alors un payload vide (sans erreur) ; l’UI affiche « Pas d’intersection valide ».
* Pour le mode plan, $\phi$ doit rester strictement dans $(-\pi/2,\,\pi/2)$ : à $\pm\pi/2$ le plan devient parallèle à l’axe du cylindre. La commande Tauri rejette ce cas avec une `Promise.reject(string)` que l’UI affiche dans la bannière d’erreur.
* Les formules de $a$, $b(\theta)$, $c_0(\theta)$ ci-dessus proviennent de l’égalité $x^2+y^2=R_1^2$ avec
  $x=R_2\cos\theta$, $y=R_2\sin\theta\cos\phi - t\sin\phi$. On a bien
  $$x^2+y^2=R_2^2\cos^2\theta+R_2^2\sin^2\theta\cos^2\phi-2R_2\sin\theta\cos\phi\sin\phi\,t+\sin^2\phi\,t^2,$$
  d’où le quadratique annoncé.
* La **branche** (`outer` / `inner`) sélectionne la lèvre supérieure ou inférieure de l’intersection. En tuyauterie la « gueule de loup » conserve presque toujours `outer` ; `inner` sert pour les pénétrations à recouvrement.
* Le développé sur le cylindre principal utilise un `unwrap` analogue à `numpy.unwrap` pour supprimer le saut de $\pm\pi$ de `atan2` et obtenir un tracé continu sur $[-\pi R_1, \pi R_1]$.

---

*Document de référence + manuel de la webapp.*

— **A. Verger**
