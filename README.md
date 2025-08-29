# Intersection de deux cylindres inclinés — théorie & mode d’emploi

*Document de référence rédigé pour expliquer la géométrie et l’usage du script Python de génération des gabarits (développés) d’intersection de deux cylindres.*

— **A. Verger**

---

## 1) Théorie (version courte mais complète)

### Modèle géométrique
- Repère orthonormé $(x,y,z)$.
- **Cylindre 1** (*vertical*), axe $Oz$, rayon $R_1>0$ :  
  
    $$x^2 + y^2 = R_1^2.$$
- **Cylindre 2** (*incliné*), rayon $R_2>0$, obtenu par rotation d’angle $\phi$ autour de l’axe $x$ d’un cylindre initialement coaxial à $Oz$.  
  Paramétrisation **avant** rotation : $(R_2\cos\theta,\ R_2\sin\theta,\ t)$.  
  **Après** rotation $R_x(\phi)$ :
  
    $$
    \begin{aligned}
    x(\theta,t) &= R_2\cos\theta,\\
    y(\theta,t) &= R_2\sin\theta\cos\phi - t\sin\phi,\\
    z(\theta,t) &= R_2\sin\theta\sin\phi + t\cos\phi.
    \end{aligned}
    $$

### Condition d’intersection et équation en $t$
Un point du cylindre 2 appartient au cylindre 1 ssi $x^2(\theta,t)+y^2(\theta,t)=R_1^2$.

On obtient, pour chaque $\theta\in[0,2\pi)$, un **quadratique en $t$** :
  
$$
a\,t^2 + b(\theta)\,t + c_0(\theta)=0.
$$

avec
  
$$
\begin{aligned}
a &= \sin^2\phi,\\
b(\theta) &= -\,2\,R_2\,\sin\theta\,\cos\phi\,\sin\phi,\\
c_0(\theta) &= R_2^2\!\left(\cos^2\theta+\sin^2\theta\,\cos^2\phi\right)-R_1^2.
\end{aligned}
$$

Le **discriminant** $\Delta(\theta)=b^2-4ac_0$ décide de l’existence de solutions réelles. Dès qu’il existe des $\theta$ avec $\Delta\ge 0$, les deux surfaces se coupent.

Les solutions sont
  
$$
t_\pm(\theta)=\frac{-\,b(\theta)\pm\sqrt{\Delta(\theta)}}{2a}\qquad(\phi\not\equiv 0).
$$

### Deux branches (« outer » / « inner »)
Pour une génératrice donnée ($\theta$ fixé), il y a en général **deux points** d’intersection — « entrée »/« sortie ». En balayant $\theta$, cela engendre **deux courbes**.  
Dans le code, je choisis l’une des deux via `branch="outer"` (racine $+$) ou `branch="inner"` (racine $-$). **En fabrication**, on n’utilise généralement **qu’une seule lèvre** (souvent `outer`).

### Développés (mises à plat)
Le développement isométrique d’un cylindre de rayon $R$ se fait via $(u,v)=(R\theta,\ v)$. On en déduit :

- **Gabarit sur le cylindre incliné (tube à couper)**  
  
    $$
    (u_2(\theta),v_2(\theta))=\big(R_2\,\theta,\ t_\star(\theta)\big).
    $$

- **« Gueule de loup » sur le cylindre vertical**  
  On passe en cylindriques du cylindre 1 : $\alpha=\mathrm{atan2}(y,x)$, $z=z$, puis  
  
    $$
    (u_1,v_1)=\big(R_1\,\alpha^\uparrow,\ z\big),
    $$
  
  où $\alpha^\uparrow$ signifie *unwrap* (suppression du saut à $2\pi$).

> Ces deux courbes 2D, **à l’échelle 1**, servent de gabarits : l’une pour **découper** le tube incliné, l’autre pour **présenter** le tube vertical (contact type « fish-mouth »).

---

## 2) Utilisation du programme

### Installation
- Dépendances Python : `numpy`, `matplotlib`.  
  ```bash
  pip install numpy matplotlib
  ```

### Script
Le fichier s’appelle par exemple `tubes_intersection_and_unwrap_v2.py`. Les points d’entrée sont :
- `intersect_cylinders(R1, R2, phi, n_samples=2000, branch="outer")`
- `plot_results(res, R1, R2, phi, **options)`
- `export_csv_cyl2(res, "gabarit_cylindre2.csv")`
- `export_csv_both(res, R1, "gabarits_developpes.csv")`

### Paramètres principaux
- `R1` : rayon du **cylindre 1** (vertical).  
- `R2` : rayon du **cylindre 2** (incliné).  
- `phi` : **angle d’inclinaison** (en radians) autour de **x**. Exemple : `np.deg2rad(30)`.  
- `branch` : `"outer"` (racine $+$) ou `"inner"` (racine $-$). Choisir la lèvre utile pour la coupe.

### Options d’affichage utiles (dans `plot_results`)
- `alpha_blue` / `alpha_orange` : transparence des surfaces (visibilité).  
- `auto_view=True` : oriente la caméra vers l’intersection pour la voir tout de suite.  
- `clip_tube2_at_cyl1=True` : **arrête le dessin du tube incliné au niveau de l’intersection** (plus lisible).  
- `figsize=(13,9)` : taille de la fenêtre ; augmenter la hauteur si nécessaire.

### Exemple minimal
```python
import numpy as np
from tubes_intersection_and_unwrap_v2 import (
    intersect_cylinders, plot_results, export_csv_both
)

R1, R2 = 1.0, 0.8
phi = np.deg2rad(30)

res = intersect_cylinders(R1, R2, phi, n_samples=2400, branch="outer")

plot_results(res, R1, R2, phi,
             auto_view=True,
             alpha_blue=0.8, alpha_orange=0.25,
             clip_tube2_at_cyl1=True,
             figsize=(13,9))

export_csv_both(res, R1, "gabarits_developpes.csv")
```

### Fichiers générés
- `gabarit_cylindre2.csv` : $(u_2,v_2)$ + points 3D — **gabarit du tube incliné**.  
- `gabarits_developpes.csv` : colonnes $(u_1,v_1)$ **gueule de loup** et $(u_2,v_2)$ **gabarit tube 2** (rééchantillonné).

### Conseils atelier
- Les courbes sont **lignes neutres** : prévoir un **jeu** si besoin (offset en DAO).  
- Unités libres mais cohérentes (mm recommandé).  
- Pour DXF/SVG : convertir la polyline $(u,v)$ en spline si exigé par la machine.

---

*Fin du document — prêt à publier / versionner.*

— **A. Verger**
