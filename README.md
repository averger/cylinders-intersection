# Intersection de deux cylindres inclinés — théorie & mode d’emploi

*Document de référence rédigé pour expliquer la géométrie et l’usage du script Python (et de la webapp Streamlit) de génération des gabarits (développés) d’intersection de deux cylindres.*

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

## 2) Utilisation du programme

### Installation rapide (script local)
Dépendances Python : `numpy`, `matplotlib` (pour le script) ou `streamlit`, `pandas`, `plotly` (pour la webapp).

```bash
# script local
pip install numpy matplotlib

# webapp
pip install streamlit numpy pandas plotly
```

### Script (CLI)
Fichier : `tubes_intersection_and_unwrap_v2.py`. Points d’entrée principaux :

- `intersect_cylinders(R1, R2, phi, n_samples=2000, branch="outer")`  
- `plot_results(res, R1, R2, phi, **options)`  
- `export_csv_cyl2(res, "gabarit_cylindre2.csv")`  
- `export_csv_both(res, R1, "gabarits_developpes.csv")`

**Exemple minimal**

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

### Webapp (facultatif)

```bash
streamlit run app.py
```

La page affiche la 3D (Plotly) et les deux développés (charts Streamlit), avec boutons d’export CSV.

---

## 3) Notes de validité et cas limites

- Pour $\phi=0$, l’équation quadratique **dégénère** ($a=0$) : on retombe sur deux cylindres coaxiaux. L’intersection est un cercle seulement si $R_1=R_2$.  
- Les formules de $a$, $b(\theta)$, $c_0(\theta)$ ci-dessus proviennent de l’égalité $x^2+y^2=R_1^2$ avec
  $x=R_2\cos\theta$, $y=R_2\sin\theta\cos\phi - t\sin\phi$. On a bien
  $$x^2+y^2=R_2^2\cos^2\theta+R_2^2\sin^2\theta\cos^2\phi-2R_2\sin\theta\cos\phi\sin\phi\,t+\sin^2\phi\,t^2,$$
  d’où le quadratique annoncé.

---

*Fin du document — prêt à publier / versionner.*

— **A. Verger**
