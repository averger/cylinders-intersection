# Théorie — Intersection de cylindres et développés de découpe

> Fondements mathématiques du moteur de calcul. Toutes les formules ci-dessous sont
> exactes (aucune approximation) ; la seule discrétisation intervient au moment de
> l'échantillonnage des courbes en polylignes, dont l'erreur est bornée en §8.
> Les longueurs sont en millimètres, les angles en radians sauf mention contraire.

---

## 1. Conventions et paramétrisations

### 1.1 Cylindre principal (cylindre 1)

Le cylindre principal, de rayon $R_1$, a son axe porté par $Oz$ :

$$
\mathcal{C}_1 : \quad x^2 + y^2 = R_1^2 .
$$

Ses coordonnées cylindriques naturelles sont $(\alpha, z)$ avec
$x = R_1\cos\alpha$, $y = R_1\sin\alpha$.

### 1.2 Cylindre secondaire (cylindre 2, « la branche »)

Le cylindre secondaire, de rayon $R_2$, est obtenu en faisant tourner d'un angle
$\varphi$ autour de l'axe $Ox$ un cylindre initialement coaxial à $Oz$. Avant
rotation, sa paramétrisation est

$$
P_0(\theta, t) = \bigl(R_2\cos\theta,\; R_2\sin\theta,\; t\bigr),
\qquad \theta \in [0, 2\pi), \; t \in \mathbb{R},
$$

où $\theta$ est l'angle autour de son axe et $t$ la coordonnée axiale locale.
La rotation autour de $Ox$ s'écrit

$$
R_x(\varphi) =
\begin{pmatrix}
1 & 0 & 0\\
0 & \cos\varphi & -\sin\varphi\\
0 & \sin\varphi & \phantom{-}\cos\varphi
\end{pmatrix},
$$

d'où la paramétrisation du cylindre incliné :

$$
\boxed{\;
\begin{aligned}
x &= R_2\cos\theta\\
y &= R_2\sin\theta\,\cos\varphi - t\,\sin\varphi\\
z &= R_2\sin\theta\,\sin\varphi + t\,\cos\varphi
\end{aligned}\;}
\tag{1}
$$

L'axe du cylindre 2 est dirigé par $R_x(\varphi)\,e_z = (0, -\sin\varphi, \cos\varphi)$ :
$\varphi = 0$ correspond à deux cylindres coaxiaux, $\varphi = \pi/2$ à deux axes
perpendiculaires. Les deux axes sont concourants à l'origine (piquage centré).

Remarque essentielle pour la suite : la rotation ayant lieu autour de $Ox$, la
coordonnée $x = R_2\cos\theta$ **ne dépend pas de $t$**. La génératrice de
$\mathcal{C}_2$ repérée par $\theta$ est donc une droite entièrement contenue
dans le plan $x = R_2\cos\theta$.

---

## 2. Courbe d'intersection : une équation du second degré exacte

Un point de $\mathcal{C}_2$ appartient à $\mathcal{C}_1$ si et seulement si
$x^2 + y^2 = R_1^2$. En substituant (1) :

$$
R_2^2\cos^2\theta + \bigl(R_2\sin\theta\cos\varphi - t\sin\varphi\bigr)^2 = R_1^2 .
$$

Développons en $t$ :

$$
\underbrace{\sin^2\varphi}_{a}\; t^2
\;\underbrace{-\,2R_2\sin\theta\cos\varphi\sin\varphi}_{b}\; t
\;+\;\underbrace{R_2^2\bigl(\cos^2\theta + \sin^2\theta\cos^2\varphi\bigr) - R_1^2}_{c}
\;=\;0 .
\tag{2}
$$

Pour $\varphi \neq 0 \pmod \pi$ (cylindres non coaxiaux), c'est une équation du
second degré en $t$ à $\theta$ fixé. Son discriminant se simplifie remarquablement :

$$
\Delta(\theta) = b^2 - 4ac
= 4\sin^2\varphi\,\Bigl(R_1^2 - R_2^2\cos^2\theta\Bigr).
\tag{3}
$$

*Vérification :* $b^2 = 4R_2^2\sin^2\theta\cos^2\varphi\sin^2\varphi$ et
$4ac = 4\sin^2\varphi\bigl(R_2^2\cos^2\theta + R_2^2\sin^2\theta\cos^2\varphi - R_1^2\bigr)$ ;
les termes en $\sin^2\theta\cos^2\varphi$ s'annulent, il reste (3). ∎

Les solutions sont donc, pour $\varphi \in (0, \pi)$ (donc $\sin\varphi > 0$,
convention adoptée dans tout ce document — pour $\sin\varphi < 0$ les deux
branches s'échangent) :

$$
\boxed{\;
t_\pm(\theta) \;=\;
\frac{R_2\sin\theta\cos\varphi \;\pm\; \sqrt{R_1^2 - R_2^2\cos^2\theta}}{\sin\varphi}
\;}
\tag{4}
$$

### 2.1 Interprétation géométrique du discriminant

La condition d'existence $\Delta \geq 0$ s'écrit

$$
R_2\,\lvert\cos\theta\rvert \;\leq\; R_1 .
\tag{5}
$$

C'est exactement la condition pour que la génératrice de $\mathcal{C}_2$ à
l'angle $\theta$ — droite du plan $x = R_2\cos\theta$ — rencontre le cylindre
$\mathcal{C}_1$, qui n'occupe que la tranche $\lvert x\rvert \leq R_1$. Les deux
racines $t_-$ et $t_+$ sont les points d'**entrée** et de **sortie** de la
génératrice à travers la paroi du cylindre principal.

### 2.2 Topologie de l'intersection

* **$R_2 \leq R_1$ (piquage classique)** : (5) est vraie pour tout $\theta$ ;
  chaque génératrice traverse. L'intersection est formée de **deux courbes
  fermées** disjointes (entrée $t_-$ et sortie $t_+$), chacune paramétrée par
  $\theta \in [0,2\pi)$. Si $R_2 < R_1$, $\Delta > 0$ strictement et chaque
  branche est une courbe analytique lisse.
* **$R_2 = R_1$** : $\Delta$ s'annule en $\theta = 0$ et $\theta = \pi$, où les
  deux branches se touchent : les courbes d'entrée et de sortie se croisent en
  deux points (les points de tangence des deux cylindres). C'est le cas du
  solide de Steinmetz (§6.1).
* **$R_2 > R_1$** : seuls les $\theta$ tels que $\lvert\cos\theta\rvert \leq R_1/R_2$
  donnent une intersection ; la branche « perce » de part en part et
  l'intersection forme deux courbes fermées faisant chacune le tour du
  cylindre **principal** (et non plus de la branche).

Le paramètre `branch` du moteur sélectionne la racine : `outer` $= t_+$,
`inner` $= t_-$.

*Remarque (branches vs composantes connexes).* Pour $R_2 \leq R_1$, chaque
branche $t_\pm$ coïncide avec une composante connexe de l'intersection. Pour
$R_2 > R_1$ en revanche, les deux composantes connexes (qui encerclent chacune
le cylindre principal) mélangent points d'entrée et de sortie : les lieux
$\{t_+\}$ et $\{t_-\}$ restent des courbes bien définies — sortie et entrée de
paroi — mais ne sont plus des composantes connexes de l'intersection. C'est le
lieu de sortie (resp. d'entrée) qui est pertinent pour le traçage.

---

## 3. Développement d'un cylindre : une isométrie exacte

Le cylindre est une surface **développable** : sa courbure de Gauss est nulle,
il s'applique donc sur le plan sans distorsion. Concrètement, la première forme
fondamentale du cylindre de rayon $R$ paramétré par $(\theta, t)$ est

$$
\mathrm{d}s^2 = R^2\,\mathrm{d}\theta^2 + \mathrm{d}t^2 ,
$$

et l'application de déroulage

$$
(\theta, t) \;\longmapsto\; (u, v) = (R\,\theta,\; t)
\tag{6}
$$

vérifie $\mathrm{d}u^2 + \mathrm{d}v^2 = \mathrm{d}s^2$ : c'est une **isométrie**.
Toute courbe tracée sur le cylindre conserve sa longueur, et les angles avec les
génératrices sont conservés. Un gabarit tracé dans le plan $(u,v)$ puis enroulé
sur le tube reproduit donc *exactement* la courbe 3D — c'est le fondement de
tous les tracés de chaudronnerie à l'échelle 1:1.

### 3.1 Développé sur le cylindre 2 — gabarit du tube coupé

En appliquant (6) au cylindre 2 avec la solution (4) :

$$
\boxed{\;
u_2 = R_2\,\theta, \qquad
v_2 = t_\pm(\theta)
\;}
\tag{7}
$$

La largeur totale du gabarit est le périmètre $2\pi R_2$ ; la courbe $v_2(u_2)$
est la ligne de coupe du tube incliné, à reporter sur la tôle roulée ou sur le
tube via une feuille enroulée.

**Symétries.** $t_\pm$ vérifie $t_\pm(\pi - \theta) = t_\pm(\theta)$ au signe du
terme en $\sin\theta$ près ; plus utile en pratique :
$t_+(\theta) + t_-(\theta) = 2R_2\sin\theta\cos\varphi/\sin\varphi$ (somme des
racines) et $t_+(-\theta) = -t_-(\theta)$. Le gabarit `outer` se déduit du
gabarit `inner` par la symétrie $(\theta, t) \mapsto (-\theta, -t)$, image du
retournement du tube.

### 3.2 Développé sur le cylindre 1 — la « gueule de loup »

La même courbe d'intersection, vue depuis le cylindre principal, définit la
découpe de la **lumière** (le trou) : la *gueule de loup*. On repère un point de
la courbe par les coordonnées cylindriques de $\mathcal{C}_1$ :

$$
\alpha = \operatorname{atan2}(y, x), \qquad v_1 = z,
\qquad u_1 = R_1\,\alpha .
\tag{8}
$$

Comme le point est sur $\mathcal{C}_1$, $x = R_1\cos\alpha$ ; or $x = R_2\cos\theta$
par (1), d'où la relation de transfert exacte entre les deux angles :

$$
\cos\alpha = \frac{R_2}{R_1}\,\cos\theta .
\tag{9}
$$

Pour $R_2 \leq R_1$, $\alpha$ reste dans une bande $[\,\alpha_{\min},\,\alpha_{\max}]$
autour de $\pm\pi/2$ : la gueule de loup n'occupe qu'une fraction du périmètre
du gros tube, comme attendu. La suite $\alpha(\theta)$ est continue mais
$\operatorname{atan2}$ la replie dans $(-\pi,\pi]$ : le moteur applique un
*unwrap* (suppression des sauts de $2\pi$) avant de poser $u_1 = R_1\alpha$, et
**conserve l'ordre de parcours en $\theta$** — les échantillons suivent ainsi le
contour physique de la lumière (un tri par $u_1$ entrelacerait les deux lèvres).
La topologie se lit sur l'excursion totale de $\alpha$ : si
$\lvert\alpha(2\pi) - \alpha(0)\rvert < \pi$, la courbe revient sur elle-même et
la lumière est un **contour fermé** (pénétration complète, $R_2 \leq R_1$) ;
si $\alpha$ a progressé de $2\pi$, la courbe **fait le tour** du cylindre
principal ($R_2 > R_1$) et le développé est une courbe ouverte sur toute la
largeur $2\pi R_1$.

### 3.3 Cas perpendiculaire : formule fermée classique

Pour $\varphi = \pi/2$ et la branche traversante, (4) donne
$t_\pm = \pm\sqrt{R_1^2 - R_2^2\cos^2\theta}$ et $z = R_2 \sin\theta$. En combinant
avec (9) ($\cos\theta = \tfrac{R_1}{R_2}\cos\alpha$) :

$$
v_1(\alpha) \;=\; z \;=\; \pm\sqrt{R_2^2 - R_1^2\cos^2\alpha} ,
\tag{10}
$$

définie pour $\lvert\cos\alpha\rvert \leq R_2/R_1$ : c'est la formule
traditionnelle de traçage de la gueule de loup perpendiculaire, retrouvée ici
comme cas particulier du modèle général. Le moteur n'utilise pas (10) : il passe
toujours par la forme paramétrique exacte (1)–(4)–(8), valable pour tout
$\varphi$, dont (10) est un contrôle de non-régression.

---

## 4. Coupe d'un cylindre par un plan (« sifflet »)

Le second mode traite un tube unique coupé par un plan incliné. Le plan est
défini par son inclinaison $\varphi$ autour de $Ox$ et son décalage $z_0$ :

$$
\mathcal{P} : \quad z = z_0 - y\,\tan\varphi ,
\qquad \varphi \in \left(-\tfrac{\pi}{2}, \tfrac{\pi}{2}\right),
\tag{11}
$$

de normale $n = (0, \sin\varphi, \cos\varphi)$ passant par $(0, 0, z_0)$.

### 4.1 Courbe d'intersection

En injectant la paramétrisation de $\mathcal{C}_1$ :

$$
\gamma(\theta) = \bigl(R_1\cos\theta,\; R_1\sin\theta,\; z_0 - R_1\tan\varphi\,\sin\theta\bigr).
\tag{12}
$$

C'est une **ellipse** : la section d'un cylindre par un plan non parallèle à son
axe est toujours elliptique, de demi-petit axe $R_1$ (direction $Ox$, contenue
dans le plan de coupe) et de demi-grand axe $R_1/\cos\varphi$ (l'axe du tube
étant vu sous l'incidence $\varphi$).

### 4.2 Développé : une sinusoïde pure

Avec (6) :

$$
\boxed{\;
v(u) \;=\; z_0 \;-\; R_1\tan\varphi\,\sin\!\frac{u}{R_1}
\;}
\tag{13}
$$

Le développé d'une coupe plane est une **sinusoïde exacte** de période
$2\pi R_1$ (le périmètre) et d'amplitude crête à crête $2R_1\tan\varphi$. Ce
résultat classique fournit deux contrôles immédiats : à $\varphi = 0$ le
développé est une droite (coupe droite), et l'amplitude diverge quand
$\varphi \to \pi/2$ (le plan devient parallèle à l'axe — cas exclu par (11)).

### 4.3 Plan pleinement orienté (deux inclinaisons)

Le cas général d'un plan orienté par deux basculements — $\varphi_x$ autour de
$Ox$ et $\varphi_y$ autour de $Oy$ — s'écrit
$z = z_0 - y\tan\varphi_x - x\tan\varphi_y$, de normale proportionnelle à
$(\tan\varphi_y,\ \tan\varphi_x,\ 1)$. En posant $a = \tan\varphi_y$ et
$b = \tan\varphi_x$, le développé reste une **sinusoïde pure, déphasée** :

$$
v(u) \;=\; z_0 \;-\; R_1\sqrt{a^2+b^2}\,
\sin\!\Bigl(\frac{u}{R_1} + \psi\Bigr),
\qquad \psi = \operatorname{atan2}(a, b),
\tag{13'}
$$

car $b\sin\theta + a\cos\theta = \sqrt{a^2+b^2}\,\sin(\theta+\psi)$.
L'amplitude crête à crête vaut $2R_1\sqrt{\tan^2\varphi_x+\tan^2\varphi_y}$
et l'angle réel entre le plan et la section droite est
$\arctan\sqrt{a^2+b^2}$ ; le second angle ne fait que **tourner la phase** du
gabarit autour du tube. (13) est le cas $a = 0$.

### 4.4 Nœud multi-tubes (N piquages sur un tube principal)

Un tube principal peut recevoir $N$ piquages, chacun décrit par son rayon
$r_i \le R_1$, sa hauteur $z_i$ sur l'axe principal, son inclinaison
$\varphi_i \in (0, \pi)$ et son **azimut** $\psi_i$ autour du tube. L'axe du
piquage $i$ passe par $C_i = (0, 0, z_i)$ avec la direction
$d_i = R_z(\psi_i)\,R_x(\varphi_i)\,e_z$, et sa génératrice d'angle $\theta$
est la droite $P(t) = C_i + r_i(\cos\theta\,u_i + \sin\theta\,w_i) + t\,d_i$.
Par défaut $z_i = 0$ : **tous les axes concourent au centre du tube
principal**, comme dans les autres modes, soit le point de travail d'un nœud
de treillis sans excentrement. $z_i \neq 0$ décale le point de percée d'un
piquage le long de la membrure et introduit l'**excentrement** du §4.5.

Chaque piquage **s'arrête au nœud** : aucun tube ne traverse un autre.
Comme deux cylindres qui se croisent ne peuvent pas *tous deux* s'arrêter
l'un sur l'autre en refermant le joint, il faut une **convention de
priorité** : l'ordre de la liste — le piquage 1 est le **membre
traversant**, les suivants viennent se poser sur lui (l'assemblage à
recouvrement des nœuds K/KT en charpente tubulaire).

Un tube réel arrive de l'extérieur : chaque génératrice ne garde que son
**segment libre le plus extérieur**, et sa lèvre basse est le *dernier*
obstacle rencontré en entrant,

$$
t_{\text{lèvre}}(\theta)
\;=\; \max\Bigl(t^{+}_{\text{principal}}(\theta),\
\max_{j<i}\ t^{\text{sortie}}_{j}(\theta)\Bigr),
$$

où $t^{+}_{\text{principal}}$ est l'atterrissage sur le tube principal
(plus grande racine du trinôme droite–cylindre) et
$t^{\text{sortie}}_{j}$ la **grande racine** du trinôme du prioritaire
$j$ — le point où la génératrice ressort de son cylindre, c'est-à-dire son
**dos**. Un prioritaire ne compte que dans sa **longueur physique** (entre
son plancher et son plafond axiaux — l'extension infinie du cylindre
n'arrête rien).

La lèvre est **une seule courbe continue** : sinusoïde d'atterrissage là où
le tube atteint la paroi du principal, **selle** là où il chevauche le dos
d'un prioritaire — les deux se raccordent exactement au bord de l'empreinte
commune, où $t^{\text{sortie}}_{j} = t^{+}_{\text{principal}}$. Chaque
point levé est, par construction, à distance $r_j$ de l'axe du prioritaire :
la couture épouse une paroi qui existe, le joint est fermé, et le pourtour
complet du tube est conservé (un tube ne perd jamais de génératrices, il
les pose ailleurs).

Sur le tube principal, chaque piquage découpe sa gueule de loup isolée
(§3.2) **translatée** de $(R_1\psi_i,\ z_i)$ dans le plan déroulé : la
rotation d'azimut se déroule en décalage horizontal, la translation le long
de l'axe en décalage vertical — les deux exactement. Lorsque deux lumières
se chevauchent, la découpe résultante est l'union des deux contours ; le
moteur le signale.

Ces énoncés sont vérifiés par les tests du moteur : un nœud à un seul
piquage reproduit (4) à $10^{-9}$ près, la lèvre du membre traversant n'est
jamais altérée, et — propriété de fermeture — chaque point levé de la lèvre
d'un chevauchant est à distance $r_j$ de l'axe d'un prioritaire (à
$10^{-9}$ près) tout en restant hors du tube principal : la couture repose
exactement sur la paroi qui l'arrête (aucun jour dans le joint).

### 4.5 Excentrement, jeu et recouvrement du nœud

Deux piquages de même azimut (ou opposés de $\pi$) sont **coplanaires** avec
l'axe de la membrure : leurs axes se coupent réellement, en un point dont la
distance à l'axe de la membrure est l'**excentrement** $e$ du nœud. Dans le
plan du nœud, l'axe du piquage $i$ s'écrit
$(a, b) = (0, z_i) + t(\sin\varphi_i, \cos\varphi_i)$, et l'élimination
donne directement

$$
e \;=\; (z_j - z_i)\,
\frac{\sin\varphi_i\,\sin\varphi_j}{\sin(\varphi_j - \varphi_i)},
\tag{16}
$$

le dénominateur devenant $\sin(\varphi_i + \varphi_j)$ pour deux piquages
opposés. Axes concourants ($z_i = z_j$) donne $e = 0$ ; axes parallèles
($\varphi_i = \varphi_j$) annule le dénominateur, les axes ne se coupant
jamais. Le facteur
$\sin\varphi_i\sin\varphi_j / \sin(\varphi_i + \varphi_j)$ est celui de
la formule d'excentrement des nœuds K de l'**EN 1993-1-8** : en prenant
$\varphi_i = \theta_1$ pour la barre montante et
$\varphi_j = \pi - \theta_2$ pour la descendante,
$\sin(\varphi_j - \varphi_i) = \sin(\theta_1 + \theta_2)$.

Le **jeu** $g$ et le **taux de recouvrement** $\lambda_{ov}$ se lisent sur la
génératrice du plan du nœud, d'abscisse $u = R_1(\psi - \pi/2)$ dans le
développé — le sommet de la lumière, puisque le piquage penche selon
$n(\psi) = (\sin\psi, -\cos\psi, 0)$, d'azimut $\psi - \pi/2$. Chaque
empreinte y découpe un segment $[v_{\min}, v_{\max}]$ ; si les deux segments
sont disjoints, leur écart est le jeu, mesuré le long de la surface (le
développé étant une isométrie, cette distance est la vraie distance sur le
tube). S'ils se chevauchent, $\lambda_{ov} = q/p$ où $q$ est le
recouvrement et $p$ l'étendue de l'empreinte du piquage **chevauchant** prise
seule, conformément à la définition de la norme.

Une conséquence utile en atelier : un piquage perce la paroi là où son propre
axe la traverse, à la distance $R_1/\sin\varphi$ du centre le long de son
axe. Deux piquages **concourants** de sens opposés laissent donc
nécessairement un jeu entre leurs lumières ; le rapprocher d'un recouvrement
donné passe obligatoirement par un excentrement, ce que (16) chiffre
directement.

Le moteur signale les deux dépassements normatifs : $|e| > 0{,}25\,\varnothing_1$
et $\lambda_{ov} < 25\,\%$.

---

## 5. Choix de la branche et points de tangence

Aux angles $\theta$ où $\Delta(\theta) = 0$, soit $R_2\lvert\cos\theta\rvert = R_1$,
les deux branches $t_+$ et $t_-$ coïncident : la génératrice est **tangente** au
cylindre principal. Le comportement local dépend de la **multiplicité** du zéro
de $\Delta$ :

* **$R_2 > R_1$ — zéro simple** en $\cos\theta^* = \pm R_1/R_2$ : le radicande
  s'annule linéairement, $\sqrt{\smash[b]{\cdot}}$ varie en
  $\sqrt{\lvert\theta - \theta^*\rvert}$ et $t_\pm'(\theta)$ **diverge** — le
  développé du tube présente une véritable tangente verticale dans le plan
  $(u_2, v_2)$. La gueule de loup atteint alors son extrémité en $u_1$ : par
  (9), $\lvert\cos\alpha\rvert = 1$ précisément quand
  $R_2\lvert\cos\theta\rvert = R_1$.
* **$R_2 = R_1 = R$ — zéro double** en $\theta^* \in \{0, \pi\}$ : le radicande
  vaut $R^2\sin^2\theta$, sa racine $R\lvert\sin\theta\rvert$ est **lipschitzienne**
  et les pentes unilatérales du développé restent finies. Un calcul direct au
  voisinage de $\theta^* = 0$ donne, pour la branche $t_+$,

$$
\left.\frac{\mathrm{d}v_2}{\mathrm{d}u_2}\right|_{0^-} = -\tan\frac{\varphi}{2},
\qquad
\left.\frac{\mathrm{d}v_2}{\mathrm{d}u_2}\right|_{0^+} = +\cot\frac{\varphi}{2} :
$$

  le gabarit présente un **point anguleux** (et non une tangente verticale) —
  à $\varphi = \pi/2$ on retrouve les pentes $\mp 1$ des arches en
  $\lvert\sin\rvert$ du cas Steinmetz (§6.1).

Pour $R_2 < R_1$, $\Delta > 0$ partout : aucune tangence, courbes analytiques
lisses. Le moteur échantillonne $\theta$ uniformément sur $[0, 2\pi)$, évalue
(4), et écarte les échantillons à discriminant négatif ; aucun traitement
spécial n'est nécessaire aux points de tangence puisque (4) y est continue.

---

## 6. Cas limites et dégénérescences

### 6.1 Solide de Steinmetz ($R_1 = R_2$, $\varphi = \pi/2$)

Pour deux cylindres égaux et perpendiculaires, (2) se factorise : l'intersection
est la réunion de **deux ellipses planes** contenues dans les plans $z = y$ et
$z = -y$ (avec nos conventions). Chaque branche $t_\pm(\theta) = \pm R\lvert\sin\theta\rvert$
suit alors une ellipse sur chaque demi-période et passe de l'une à l'autre aux
points de tangence $\theta \in \{0, \pi\}$ : le développé `outer` présente les
arches en $\lvert\sin\rvert$ caractéristiques, avec points anguleux — c'est un
comportement géométrique réel, pas un artefact numérique.

### 6.2 Cylindres coaxiaux ($\varphi = 0$)

L'équation (2) perd son terme quadratique ($a = 0$) et $b = 0$ : elle n'a pas de
solution isolée ($c \neq 0$) ou en a une infinité ($R_1 = R_2$, cylindres
confondus). Le cas est rejeté par validation d'entrée ($\varphi$ borné loin
de $0$) et le moteur retourne proprement une courbe vide si tous les
échantillons sont invalides.

### 6.3 Plan quasi parallèle à l'axe ($\lvert\varphi\rvert \to \pi/2$, mode plan)

$\tan\varphi$ diverge ; l'API borne $\lvert\varphi\rvert < \pi/2 - \varepsilon$
et renvoie une erreur explicite au-delà.

---

## 7. Grandeurs dérivées utiles à la fabrication

Pour un développé $v(u)$ échantillonné en $(u_k, v_k)$ :

* **Largeur de tôle** : $2\pi R$ (périmètre exact du tube développé) ;
* **Hauteur de gabarit** : $\max_k v_k - \min_k v_k$ ; pour le sifflet,
  exactement $2R_1\lvert\tan\varphi\rvert$ ;
* **Longueur de la ligne de coupe** :
  $L = \sum_k \sqrt{(u_{k+1}-u_k)^2 + (v_{k+1}-v_k)^2}$, qui converge vers la
  longueur vraie de la courbe 3D (l'isométrie (6) préserve les longueurs) —
  utile pour estimer le temps de découpe et la longueur de cordon de soudure ;
* **Repère $\theta = 0$** : la génératrice origine est tracée sur les gabarits
  pour caler la feuille sur le tube (alignement du zéro angulaire).

---

## 8. Discrétisation et bornes d'erreur

Les formules (4), (7), (12), (13) sont exactes. La seule approximation est le
remplacement de la courbe par la polyligne de ses échantillons. Pour un pas
angulaire $h = 2\pi/N$, l'écart maximal entre une corde et l'arc qu'elle
sous-tend est borné par la flèche

$$
e \;\leq\; \frac{\kappa_{\max}}{8}\,\ell^2 ,
$$

où $\ell$ est la longueur de corde et $\kappa_{\max}$ la courbure maximale de la
courbe développée sur le segment. Avec le réglage par défaut $N = 1440$
($h = 0{,}25^\circ$) et des tubes de diamètre $\leq 500$ mm, l'écart reste très
inférieur à $0{,}01$ mm partout où la courbure est bornée — sans objet pratique
devant le trait de coupe. Près des points de tangence (§5) la courbure diverge
mais l'erreur reste contrôlée par l'échantillonnage dense et la continuité
de (4) ; augmenter $N$ resserre localement la polyligne.

Le *unwrap* angulaire (§3.2) est numériquement robuste : les échantillons étant
denses, les sauts de $2\pi$ de $\operatorname{atan2}$ sont détectés par simple
comparaison au demi-tour, exactement comme `numpy.unwrap`.

---

## 9. Validation

Chaque affirmation de ce document est vérifiée mécaniquement — les scripts sont
versionnés dans [`docs/verification/`](verification/) :

* **Vérification symbolique et numérique**
  (`verify_theory.py`, SymPy + NumPy, 25 contrôles) : coefficients (2),
  discriminant (3), racines (4) sur 200 tirages aléatoires, condition
  d'existence (5), première forme fondamentale (6), symétries §3.1, relation de
  transfert (9), formule fermée (10), appartenance au plan (12), demi-axes de
  l'ellipse §4.1, sinusoïde (13) et son amplitude, factorisation de Steinmetz
  §6.1, divergence de pente au zéro simple et pentes unilatérales
  $-\tan(\varphi/2)$ / $+\cot(\varphi/2)$ au zéro double (§5) ;
* **Validation croisée du moteur** (`verify_engine.py`) : le moteur Rust et
  l'implémentation NumPy de référence (`tubes_intersection_and_unwrap.py`)
  coïncident à **moins de $2\times10^{-13}$ mm** (courbe 3D, développé branche,
  gueule de loup) sur 6 configurations couvrant $R_2 < R_1$, $R_2 = R_1$,
  $R_2 > R_1$, les deux branches et $\varphi \in \{20°, 45°, 60°, 75°, 90°\}$ ;
* **Tests unitaires Rust** (`cargo test`) : amplitude du cas perpendiculaire,
  période complète du sifflet, continuité du contour de la gueule de loup
  (ordre de parcours, détection de fermeture), plans de tuilage et writers
  PDF/DXF ;
* **Contrôle métrologique à l'impression** : chaque export porte une règle de
  référence de 100 mm à vérifier au réglet après impression « taille réelle ».

---

## Références

1. Géométrie descriptive classique : intersections de surfaces de révolution et
   développements (traités de chaudronnerie et de traçage industriel).
2. Solide de Steinmetz — intersection de deux cylindres égaux perpendiculaires.
3. do Carmo, *Differential Geometry of Curves and Surfaces* — surfaces
   développables, première forme fondamentale, isométries.
