# Assets de communication

## Démo vidéo 15 s — cas cylindre / cylindre, avant / après

`demo-video.mjs` enregistre la démonstration avec Playwright, `retime.py` la
monte à 15 s.

Scénario : les deux gabarits 2D **avant** modification (Ø₂ = Ø₁, 45°), passage
en 3D et rotation de la pièce, Ø₂ de 100 à 80 mm (0,8·Ø₁), angle des axes de
45° à 30°, retour en 2D pour montrer les gabarits **après**. Thème clair,
**aucun sous-titre** : la démo se lit sur les valeurs du panneau et la forme
des gabarits. Filmée en 1920×1080, la taille pour laquelle l'app est dessinée,
et réduite à 1280×720 à l'encodage — filmer directement en 720p grossit tout
d'un tiers, les éléments d'interface étant à taille fixe.

```bash
cargo build --release && ./target/release/cylix &   # serveur sur :8787
cd docs/marketing
node demo-video.mjs                                 # -> video-raw/*.webm + marks.json
python retime.py cylix-demo-cylcyl.mp4              # -> 15 s, H.264 faststart
```

`retime.py` attend le chemin d'un binaire ffmpeg dans un fichier `.ffmpeg`
(par exemple celui fourni par le paquet Python `imageio-ffmpeg`).

### Pourquoi une étape de montage

Le rendu WebGL logiciel d'une machine sans GPU tourne à quelques images par
seconde, et chaque événement d'entrée déclenche en plus un recalcul complet :
l'enregistrement brut de ces mêmes gestes dure environ 80 s. Chaque phase est
donc ré-accélérée séparément — les plans tenus sur les gabarits sont même
légèrement ralentis, les glissements de curseurs accélérés d'un facteur 9 —
ce qui donne un rythme lisible que l'uniforme ne donnerait pas. Avec un vrai
GPU l'enregistrement est déjà fluide et le montage devient inutile.

Le script masque la pastille « calcul… » pendant l'enregistrement : c'est un
état réel de l'app, mais le rendu logiciel la fait clignoter en permanence
alors qu'avec un GPU le calcul est imperceptible. La masquer rend la démo
plus représentative de l'expérience réelle, pas moins.

### Trois pièges rencontrés, à ne pas reproduire

* `-ss` / `-to` doivent être placés **avant** `-i`. Placés après, ce sont des
  options de sortie : le découpage s'applique après le retiming et la durée
  obtenue n'a plus de sens.
* Dans le filtre, normaliser en cadence constante **avant** de retimer :
  `fps=30,setpts=PTS/k,fps=30`. L'inverse, sur une source à cadence
  irrégulière, donne une durée imprévisible.
* La vidéo Playwright n'enregistre pas le pointeur système : le script dessine
  un curseur factice, suivi **en page** par un écouteur `mousemove` pour ne pas
  payer un aller-retour CDP par déplacement.
