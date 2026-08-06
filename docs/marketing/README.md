# Assets de communication

## `demo-video.mjs` — démo vidéo 15 s (cas cylindre / cylindre)

Enregistre une démonstration de l'app avec Playwright : rotation de la 3D,
Ø₂ de 100 à 80 mm (0,8·Ø₁), angle de 45° à 30°, puis bascule en 2D pour
montrer que les deux gabarits ont suivi. Thème clair, 1280×720.

```bash
cargo build --release && ./target/release/cylix &   # serveur sur :8787
cd docs/marketing && node demo-video.mjs            # -> video-raw/*.webm
```

Le script écrit dans la console les frontières de chaque phase
(`orbite`, `diamètre`, `angle`, `2D`) et l'`OFFSET_MS` du chargement.

### Montage à 15 s

Le rendu WebGL logiciel d'un environnement sans GPU tourne à quelques
images par seconde, donc l'enregistrement brut dure une minute et demie.
Chaque phase est ensuite ré-accélérée séparément pour tenir 15 s avec un
rythme lisible (l'orbite ×4, les curseurs ×8, la révélation des gabarits
au ralenti). Avec un vrai GPU, l'enregistrement est déjà fluide et
l'étape de montage devient inutile.

Deux pièges rencontrés, à ne pas reproduire :

* `-ss` / `-to` doivent être placés **avant** `-i`. Placés après, ce sont
  des options de sortie : le découpage s'applique après le retiming et la
  durée obtenue n'a plus de sens.
* dans le filtre, normaliser en cadence constante **avant** de retimer :
  `fps=30,setpts=PTS/k,fps=30`. L'inverse, sur une source à cadence
  irrégulière, donne une durée imprévisible.

```bash
ffmpeg -ss <début> -to <fin> -i video-raw/*.webm \
  -vf "fps=30,setpts=PTS/<k>,fps=30,format=yuv420p" -an \
  -c:v libx264 -preset slow -crf 19 seg/<phase>.mp4
ffmpeg -f concat -safe 0 -i seg/list.txt -c copy \
  -movflags +faststart cylix-demo-cylcyl.mp4
```
