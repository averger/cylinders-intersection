# Post LinkedIn — lancement Cylix

> Prêt à coller. Visuels conseillés : `docs/screenshots/01-cylix-3d.png` (3D clair),
> `03-cylix-2d.png` (les deux mises à plat), `04-cylix-3d-dark.png` (sombre) —
> LinkedIn met bien en valeur un carrousel de 3–4 images + le PDF du papier en document.

---

Découper deux tubes qui se croisent, ça a l'air simple. Jusqu'au jour où il faut le faire. ✂️

En chaudronnerie, ce raccord s'appelle une « gueule de loup » : le tube incliné doit être
coupé selon une courbe gauche, et le tube principal percé d'une lumière exacte. Sur le
terrain, on trace encore souvent ça aux abaques, ou à l'œil.

J'ai construit **Cylix** pour faire les choses proprement — et j'en ai profité pour écrire
le papier qui va avec.

𝗟'𝗼𝘂𝘁𝗶𝗹 🛠️
→ Vue 3D de l'assemblage : chaque cylindre s'arrête exactement à l'intersection
→ Les deux développés calculés en formules fermées — zéro approximation
→ Gabarits annotables, export **PDF vectoriel à l'échelle 1:1** tuilé sur A4/A3/A2 avec
repères de collage : on imprime, on enroule sur le tube, on coupe
→ Export **DXF** pour la découpe laser/plasma, règle de contrôle 100 mm sur chaque planche
→ Moteur **Rust**, interface **Svelte**, études en onglets, thème clair/sombre

𝗟𝗲 𝗽𝗮𝗽𝗶𝗲𝗿 📐
La géométrie est élémentaire mais pleine de pièges (choix de branche, topologie selon le
rapport des rayons, points de tangence…). Le discriminant de l'intersection tient en une
ligne : Δ = 4·sin²φ·(R₁² − R₂²·cos²θ).

Chaque équation du papier est **vérifiée mécaniquement** — 25 contrôles symboliques et
numériques (SymPy/NumPy) — et le moteur Rust coïncide avec une implémentation NumPy
indépendante à moins de 2×10⁻¹³ mm. Autrement dit : l'epsilon machine.

Mon préféré ? Quand les deux tubes ont le même diamètre, le gabarit ne présente pas une
tangente verticale mais un point anguleux de pentes −tan(φ/2) et +cot(φ/2). Le genre de
détail qui fait la différence entre « ça a l'air juste » et « c'est juste ».

La suite : intersection cylindre / plan orienté quelconque, et export G-code direct.

Curieux d'avoir les retours des gens de la tuyauterie, de la CAO et du calcul — le papier
est en commentaire. 👇

#chaudronnerie #tuyauterie #CAO #CNC #Rust #géométrie #fabrication #openSource #svelte

---

**Commentaire d'accompagnement (à poster juste après) :**

Le papier complet (4 pages, PDF) : dérivation, cas limites, bornes d'erreur et validation
croisée — [lien vers le PDF ou le dépôt]. Les scripts de vérification sont versionnés avec
le code : chaque affirmation est rejouable en une commande.
