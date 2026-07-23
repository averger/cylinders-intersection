# Post LinkedIn — lancement Cylix

> Prêt à coller. Visuels conseillés : `docs/screenshots/01-cylix-3d.png` (3D clair),
> `03-cylix-2d.png` (les deux mises à plat), `04-cylix-3d-dark.png` (sombre) —
> carrousel de 3–4 images + le PDF du papier en document.

---

Cylix calcule les développés d'intersections de cylindres : gueule de loup, piquage
incliné, coupe en sifflet. Formules fermées, sans approximation.

Deux gabarits par intersection : le contour du tube sécant et la lumière à percer dans
le tube principal. Export PDF à l'échelle 1:1, à enrouler sur le tube et couper
directement. Export DXF par calques pour la découpe plasma. Repères d'alignement et
génératrice de référence tracés sur chaque gabarit.

Mode plan incliné à deux basculements pour les coupes d'onglet composées.

La théorie est dérivée intégralement et vérifiée mécaniquement — 26 contrôles
symboliques et numériques, publiés avec le code.

Moteur Rust, interface web.

#chaudronnerie #tuyauterie #métallerie #ingenierie

---

**Commentaire d'accompagnement (à poster juste après) :**

Le papier complet (PDF) : dérivation, cas limites, bornes d'erreur et validation
croisée — [lien vers le PDF ou le dépôt].
