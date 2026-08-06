"""Monte la démo à 15 s : une vitesse par phase, pour garder un rythme lisible.

Lit marks.json (écrit par demo-video.mjs) et la table de durées cibles
ci-dessous.  Deux règles à ne pas oublier :
  * -ss / -to AVANT -i (après, ce sont des options de sortie appliquées
    après le retiming) ;
  * fps=30 AVANT setpts (normaliser en cadence constante d'abord, sinon la
    durée obtenue est imprévisible sur une source à cadence irrégulière).
"""
import json
import os
import subprocess
import sys

FF = open(".ffmpeg").read().strip()
RAW = os.path.join("video-raw", [f for f in os.listdir("video-raw") if f.endswith(".webm")][0])
M = json.load(open("marks.json"))
OUT = sys.argv[1] if len(sys.argv) > 1 else "cylix-demo-cylcyl.mp4"

# phase -> durée finale (s).  La somme fait la durée de la vidéo.
TARGET = {
    "avant": 2.72,      # les gabarits avant modification (plan tenu)
    "vers3d": 0.6,     # bascule vers la 3D
    "orbite": 2.4,     # on tourne la pièce
    "diametre": 3.2,   # Ø₂ : 100 -> 80 mm
    "angle": 2.8,      # 45° -> 30°
    "clic2d": 0.6,     # retour en 2D
    "apres": 2.52,      # les gabarits après (plan tenu)
}
ORDER = ["avant", "vers3d", "orbite", "diametre", "angle", "clic2d", "apres"]
NEXT = {"avant": "vers3d", "vers3d": "orbite", "orbite": "diametre", "diametre": "angle",
        "angle": "clic2d", "clic2d": "apres", "apres": "fin"}


def frames(path: str) -> int:
    err = subprocess.run([FF, "-i", path, "-f", "null", "-"], capture_output=True, text=True).stderr
    for line in err.splitlines()[::-1]:
        if "frame=" in line:
            return int(line.split("frame=")[1].split()[0])
    return 0


os.makedirs("seg", exist_ok=True)
parts = []
for name in ORDER:
    a, b = M[name], M[NEXT[name]]
    k = (b - a) / TARGET[name]
    out = f"seg/{name}.mp4"
    subprocess.run(
        [FF, "-v", "error", "-y", "-ss", f"{a:.3f}", "-to", f"{b:.3f}", "-i", RAW,
         "-vf", f"fps=30,setpts=PTS/{k:.6f},fps=30,format=yuv420p", "-an",
         "-c:v", "libx264", "-preset", "slow", "-crf", "19", out],
        check=True,
    )
    n = frames(out)
    parts.append(out)
    print(f"{name:9s} source {b - a:6.2f}s ×{k:5.2f} -> {n / 30:5.2f}s ({n} images)")

with open("seg/list.txt", "w") as f:
    for p in parts:
        f.write(f"file '{os.path.basename(p)}'\n")
subprocess.run(
    [FF, "-v", "error", "-y", "-f", "concat", "-safe", "0", "-i", "seg/list.txt",
     "-c", "copy", "-movflags", "+faststart", OUT],
    check=True,
)
n = frames(OUT)
print(f"TOTAL {n / 30:.2f}s ({n} images) · {os.path.getsize(OUT) / 1e6:.2f} Mo -> {OUT}")
subprocess.run(
    [FF, "-v", "error", "-y", "-i", OUT,
     "-vf", "select=not(mod(n\\,50)),scale=426:-1,tile=3x3", "-frames:v", "1", "demo-contact.png"],
    check=True,
)
