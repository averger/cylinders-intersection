# tubes_intersection_and_unwrap_v2.py
# Intersection de deux cylindres + doubles développés (gabarits)
# - Cylindre 1 : axe Oz (vertical), rayon R1
# - Cylindre 2 : cylindre de rayon R2 incliné d’un angle phi autour de l’axe x
#   Param avant rotation : (R2*cosθ, R2*sinθ, t)
#   Rotation Rx(phi) :
#       x = R2*cosθ
#       y = R2*sinθ*cosφ - t*sinφ
#       z = R2*sinθ*sinφ + t*cosφ
#
# Sorties :
#   - Courbe d'intersection 3D
#   - Développé sur cylindre 2 (u2,v2) = (R2*θ, t)
#   - Développé sur cylindre 1 (u1,v1) = (R1*α, z)  [« gueule de loup »]
#
# Usage :
#   python tubes_intersection_and_unwrap_v2.py
#
# Dépendances : numpy, matplotlib

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


# ============================
# Structures de données
# ============================

@dataclass
class IntersectionResult:
    theta: np.ndarray        # paramètre angulaire autour du cylindre incliné
    t: np.ndarray            # coordonnée axiale locale (avant rotation) sur le cylindre incliné
    xyz: np.ndarray          # points 3D de l'intersection (N,3)
    uv: np.ndarray           # développé sur cylindre 2 : (u=R2*theta, v=t) (N,2)
    valid_mask: np.ndarray   # masque des θ valides (discriminant >= 0)


# ============================
# Calcul de l'intersection
# ============================

def intersect_cylinders(R1: float,
                        R2: float,
                        phi: float,
                        n_samples: int = 2000,
                        branch: str = "outer") -> IntersectionResult:
    """
    Calcule la courbe d'intersection de deux cylindres :
      - Cylindre 1 (vertical, axe Oz) : x^2 + y^2 = R1^2
      - Cylindre 2 (incliné d'un angle phi autour de l'axe x), rayon R2
        Paramétrisation (avant rotation autour de x) : (R2*cosθ, R2*sinθ, t)
        Après rotation Rx(phi) :
          x = R2*cosθ
          y = R2*sinθ*cosφ - t*sinφ
          z = R2*sinθ*sinφ + t*cosφ

    En imposant x^2+y^2=R1^2 on obtient un quadratique en t(θ).
    branch : "outer" | "inner" pour choisir la racine (+ ou -).
    """
    theta = np.linspace(0.0, 2.0*np.pi, n_samples, endpoint=False)

    s = np.sin(phi)
    c = np.cos(phi)

    # Quadratique : a t^2 + b t + c0 = 0
    a = s**2 * np.ones_like(theta)
    b = -2.0 * R2 * np.sin(theta) * c * s
    c0 = R2**2 * (np.cos(theta)**2 + (np.sin(theta)**2) * c**2) - R1**2

    disc = b*b - 4.0*a*c0
    valid = disc >= 0.0

    t = np.full_like(theta, np.nan)
    sqrt_disc = np.zeros_like(theta)
    sqrt_disc[valid] = np.sqrt(disc[valid])

    t_plus  = (-b + sqrt_disc) / (2.0*a + 1e-300)
    t_minus = (-b - sqrt_disc) / (2.0*a + 1e-300)

    if branch.lower() in ("outer", "max", "+"):
        t_sel = t_plus
    elif branch.lower() in ("inner", "min", "-"):
        t_sel = t_minus
    else:
        raise ValueError("branch must be 'outer' or 'inner'.")

    t[valid] = t_sel[valid]

    # Coordonnées 3D (après rotation Rx(phi))
    x = R2 * np.cos(theta)
    y = R2 * np.sin(theta) * c - t * s
    z = R2 * np.sin(theta) * s + t * c
    xyz = np.stack([x, y, z], axis=1)

    # Développé cylindre 2 : (u2, v2) = (R2*theta, t)
    u2 = R2 * theta
    v2 = t
    uv = np.stack([u2, v2], axis=1)

    return IntersectionResult(theta=theta, t=t, xyz=xyz, uv=uv, valid_mask=valid)


# ============================
# Développé « gueule de loup » (cylindre 1)
# ============================

def unwrap_on_cylinder1(res: IntersectionResult, R1: float):
    """
    Déroulé de l'intersection sur le cylindre 1 (vertical, axe Oz).
    Coord cylindriques locales : (alpha, z) -> (u1, v1) = (R1*alpha, z)
    Retourne des vecteurs triés par u1 pour un gabarit « z(u) » propre.
    """
    xyz = res.xyz[res.valid_mask]
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

    alpha = np.arctan2(y, x)     # angle autour de l'axe Oz
    alpha_u = np.unwrap(alpha)   # supprime le saut 2π
    u1 = R1 * alpha_u
    v1 = z

    order = np.argsort(u1)
    return u1[order], v1[order]


# ============================
# Exports CSV
# ============================

def export_csv_cyl2(res: IntersectionResult, path_csv: str):
    """Exporte le gabarit à plat sur cylindre 2 (u2, v2) + la courbe 3D (x,y,z)."""
    valid = res.valid_mask
    uv2 = res.uv[valid]
    xyz = res.xyz[valid]
    header = "u2,v2,x,y,z"
    data = np.column_stack([uv2, xyz])
    np.savetxt(path_csv, data, delimiter=",", header=header, comments='', fmt="%.8f")


def export_csv_both(res: IntersectionResult, R1: float, path_csv: str):
    """
    Exporte un CSV combinant :
      - u1,v1 : gueule de loup (cylindre 1)
      - u2,v2 : gabarit cylindre 2 (rééchantillonné sur la longueur de u1 pour une table commune)
    """
    valid = res.valid_mask
    u2 = res.uv[valid, 0]
    v2 = res.uv[valid, 1]
    u1, v1 = unwrap_on_cylinder1(res, R1)

    # Rééchantillonne u2,v2 sur la longueur de u1 pour colonnes alignées
    # (Optionnel ; sinon, écrire deux fichiers séparés)
    if len(u1) < 2:
        raise RuntimeError("Développé cylindrique 1 insuffisant pour rééchantillonner.")
    u2q = np.linspace(u2.min(), u2.max(), len(u1))
    # Tri (u2, v2) par u2 pour interpolation monotone
    order2 = np.argsort(u2)
    u2_sorted = u2[order2]
    v2_sorted = v2[order2]
    v2q = np.interp(u2q, u2_sorted, v2_sorted)

    header = "u1,v1,u2,v2"
    data = np.column_stack([u1, v1, u2q, v2q])
    np.savetxt(path_csv, data, delimiter=",", header=header, comments='', fmt="%.8f")


# ============================
# Visualisation
# ============================

# --- MISE À JOUR : pas de remplissage « gueule de loup » + cadrage auto de la vue 3D

# --- STOPPER LE TUBE 2 À L'INTERSECTION (CLIPPING PAR LE CYLINDRE 1)
# Remplace ta fonction plot_results par celle-ci.

def plot_results(res: IntersectionResult,
                 R1: float,
                 R2: float,
                 phi: float,
                 title_prefix: str = "",
                 pad_factor: float = 0.35,
                 n_ang: int = 120,
                 n_ax: int = 140,
                 alpha_blue: float = 0.8,
                 alpha_orange: float = 0.25,
                 auto_view: bool = True,
                 elev_default: float = 22,
                 azim_default: float = 35,
                 figsize: tuple = (13, 9),
                 clip_tube2_at_cyl1: bool = True,
                 clip_side: str = "outside",
                 clip_eps: float = 1e-9):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    valid_pts = res.xyz[res.valid_mask]
    if valid_pts.size == 0:
        raise RuntimeError("Pas de points d'intersection valides : ajustez R1, R2, phi.")

    # --- bornes géométriques
    x_min, y_min, z_min = valid_pts.min(axis=0)
    x_max, y_max, z_max = valid_pts.max(axis=0)
    t_valid = res.t[res.valid_mask]
    tmin, tmax = np.nanmin(t_valid), np.nanmax(t_valid)

    diag = max(R1, R2)
    pad_z = pad_factor * max(z_max - z_min, diag)
    pad_t = pad_factor * max(tmax - tmin, diag)

    # --- maillages
    th = np.linspace(0, 2*np.pi, n_ang)
    zc = np.linspace(z_min - pad_z, z_max + pad_z, n_ax)
    TH, ZC = np.meshgrid(th, zc)

    # cylindre 1
    X1 = R1*np.cos(TH); Y1 = R1*np.sin(TH); Z1 = ZC

    # cylindre 2 incliné (Rx(phi))
    s, c = np.sin(phi), np.cos(phi)
    TT, ZZ = np.meshgrid(th, np.linspace(tmin - pad_t, tmax + pad_t, n_ax))
    X2 = R2*np.cos(TT)
    Y2 = R2*np.sin(TT)*c - ZZ*s
    Z2 = R2*np.sin(TT)*s + ZZ*c

    # clipping visuel du cylindre 2
    if clip_tube2_at_cyl1:
        r2_sq = X2**2 + Y2**2
        if clip_side.lower() == "outside":
            mask = r2_sq < (R1**2 - clip_eps)
        elif clip_side.lower() == "inside":
            mask = r2_sq > (R1**2 + clip_eps)
        else:
            raise ValueError("clip_side doit valoir 'outside' ou 'inside'")
        X2, Y2, Z2 = np.ma.array(X2, mask=mask), np.ma.array(Y2, mask=mask), np.ma.array(Z2, mask=mask)

    # --- figure
    fig = plt.figure(figsize=figsize)
    fig.set_constrained_layout(True)
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.1, 1.0], height_ratios=[1, 1],
                           wspace=0.25, hspace=0.22)

    # 3D
    ax3d = fig.add_subplot(gs[:, 0], projection='3d')
    ax3d.plot_surface(X1, Y1, Z1, alpha=alpha_blue, color="#7ec8e3", linewidth=0)
    ax3d.plot_surface(X2, Y2, Z2, alpha=alpha_orange, color="#f0a35e", linewidth=0)

    pts = valid_pts
    ax3d.plot(pts[:,0], pts[:,1], pts[:,2], lw=3, color="red", label="Intersection")

    xr = (x_min - diag, x_max + diag)
    yr = (y_min - diag, y_max + diag)
    zr = (z_min - pad_z, z_max + pad_z)
    ax3d.set_xlim(*xr); ax3d.set_ylim(*yr); ax3d.set_zlim(*zr)
    ax3d.set_title(f"{title_prefix}Intersection 3D\nR1={R1}, R2={R2}, phi={phi:.3f} rad")
    ax3d.set_xlabel("x"); ax3d.set_ylabel("y"); ax3d.set_zlabel("z")
    ax3d.set_box_aspect([1,1,1]); ax3d.legend(loc="upper left")

    if auto_view:
        cx, cy, cz = pts.mean(axis=0)
        azim = np.degrees(np.arctan2(cy, cx)) + 35.0
        elev = 20.0 + 10.0*np.clip((cz - zr[0])/(zr[1]-zr[0]), 0, 1)
        ax3d.view_init(elev=elev, azim=azim)
    else:
        ax3d.view_init(elev=elev_default, azim=azim_default)

    # Développé cylindre 2 (haut droite)
    ax2_top = fig.add_subplot(gs[0, 1])
    uv2 = res.uv[res.valid_mask]
    ax2_top.plot(uv2[:,0], uv2[:,1], lw=2, color="red")
    ax2_top.set_title("Développé sur cylindre 2 (incliné) – gabarit tube 2")
    ax2_top.set_xlabel("u₂ = R2·θ (longueur déroulée)")
    ax2_top.set_ylabel("v₂ = t (axe cyl. 2)")
    ax2_top.grid(True, linestyle=':')

    # Développé cylindre 1 (bas droite) — CONTOUR UNIQUEMENT
    ax2_bot = fig.add_subplot(gs[1, 1])
    u1, v1 = unwrap_on_cylinder1(res, R1)
    # on trace un contour fermé SANS remplissage
    u_out = np.r_[u1, u1[0]]
    v_out = np.r_[v1, v1[0]]
    ax2_bot.plot(u_out, v_out, lw=2, color="red")
    ax2_bot.set_title("Développé sur cylindre 1 (vertical) – « gueule de loup »")
    ax2_bot.set_xlabel("u₁ = R1·α (longueur déroulée)")
    ax2_bot.set_ylabel("v₁ = z (axe cyl. 1)")
    ax2_bot.grid(True, linestyle=':')

    plt.show()

# ============================
# Exemple d'exécution
# ============================

if __name__ == "__main__":
    # Paramètres
    R1 = 1.0                 # rayon cylindre 1 (vertical)
    R2 = 1.0                # rayon cylindre 2 (incliné)
    phi = np.deg2rad(90.0)   # angle d'inclinaison autour de l'axe x (en radians)
    branch = "outer"         # "outer" ou "inner" (choix de la racine)

    # Calcul de l'intersection
    res = intersect_cylinders(R1=R1, R2=R2, phi=phi, n_samples=2400, branch=branch)

    # Visualisation : 3D + deux développés
    plot_results(res, R1, R2, phi, alpha_blue=0.8, alpha_orange=0.25,
             auto_view=True)

    # Exports CSV
    export_csv_cyl2(res, "gabarit_cylindre2.csv")
    export_csv_both(res, R1=R1, path_csv="gabarits_developpes.csv")
    print("Exports : gabarit_cylindre2.csv, gabarits_developpes.csv")
