# Symbolic + numeric verification of every claim in docs/THEORY.md.
# Each check prints PASS/FAIL; the script exits non-zero on any failure.

import sys
import sympy as sp
import numpy as np

FAIL = []

def check(name, ok):
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)

R1, R2, t, theta, phi, z0, alpha = sp.symbols("R1 R2 t theta phi z0 alpha", real=True, positive=False)
R1p, R2p = sp.symbols("R1p R2p", positive=True)

# --- (1) parameterisation: rotation Rx(phi) applied to (R2 cosθ, R2 sinθ, t)
Rx = sp.Matrix([[1, 0, 0],
                [0, sp.cos(phi), -sp.sin(phi)],
                [0, sp.sin(phi),  sp.cos(phi)]])
P0 = sp.Matrix([R2*sp.cos(theta), R2*sp.sin(theta), t])
P = sp.simplify(Rx * P0)
x, y, z = P
check("(1) x = R2 cosθ", sp.simplify(x - R2*sp.cos(theta)) == 0)
check("(1) y = R2 sinθ cosφ − t sinφ", sp.simplify(y - (R2*sp.sin(theta)*sp.cos(phi) - t*sp.sin(phi))) == 0)
check("(1) z = R2 sinθ sinφ + t cosφ", sp.simplify(z - (R2*sp.sin(theta)*sp.sin(phi) + t*sp.cos(phi))) == 0)
axis = Rx * sp.Matrix([0, 0, 1])
check("(1) axe cyl2 = (0, −sinφ, cosφ)", sp.simplify(axis - sp.Matrix([0, -sp.sin(phi), sp.cos(phi)])) == sp.zeros(3, 1))

# --- (2) quadratic in t and its coefficients
expr = sp.expand(x**2 + y**2 - R1**2)
poly = sp.Poly(expr, t)
a_c, b_c, c_c = poly.all_coeffs()
check("(2) a = sin²φ", sp.simplify(a_c - sp.sin(phi)**2) == 0)
check("(2) b = −2 R2 sinθ cosφ sinφ", sp.simplify(b_c + 2*R2*sp.sin(theta)*sp.cos(phi)*sp.sin(phi)) == 0)
check("(2) c = R2²(cos²θ + sin²θ cos²φ) − R1²",
      sp.simplify(c_c - (R2**2*(sp.cos(theta)**2 + sp.sin(theta)**2*sp.cos(phi)**2) - R1**2)) == 0)

# --- (3) discriminant closed form
disc = sp.simplify(b_c**2 - 4*a_c*c_c)
check("(3) Δ = 4 sin²φ (R1² − R2² cos²θ)",
      sp.simplify(disc - 4*sp.sin(phi)**2*(R1**2 - R2**2*sp.cos(theta)**2)) == 0)

# --- (4) roots (for sinφ > 0)
phi_pos = sp.Symbol("phi", positive=True)  # 0 < φ < π ⇒ sinφ > 0 checked numerically below
tp = (R2*sp.sin(theta)*sp.cos(phi) + sp.sqrt(R1**2 - R2**2*sp.cos(theta)**2)) / sp.sin(phi)
tm = (R2*sp.sin(theta)*sp.cos(phi) - sp.sqrt(R1**2 - R2**2*sp.cos(theta)**2)) / sp.sin(phi)
resid_p = expr.subs(t, tp)
resid_m = expr.subs(t, tm)
ok_p = ok_m = True
rng = np.random.default_rng(1)
for _ in range(200):
    vals = {R1: float(rng.uniform(5, 100)), R2: float(rng.uniform(5, 100)),
            theta: float(rng.uniform(0, 2*np.pi)), phi: float(rng.uniform(0.05, np.pi - 0.05))}
    if vals[R2]*abs(np.cos(vals[theta])) > vals[R1]:
        continue  # outside existence domain
    ok_p &= abs(complex(resid_p.subs(vals)).real) < 1e-6
    ok_m &= abs(complex(resid_m.subs(vals)).real) < 1e-6
check("(4) t₊ vérifie x²+y²=R1² (200 tirages)", ok_p)
check("(4) t₋ vérifie x²+y²=R1² (200 tirages)", ok_m)

# --- (5) existence condition Δ ≥ 0 ⟺ R2|cosθ| ≤ R1 (numeric scan)
ok = True
for _ in range(500):
    r1v, r2v, thv, phv = rng.uniform(1, 100), rng.uniform(1, 100), rng.uniform(0, 2*np.pi), rng.uniform(0.05, np.pi-0.05)
    d = 4*np.sin(phv)**2*(r1v**2 - r2v**2*np.cos(thv)**2)
    ok &= (d >= 0) == (r2v*abs(np.cos(thv)) <= r1v)
check("(5) Δ≥0 ⟺ R2|cosθ| ≤ R1", ok)

# --- §3 isometry: first fundamental form of the cylinder chart (θ, t)
Rc = sp.Symbol("R", positive=True)
S = sp.Matrix([Rc*sp.cos(theta), Rc*sp.sin(theta), t])
Su, Sv = S.diff(theta), S.diff(t)
E, F, G = sp.simplify(Su.dot(Su)), sp.simplify(Su.dot(Sv)), sp.simplify(Sv.dot(Sv))
check("(6) ds² = R²dθ² + dt² (E=R², F=0, G=1)", (E, F, G) == (Rc**2, 0, 1))

# --- §3.1 symmetry claims
s_sum = sp.simplify(tp + tm - 2*R2*sp.sin(theta)*sp.cos(phi)/sp.sin(phi))
check("(§3.1) t₊+t₋ = 2R2 sinθ cosφ / sinφ", s_sum == 0)
sym = sp.simplify(tp.subs(theta, -theta) + tm)
check("(§3.1) t₊(−θ) = −t₋(θ)", sym == 0)

# --- (9) transfer relation cosα = (R2/R1) cosθ  (x = R1 cosα on C1, x = R2 cosθ on C2)
# Direct: on the intersection, x is shared; C1 gives x = R1 cosα by definition of α.
# Verified numerically end-to-end with (4):
ok = True
for _ in range(300):
    r1v, r2v = rng.uniform(5, 100), rng.uniform(5, 100)
    thv, phv = rng.uniform(0, 2*np.pi), rng.uniform(0.05, np.pi-0.05)
    if r2v*abs(np.cos(thv)) > r1v:
        continue
    tv = (r2v*np.sin(thv)*np.cos(phv) + np.sqrt(r1v**2 - r2v**2*np.cos(thv)**2))/np.sin(phv)
    xv = r2v*np.cos(thv)
    yv = r2v*np.sin(thv)*np.cos(phv) - tv*np.sin(phv)
    al = np.arctan2(yv, xv)
    ok &= abs(np.cos(al) - (r2v/r1v)*np.cos(thv)) < 1e-9
check("(9) cosα = (R2/R1) cosθ sur la courbe", ok)

# --- (10) perpendicular case closed form v1(α) = ±sqrt(R2² − R1² cos²α)
ok = True
for _ in range(300):
    r1v, r2v = rng.uniform(5, 100), rng.uniform(5, 100)
    thv = rng.uniform(0, 2*np.pi)
    if r2v*abs(np.cos(thv)) > r1v:
        continue
    tv = np.sqrt(r1v**2 - r2v**2*np.cos(thv)**2)  # outer root, φ = π/2
    xv = r2v*np.cos(thv)
    yv = -tv                                       # y = −t sinφ = −t
    zv = r2v*np.sin(thv)                           # z = R2 sinθ sinφ = R2 sinθ
    al = np.arctan2(yv, xv)
    lhs = zv
    inside = r2v**2 - r1v**2*np.cos(al)**2
    ok &= inside >= -1e-9 and abs(abs(lhs) - np.sqrt(max(inside, 0))) < 1e-6
check("(10) |v₁| = sqrt(R2² − R1² cos²α) à 90°", ok)

# --- §4 plane cut: curve, ellipse axes, sinusoid
zc = z0 - R1*sp.tan(phi)*sp.sin(theta)
gamma = sp.Matrix([R1*sp.cos(theta), R1*sp.sin(theta), zc])
n = sp.Matrix([0, sp.sin(phi), sp.cos(phi)])
plane_resid = sp.simplify(n.dot(gamma) - z0*sp.cos(phi))
check("(11/12) γ(θ) est dans le plan", plane_resid == 0)

# Ellipse semi-axes: project the curve into the cutting plane basis
e1 = sp.Matrix([1, 0, 0])                            # in-plane, ⟂ tube axis view
e2 = sp.Matrix([0, sp.cos(phi), -sp.sin(phi)])       # in-plane, steepest direction
center = sp.Matrix([0, 0, z0])
u_coord = sp.simplify(e1.dot(gamma - center))
v_coord = sp.simplify(e2.dot(gamma - center))
# u = R1 cosθ ; v should be R1 sinθ / cosφ (semi-axis R1/cosφ)
check("(§4.1) demi-petit axe = R1", sp.simplify(u_coord - R1*sp.cos(theta)) == 0)
check("(§4.1) demi-grand axe = R1/cosφ",
      sp.simplify(v_coord - R1*sp.sin(theta)/sp.cos(phi)) == 0)
check("(§4.1) (u/R1)² + (v cosφ/R1)² = 1",
      sp.simplify((u_coord/R1)**2 + (v_coord*sp.cos(phi)/R1)**2 - 1) == 0)

# (13) developed curve of the plane cut
u = R1*theta
v13 = z0 - R1*sp.tan(phi)*sp.sin(u/R1)
check("(13) v(u) = z0 − R1 tanφ sin(u/R1)", sp.simplify(v13 - zc.subs(theta, u/R1)) == 0)
check("(13) amplitude crête-à-crête = 2 R1 |tanφ| (numérique)",
      all(abs((z0v := 0) + r1v*abs(np.tan(p)) - max(abs(-r1v*np.tan(p)*np.sin(np.linspace(0, 2*np.pi, 20001))))) < 1e-6
          for r1v, p in [(30.0, 0.6), (75.0, -1.2), (12.5, 1.4)]))

# --- §6.1 Steinmetz: R1 = R2 = R, φ = π/2 → intersection in planes z = ±y
Rs = sp.Symbol("R", positive=True)
tps = tp.subs({R1: Rs, R2: Rs, phi: sp.pi/2})
xs = (Rs*sp.cos(theta))
ys = (Rs*sp.sin(theta)*sp.cos(sp.pi/2) - tps*sp.sin(sp.pi/2))
zs = (Rs*sp.sin(theta)*sp.sin(sp.pi/2) + tps*sp.cos(sp.pi/2))
prod = sp.simplify((zs - ys)*(zs + ys))   # must vanish on the curve
check("(§6.1) Steinmetz : (z−y)(z+y) = 0 sur la courbe", sp.simplify(prod) == 0)

# --- §5 tangency, case R2 > R1: Δ has a SIMPLE zero at cosθ* = R1/R2
#     → dt±/dθ diverges (vertical tangent in the developed plane).
r1v, r2v, phv = 30.0, 50.0, np.pi/3
th_star = np.arccos(r1v/r2v)
def tplus(th, r1=r1v, r2=r2v, ph=phv):
    return (r2*np.sin(th)*np.cos(ph) + np.sqrt(max(r1**2 - r2**2*np.cos(th)**2, 0)))/np.sin(ph)
d_near = (tplus(th_star + 2e-7) - tplus(th_star + 1e-7)) / 1e-7
d_far = (tplus(th_star + 2e-3) - tplus(th_star + 1e-3)) / 1e-3
check("(§5) R2>R1 : pente t₊ diverge au point de tangence (zéro simple)",
      abs(d_near) > 10*abs(d_far) > 0)

# --- §5/§6.1, case R1 = R2 = R: Δ has a DOUBLE zero at θ = 0, π
#     → finite one-sided slopes; the developed curve has a corner with
#       slopes −tan(φ/2) and +cot(φ/2) (in dv/du = t'/R units).
Rv, phv2 = 40.0, 0.9
h = 1e-4  # small enough for the limit, large enough to avoid float cancellation
left = (tplus(0 - h, Rv, Rv, phv2) - tplus(0 - 2*h, Rv, Rv, phv2)) / h / Rv
right = (tplus(0 + 2*h, Rv, Rv, phv2) - tplus(0 + h, Rv, Rv, phv2)) / h / Rv
check("(§5) R1=R2 : pentes unilatérales finies −tan(φ/2) | +cot(φ/2)",
      abs(left + np.tan(phv2/2)) < 1e-6 and abs(right - 1/np.tan(phv2/2)) < 1e-6)

print()
if FAIL:
    print(f"{len(FAIL)} ÉCHEC(S) :", *FAIL, sep="\n  - ")
    sys.exit(1)
print("Toutes les vérifications symboliques/numériques passent.")
