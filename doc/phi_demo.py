import numpy as np

# ---------------------------------------------------------------
# Worked 3-generation tree for the phi-tolerance illustration.
# Hagen-Poiseuille conductance k = pi R^4 / (8 mu L).
# We build a small binary-ish arterial path from an inlet (root)
# down to a tumor terminal, break one segment at a time, and
# compute phi = (1/k_break) / R_tot(path) and the resulting
# first-order tolerance.
#
# Units: work in consistent SI-ish; only ratios (phi) matter, so
# absolute scale of mu, pressures cancels. We keep mu explicit.
# ---------------------------------------------------------------

mu = 3.5e-3            # Pa.s, blood
def k(R, L):           # conductance, R,L in meters
    return np.pi * R**4 / (8.0 * mu * L)

# Path from injection to tumor: 3 generations, radius tapers,
# length per segment roughly constant. Radii in mm -> m.
# Gen boundaries chosen to reflect: lobar artery (large) ->
# segmental -> subsegmental feeding the tumor bed.
mm = 1e-3
segments = [
    # name,           R(mm), L(mm)
    ("g1_lobar",       1.20,  15.0),   # proximal, large radius
    ("g2_segmental",   0.70,  12.0),
    ("g3_subseg",      0.40,  10.0),   # distal, small radius, near tumor
]

# terminal tumor-bed lumped resistance (downstream of last node),
# modeled as a sink conductance gamma_a/mu style closure.
# Pick so the tumor bed is a modest fraction of total (realistic:
# distal vessels + bed dominate). We'll express as a resistance.
R_bed = None  # set below relative to computed segment resistances

# series resistances of each path segment
res = []
for name,Rmm,Lmm in segments:
    ki = k(Rmm*mm, Lmm*mm)
    res.append((name, 1.0/ki, ki))
    print(f"{name:14s} R={Rmm:.2f}mm L={Lmm:.1f}mm  k={ki:.3e}  1/k={1.0/ki:.3e}")

# choose bed resistance = 40% of the summed segment resistance,
# so distal vessels + bed dominate the path (typical: small distal
# vessels carry most resistance).
Rseg_sum = sum(r for _,r,_ in res)
R_bed = 0.40 * Rseg_sum
print(f"\nsum segment R = {Rseg_sum:.3e},  R_bed = {R_bed:.3e}")

# Upstream resistance feeding the path root: proximal aorta/celiac
# is large-bore -> very low resistance; set to 5% of seg sum.
R_up_root = 0.05 * Rseg_sum

R_tot = R_up_root + Rseg_sum + R_bed
print(f"R_up_root = {R_up_root:.3e},  R_tot(path) = {R_tot:.3e}\n")

# For a break in segment i, phi_i = (1/k_i)/R_tot
print("Per-break resistance fraction phi and first-order tolerance:")
print(f"{'break at':16s} {'phi':>8s}  {'delta for 10% flow err':>24s}")
rows=[]
for name, r_i, k_i in res:
    phi = r_i / R_tot
    # flow err ~ phi*delta  =>  delta for 10% err:
    delta_10 = 0.10/phi
    rows.append((name, phi, delta_10))
    print(f"{name:16s} {phi:8.3f}  {delta_10*100:20.0f}% ")

print("\nInterpretation check:")
print(f"  proximal (g1) phi = {rows[0][1]:.3f}  -> tolerance wide")
print(f"  distal   (g3) phi = {rows[2][1]:.3f}  -> tolerance narrow")
print(f"  ratio of R^4: (1.20/0.40)^4 = {(1.20/0.40)**4:.1f}")

# Verify exact (nonlinear) flow ratio for a mis-estimated bridge:
# q_broken/q_intact = R_tot(k_ab)/R_tot(k_br), k_br=k_ab*(1+delta)
def flow_ratio(r_i, delta):
    # R_tot with segment i replaced by 1/(k_i(1+delta)) = r_i/(1+delta)
    R_new = R_tot - r_i + r_i/(1+delta)
    return R_tot / R_new

print("\nExact flow ratio q_broken/q_intact for delta=+0.5 (bridge 50% too conductive):")
for name, r_i, k_i in res:
    phi = r_i/R_tot
    fr = flow_ratio(r_i, 0.5)
    approx = 1 + phi*0.5/(1+0.5)
    print(f"  {name:16s} exact={fr:.4f}  1st-order={approx:.4f}  phi={phi:.3f}")

print("\nExact flow ratio for delta=-0.9 (bridge with near-tissue, 90% too resistive):")
for name, r_i, k_i in res:
    phi = r_i/R_tot
    fr = flow_ratio(r_i, -0.9)
    print(f"  {name:16s} q_broken/q_intact={fr:.4f}  (phi={phi:.3f})")
