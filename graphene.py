import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# 1. Parameters and geometry

a = 0.246            # lattice constant (nm)
d = a / np.sqrt(3)   # nearest-neighbour distance (nm), ~0.142 nm
t = 2.7              # hopping energy (eV), a typical fitted value
hbar = 6.582e-16     # eV * s

# vectors from an A atom to its three B neighbours (nm)
deltas = d * np.array([[0.0, 1.0],
                       [np.sqrt(3) / 2, -0.5],
                       [-np.sqrt(3) / 2, -0.5]])

# lattice vectors (each has length a, 60 degrees apart)
a1 = deltas[0] - deltas[1]
a2 = deltas[0] - deltas[2]

# high-symmetry points in the Brillouin zone (1/nm)
G = np.array([0.0, 0.0])
K = np.array([4 * np.pi / (3 * a), 0.0])
M = np.array([np.pi / a, np.pi / (np.sqrt(3) * a)])



# 2. Hamiltonian and eigenvalues

def f(kx, ky):
    """Structure factor f(k) = sum_j exp(i k . d_j)."""
    k = np.stack([np.asarray(kx, float), np.asarray(ky, float)], axis=-1)
    return np.sum(np.exp(1j * (k @ deltas.T)), axis=-1)


def hamiltonian(kx, ky, t=t, delta=0.0):
    """2x2 Bloch Hamiltonian at every k. Shape: (..., 2, 2)."""
    fk = f(kx, ky)
    H = np.zeros(fk.shape + (2, 2), dtype=complex)
    H[..., 0, 0] = delta
    H[..., 1, 1] = -delta
    H[..., 0, 1] = -t * fk
    H[..., 1, 0] = -t * np.conj(fk)
    return H


def bands(kx, ky, **kw):
    """Eigenvalues sorted ascending: [..., 0] valence, [..., 1] conduction."""
    return np.linalg.eigvalsh(hamiltonian(kx, ky, **kw))



# 3. Path through the Brillouin zone: Gamma -> M -> K -> Gamma

def k_path(points, n=300):
    segs = [np.linspace(p, q, n, endpoint=False)
            for p, q in zip(points[:-1], points[1:])]
    ks = np.vstack(segs + [points[-1][None, :]])
    dist = np.concatenate([[0.0],
                           np.cumsum(np.linalg.norm(np.diff(ks, axis=0), axis=1))])
    nodes = [0.0]
    for p, q in zip(points[:-1], points[1:]):
        nodes.append(nodes[-1] + np.linalg.norm(q - p))
    return ks, dist, nodes


ks, dist, nodes = k_path([G, M, K, G])
E_path = bands(ks[:, 0], ks[:, 1])


# 4. Sanity checks against analytic results

print("Sanity checks (nearest-neighbour model, t = %.2f eV)" % t)
for name, kp, expected in [("Gamma", G, "+/- 3t = +/- %.2f" % (3 * t)),
                           ("M", M, "+/- t  = +/- %.2f" % t),
                           ("K", K, "0")]:
    e = bands(kp[0], kp[1])
    print("  E(%-5s) = %+.4f, %+.4f eV   expected %s" % (name, e[0], e[1], expected))

# Fermi velocity: slope of E vs |q| near K, compared with analytic (3/2) t d / hbar 

q = 0.05                                  # 1/nm, small displacement from K
slope = bands(K[0] + q, K[1])[1] / q      # eV * nm
v_num = slope / hbar * 1e-9               # m/s
v_ana = 1.5 * t * d / hbar * 1e-9
print("  Fermi velocity: numeric %.3e m/s, analytic %.3e m/s (c/%.0f)"
      % (v_num, v_ana, 2.998e8 / v_ana))

# Is the gap zero only at K? Scan a grid over the zone.
gx, gy = np.meshgrid(np.linspace(-20, 20, 401), np.linspace(-20, 20, 401))
Eg = bands(gx, gy)
gap = Eg[..., 1] - Eg[..., 0]
i = np.unravel_index(np.argmin(gap), gap.shape)
print("  Smallest gap on grid: %.4f eV at k = (%.2f, %.2f) 1/nm" % (gap[i], gx[i], gy[i]))


# 5. Figures

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

# (a) bands along the high-symmetry path

ax[0].plot(dist, E_path[:, 0], color="C0")
ax[0].plot(dist, E_path[:, 1], color="C0")
for x in nodes:
    ax[0].axvline(x, color="gray", lw=0.6)
ax[0].set_xticks(nodes)
ax[0].set_xticklabels([r"$\Gamma$", "M", "K", r"$\Gamma$"])
ax[0].set_xlim(0, dist[-1])
ax[0].set_ylabel("Energy (eV)")
ax[0].set_title("(a) Bands along $\\Gamma$-M-K-$\\Gamma$")
ax[0].axhline(0, color="k", lw=0.5, ls=":")

# (b) conduction band over the plane, with the Brillouin zone 

lim = 1.25 * np.linalg.norm(K)
kx, ky = np.meshgrid(np.linspace(-lim, lim, 500), np.linspace(-lim, lim, 500))
Ec = bands(kx, ky)[..., 1]
im = ax[1].pcolormesh(kx, ky, Ec, shading="auto", cmap="viridis")
ang = np.deg2rad(np.arange(0, 361, 60))
ax[1].plot(np.linalg.norm(K) * np.cos(ang), np.linalg.norm(K) * np.sin(ang),
           "w-", lw=1.5)
ax[1].set_aspect("equal")
ax[1].set_xlabel("$k_x$ (1/nm)")
ax[1].set_ylabel("$k_y$ (1/nm)")
ax[1].set_title("(b) Conduction band; white = Brillouin zone")
fig.colorbar(im, ax=ax[1], label="E (eV)", fraction=0.046)

# (c) zoom at K: linear Dirac cone vs gapped version

qs = np.linspace(-1.5, 1.5, 400)
for delta, style, lab in [(0.0, "-", "$\\Delta$ = 0 (graphene)"),
                          (0.3, "--", "$\\Delta$ = 0.3 eV (gapped)")]:
    e = bands(K[0] + qs, np.full_like(qs, K[1]), delta=delta)
    ax[2].plot(qs, e[:, 1], "C1" + style, label=lab)
    ax[2].plot(qs, e[:, 0], "C1" + style)
ax[2].plot(qs, 1.5 * t * d * np.abs(qs), "k:", lw=1, label="linear: $\\hbar v_F |q|$")
ax[2].plot(qs, -1.5 * t * d * np.abs(qs), "k:", lw=1)
ax[2].set_xlabel("q along $k_x$ from K (1/nm)")
ax[2].set_ylabel("Energy (eV)")
ax[2].set_title("(c) Zoom at K: Dirac cone and a gap")
ax[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig("graphene_bands.png", dpi=150)
print("Saved graphene_bands.png")
