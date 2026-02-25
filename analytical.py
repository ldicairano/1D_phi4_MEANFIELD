import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import logsumexp
import matplotlib.pyplot as plt


# =========================
# PARAMETERS (YOUR VALUES)
# =========================
N   = 4000
c   = 1.0
mu  = 2.0
lam = 3.0/5.0

a4   = lam / 24.0
a2   = 2.0 * c * (N / (N - 1.0)) - mu / 2.0
facN = N / (N - 1.0)


# =========================
# FULL RANGE IN T
# =========================
TMIN = 1e-1
TMAX = 2e2
NPTS = 450


# =========================
# QUADRATURE
# =========================
NGL = 240
x_gl, w_gl = leggauss(NGL)
TAIL = 30.0


# =========================
# ROOT SEARCH IN SIGMA
# =========================
SIGMAX = 120.0
NSIG   = 2001
XTOL_SIG = 1e-11


# =========================
# OUTPUT
# =========================
OUT_CAL = "caloric_full_selected.dat"
OUT_CV  = "cv_full_selected_fd.dat"
OUT_ALL = "all_roots_full.dat"


def cutoff_L(beta, sigma_scale):
    """
    Choose a HALF-WIDTH L for integration that is safe for all sigma in [-sigma_scale, +sigma_scale].
    IMPORTANT: for branch comparison at fixed beta, L MUST NOT depend on sigma.
    """
    beta = float(beta)
    # characteristic sigma scale (use SIGMAX for full scan)
    sigma = float(sigma_scale)

    L4 = (TAIL / max(beta * a4, 1e-14)) ** 0.25
    L2 = (TAIL / max(beta * abs(a2), 1e-14)) ** 0.5
    shift = abs(sigma) / max(2.0 * beta * abs(a2), 1e-14)

    L = 2.0 * max(L2, L4) + 1.8 * shift
    return float(max(L, 12.0))


def moments_and_logz(beta, sigma, L, phi_grid, logw_base):
    """
    Compute logZ and moments using fixed L and precomputed phi_grid.
    Weight: exp[-beta(a2 phi^2 + a4 phi^4) + sigma phi]
    """
    beta = float(beta); sigma = float(sigma)

    expo = -beta * (a2 * phi_grid**2 + a4 * phi_grid**4) + sigma * phi_grid

    # log integrand weights: log(w_i) + expo_i + log(Jacobian=L)
    # logw_base = log(w_gl) already included, Jacobian handled outside
    logw = logw_base + expo
    logS = logsumexp(logw)
    logZ = np.log(L) + logS

    p = np.exp(logw - logS)
    m1 = np.sum(p * phi_grid)
    m2 = np.sum(p * phi_grid**2)
    m4 = np.sum(p * phi_grid**4)
    return logZ, m1, m2, m4


def g_sigma(beta, sigma, L, phi_grid, logw_base):
    _, m1, _, _ = moments_and_logz(beta, sigma, L, phi_grid, logw_base)
    return sigma - 4.0 * beta * c * facN * m1


def energy_per_particle(beta, sigma_star, L, phi_grid, logw_base):
    _, m1, m2, m4 = moments_and_logz(beta, sigma_star, L, phi_grid, logw_base)
    return 1.0/(2.0*beta) + a2*m2 + a4*m4 - 2.0*c*facN*(m1**2)


def free_energy_per_particle(beta, sigma_star, L, phi_grid, logw_base):
    logZ, _, _, _ = moments_and_logz(beta, sigma_star, L, phi_grid, logw_base)
    hs  = (sigma_star**2) / (8.0 * beta * c * facN)
    kin = 0.5 * np.log(2.0*np.pi/beta)
    return -(kin + logZ - hs) / beta


def find_all_sigma_roots(beta, L, phi_grid, logw_base):
    sig = np.linspace(-SIGMAX, SIGMAX, NSIG)
    gvals = np.empty_like(sig)

    for i, s in enumerate(sig):
        gvals[i] = g_sigma(beta, s, L, phi_grid, logw_base)

    roots = []
    for i in range(NSIG - 1):
        a, b = sig[i], sig[i+1]
        fa, fb = gvals[i], gvals[i+1]
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0.0:
            roots.append(a)
        elif fa * fb < 0:
            r = brentq(lambda x: g_sigma(beta, x, L, phi_grid, logw_base),
                       a, b, xtol=XTOL_SIG, maxiter=200)
            roots.append(r)

    roots = np.array(sorted(roots))
    if roots.size == 0:
        return roots, sig, gvals

    uniq = [roots[0]]
    for r in roots[1:]:
        if abs(r - uniq[-1]) > 1e-6:
            uniq.append(r)
    return np.array(uniq), sig, gvals


def cV_finite_diff(T, e):
    T = np.asarray(T); e = np.asarray(e)
    Tmid = 0.5*(T[1:]+T[:-1])
    return Tmid, np.diff(e)/np.diff(T)


def main():
    print("Params:")
    print(f"  N={N}  c={c}  mu={mu}  lambda={lam}")
    print(f"  a2={a2:.10f}  a4={a4:.10f}")
    print(f"  T in [{TMIN}, {TMAX}] on log grid (NPTS={NPTS})")
    print(f"  sigma scan in [-{SIGMAX},{SIGMAX}] with NSIG={NSIG}")

    Tgrid = np.logspace(np.log10(TMIN), np.log10(TMAX), NPTS)

    T_sel, e_sel, m_sel, s_sel, f_sel = [], [], [], [], []
    nrootsT = []
    all_rows = []

    # diagnostic snapshots
    T_tests = [0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0]
    dumped = set()

    for T in Tgrid:
        beta = 1.0/float(T)

        # FIX: use SAME L for ALL sigma at this beta
        L = cutoff_L(beta, SIGMAX)
        phi_grid = L * x_gl
        logw_base = np.log(w_gl)

        roots, sig_grid, gvals = find_all_sigma_roots(beta, L, phi_grid, logw_base)
        nrootsT.append(len(roots))

        # dump g(sigma)
        t0 = min(T_tests, key=lambda x: abs(T-x))
        if abs(T - t0)/t0 < 0.03 and t0 not in dumped:
            plt.figure()
            plt.plot(sig_grid, gvals)
            plt.axhline(0, ls="--")
            plt.title(f"g(sigma) at T={T:.4g} (beta={beta:.4g}), roots={len(roots)}")
            plt.xlabel("sigma"); plt.ylabel("g(sigma)")
            plt.tight_layout()
            plt.savefig(f"g_sigma_T{t0}.png", dpi=200)
            dumped.add(t0)

        if len(roots) == 0:
            T_sel.append(T); e_sel.append(np.nan); m_sel.append(np.nan); s_sel.append(np.nan); f_sel.append(np.nan)
            continue

        # evaluate candidates and select min free energy
        candidates = []
        for s in roots:
            logZ, m1, _, _ = moments_and_logz(beta, s, L, phi_grid, logw_base)
            e = energy_per_particle(beta, s, L, phi_grid, logw_base)
            f = free_energy_per_particle(beta, s, L, phi_grid, logw_base)
            candidates.append((f, e, m1, s))
            all_rows.append((float(T), float(beta), float(s), float(m1), float(e), float(f)))

        candidates.sort(key=lambda x: x[0])
        f0, e0, m0, s0 = candidates[0]
        T_sel.append(T); e_sel.append(e0); m_sel.append(m0); s_sel.append(s0); f_sel.append(f0)

    T_sel = np.array(T_sel)
    e_sel = np.array(e_sel)
    m_sel = np.array(m_sel)
    s_sel = np.array(s_sel)
    f_sel = np.array(f_sel)
    nrootsT = np.array(nrootsT)

    mask = np.isfinite(e_sel)
    T_use = T_sel[mask]
    e_use = e_sel[mask]
    m_use = m_sel[mask]
    s_use = s_sel[mask]
    f_use = f_sel[mask]

    Tmid, cV = cV_finite_diff(T_use, e_use)

    np.savetxt(OUT_CAL, np.column_stack([T_use, e_use, m_use, s_use, f_use]),
               header="T   e(T)   m(T)   sigma*(T)   f(T)  (selected by min free energy)")
    np.savetxt(OUT_CV, np.column_stack([Tmid, cV]),
               header="Tmid   cV(T)=de/dT  (finite-diff on selected branch)")
    np.savetxt(OUT_ALL, np.array(all_rows),
               header="T  beta  sigma_root  m  e  f  (ALL roots before selection)")

    # plots
    plt.figure()
    plt.loglog(T_use, e_use)
    plt.xlabel("T"); plt.ylabel("e(T)=<H>/N (selected)")
    plt.tight_layout()
    plt.savefig("caloric_full_selected.png", dpi=200)

    plt.figure()
    plt.semilogx(Tmid, cV)
    plt.xlabel("T"); plt.ylabel("cV(T)=de/dT (FD, selected)")
    plt.tight_layout()
    plt.savefig("cv_full_selected.png", dpi=200)

    plt.figure()
    plt.semilogx(T_use, np.abs(m_use))
    plt.xlabel("T"); plt.ylabel("|m(T)| (selected)")
    plt.tight_layout()
    plt.savefig("m_full_selected.png", dpi=200)

    plt.figure()
    plt.semilogx(T_sel, nrootsT, marker="o", ms=2, lw=1)
    plt.xlabel("T"); plt.ylabel("number of sigma roots of g(sigma)")
    plt.tight_layout()
    plt.savefig("nroots_vs_T.png", dpi=200)

    print("Saved:")
    print(f"  {OUT_CAL}")
    print(f"  {OUT_CV}")
    print(f"  {OUT_ALL}")
    print("  caloric_full_selected.png, cv_full_selected.png, m_full_selected.png, nroots_vs_T.png")
    print("  plus some g_sigma_T*.png diagnostics")

    plt.show()


if __name__ == "__main__":
    main()
