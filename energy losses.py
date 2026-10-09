from atmosphere_path_lenght_X import Ray
import numpy as np
import matplotlib.pyplot as plt

M_MU = 105.65839  # muon mass [MeV]

#################### p <-> T ######################
def T_from_p(p_MeV):
    p = np.asarray(p_MeV, dtype=float)
    return p ** 2 / (np.sqrt(p ** 2 + M_MU ** 2) + M_MU)

def p_from_T(T_MeV):
    T = np.asarray(T_MeV, dtype=float)
    return np.sqrt(T ** 2 + 2.0 * T * M_MU)


#################### GROOM TABLE ######################
def load_groom(path="table_Groom.txt", T_min=10.0):
    T, p, dEdx = [], [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 12:
                continue
            try:
                v = [float(x) for x in parts]
            except ValueError:  # title
                continue
            if v[0] < T_min:
                continue
            T.append(v[0])
            p.append(v[1])
            dEdx.append(v[7])  # (dE/dx)
    return np.array(T), np.array(p), np.array(dEdx)

T_TAB, P_TAB, DEDX_TAB = load_groom()

def print_groom_table():
    print(f"{'i':>4} {'T [MeV]':>12} {'p [MeV/c]':>12} {'dE/dx [MeV cm2/g]':>18}")
    for i, (t, p, d, r) in enumerate(zip(T_TAB, P_TAB, DEDX_TAB)):
        print(f"{i:>4} {t:>12.4e} {p:>12.4e} {d:>18.4e}")

def plot_dedx_vs_p(save="dedx_vs_p.png"):

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 16,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "legend.fontsize": 13,
        "mathtext.fontset": "cm",
        "mathtext.rm": "serif",
        "mathtext.it": "serif:italic",
        "mathtext.bf": "serif:bold"
    })

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(P_TAB / 1000, DEDX_TAB, "k.-", ms=4, lw=1)

    i_min = np.argmin(DEDX_TAB)

    ax.plot( P_TAB[i_min] / 1000,DEDX_TAB[i_min],"ro",
             label=f"minimum: {DEDX_TAB[i_min]:.2f} MeV·cm²/g "f"at p={P_TAB[i_min] / 1000:.2f} GeV/c")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$p$ [GeV/c]")
    ax.set_ylabel(r"$dE/dx$ [MeV cm$^2$/g]")
    ax.set_title("Total muon energy loss in air (Groom)")

    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(save, dpi=150)
    plt.show()

#################### ENERGY LOSSES X ######################
def dedx_groom(T_MeV):
    T = np.asarray(T_MeV, dtype=float)
    logT = np.log(T)
    logT_tab = np.log(T_TAB)
    logD_tab = np.log(DEDX_TAB)

    return np.exp(np.interp(logT, logT_tab, logD_tab))

def energy_loss_dedx(T0_MeV, X, steps=20000):
    T = float(T0_MeV)
    if X <= 0:
        return T
    dX = X / steps

    for _ in range(steps):

        if T <= T_TAB[0]:
            return np.nan

        dEdX = dedx_groom(T)
        T -= dEdX * dX
    if T <= T_TAB[0]:
        return np.nan
    return T

def loss_for_ray(theta_deg, T0_MeV, h_prod=15.0, h_obs=0.0):
    ray = Ray(
        theta_deg,
        h_obs=h_obs,
        h_prod=h_prod
    )

    X = ray.X_total
    T_final = energy_loss_dedx(T0_MeV,X)
    delta_T = T0_MeV - T_final
    return X, T_final, delta_T

def energy_at_production(theta_deg, T_ground_MeV, h_obs=0.0, h_prod=15.0, steps=20000, X=None):

    if X is None:
        ray = Ray(theta_deg, h_obs=h_obs, h_prod=h_prod)
        X = ray.X_total

    T = float(T_ground_MeV)
    if X <= 0:
        return T
    dX = -X / steps

    for _ in range(steps):
        dEdX = dedx_groom(T)
        T -= dEdX * dX

    return T

def dTprod_dTground(theta_deg, T_ground_MeV, h_obs=0.0, h_prod=15.0, steps=20000, rel_step=1e-3):
        #dTp/dTg ≈ [Tp(Tg(1+h)) - Tp(Tg(1-h))] / (2 h Tg),   h = rel_step
    X = Ray(theta_deg, h_obs=h_obs, h_prod=h_prod).X_total
    T_hi = T_ground_MeV * (1.0 + rel_step)
    T_lo = T_ground_MeV * (1.0 - rel_step)
    Tp_hi = energy_at_production(theta_deg, T_hi, steps=steps, X=X)
    Tp_lo = energy_at_production(theta_deg, T_lo, steps=steps, X=X)
    return (Tp_hi - Tp_lo) / (T_hi - T_lo)

T0 = 10e3  #  [MeV]
#print_groom_table()
#plot_dedx_vs_p()

#X, T_final, delta_T = loss_for_ray(theta_deg=0, T0_MeV=T0, h_prod=15)
#T_at_prod = energy_at_production(theta_deg=0, T_ground_MeV=T_final, h_prod=15)

#print(f"X           = {X:.1f} g/cm²")
#print(f"T_final     = {T_final/1000:.5f} GeV")
#print(f"dE          = {delta_T/1000:.5f} GeV")
#print(f"T_production = {T_at_prod/1000:.5f} GeV")
