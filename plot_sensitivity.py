"""Plot 3: how clearance, volume and absorption rate change Cmax, Tmax, AUC and half-life."""
import os
import numpy as np
import matplotlib.pyplot as plt
from pk_model import oral_analytical, DOSE, V, CL, KA, F

os.makedirs("plots", exist_ok=True)

FACTORS = [0.5, 1.0, 1.5]          # 50%, baseline, 150%
T = np.linspace(0, 400, 40001)     # long horizon so AUC reaches ~infinity


def metrics(cl, v, ka):
    ke = cl / v
    c = oral_analytical(T, v=v, cl=cl, ka=ka)
    tmax = np.log(ka / ke) / (ka - ke)
    cmax = oral_analytical(np.array([tmax]), v=v, cl=cl, ka=ka)[0]
    auc = np.sum((c[1:] + c[:-1]) / 2 * np.diff(T))      # trapezoid rule
    half_life = np.log(2) / ke
    return cmax, tmax, auc, half_life


params = {
    "Clearance (CL)":       lambda f: dict(cl=CL * f, v=V, ka=KA),
    "Volume (V)":           lambda f: dict(cl=CL, v=V * f, ka=KA),
    "Absorption rate (ka)": lambda f: dict(cl=CL, v=V, ka=KA * f),
}

print(f"Baseline: CL={CL} L/h, V={V} L, ka={KA} 1/h, dose={DOSE} mg")
print(f"Analytical AUC = F*Dose/CL = {F * DOSE / CL:.1f} mg*h/L\n")
print(f"{'Parameter':<22}{'x':>5}{'Cmax':>9}{'Tmax':>8}{'AUC':>9}{'t1/2':>8}")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
t_plot = np.linspace(0, 48, 481)

for ax, (name, make) in zip(axes, params.items()):
    for f in FACTORS:
        kw = make(f)
        cmax, tmax, auc, hl = metrics(**kw)
        print(f"{name:<22}{f:>5.1f}{cmax:>9.2f}{tmax:>8.2f}{auc:>9.1f}{hl:>8.2f}")
        ax.plot(t_plot, oral_analytical(t_plot, **kw), lw=2, label=f"{f:.1f}x")
    ax.set_title(name)
    ax.set_xlabel("Time (h)")
    ax.grid(alpha=0.3)
    ax.legend(title="Change")
    print()

axes[0].set_ylabel("Concentration (mg/L)")
fig.suptitle("Sensitivity of oral PK profile to each parameter (others held at baseline)")
plt.tight_layout()
plt.savefig("plots/03_sensitivity.png", dpi=150)
print("Saved plots/03_sensitivity.png")
plt.show()