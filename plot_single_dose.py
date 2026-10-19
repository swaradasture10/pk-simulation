"""Plot 1: single-dose concentration-time curves (IV vs oral), ODE vs analytical."""
import os
import numpy as np
import matplotlib.pyplot as plt
from pk_model import (iv_analytical, oral_analytical, iv_ode, oral_ode,
                      DOSE, V, CL, KA, F, elimination_rate)

os.makedirs("plots", exist_ok=True)

t = np.linspace(0, 48, 481)
ke = elimination_rate()

# Oral Tmax and Cmax (closed form)
tmax = np.log(KA / ke) / (KA - ke)
cmax = oral_analytical(np.array([tmax]))[0]

fig, ax = plt.subplots(figsize=(8, 5))

# ODE solutions (solid lines)
ax.plot(t, iv_ode(t), color="tab:blue", lw=2, label="IV bolus (ODE)")
ax.plot(t, oral_ode(t), color="tab:orange", lw=2, label="Oral (ODE)")

# Analytical solutions (dots, every 20th point so they don't hide the lines)
ax.plot(t[::20], iv_analytical(t[::20]), "k.", ms=6, label="Analytical")
ax.plot(t[::20], oral_analytical(t[::20]), "k.", ms=6)

# Mark oral Cmax
ax.plot(tmax, cmax, "o", color="tab:red", ms=8, zorder=5)
ax.annotate(f"Oral Cmax = {cmax:.2f} mg/L\nTmax = {tmax:.2f} h",
            xy=(tmax, cmax), xytext=(tmax + 6, cmax + 1.2),
            arrowprops=dict(arrowstyle="->"))

ax.set_xlabel("Time (h)")
ax.set_ylabel("Concentration (mg/L)")
ax.set_title(f"One-compartment PK, {DOSE:.0f} mg single dose\n"
             f"V = {V:.0f} L, CL = {CL:.0f} L/h, ka = {KA:.1f} 1/h")
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("plots/01_single_dose.png", dpi=150)
print(f"Oral Tmax = {tmax:.2f} h, Cmax = {cmax:.2f} mg/L")
print("Saved plots/01_single_dose.png")
plt.show()