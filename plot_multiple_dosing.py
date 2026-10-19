"""Plot 2: 500 mg oral every 12 h for 7 days, steady-state Cmax/Cmin vs therapeutic window."""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pk_model import DOSE, V, CL, KA, F, elimination_rate

os.makedirs("plots", exist_ok=True)

TAU = 12.0            # dosing interval (h)
N_DOSES = 14          # 7 days x 2 doses/day
WINDOW = (3.0, 15.0)  # illustrative therapeutic window (mg/L), NOT a real drug

ke = elimination_rate()


def rhs(t, y):
    gut, central = y
    return [-KA * gut, KA * gut - ke * central]


# Simulate dose by dose: add each dose to the gut, then integrate to the next dose
state = np.array([0.0, 0.0])
t_parts, c_parts = [], []
for i in range(N_DOSES):
    state[0] += F * DOSE
    t_eval = np.linspace(i * TAU, (i + 1) * TAU, 241)
    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), state,
                    t_eval=t_eval, rtol=1e-10, atol=1e-12)
    t_parts.append(sol.t)
    c_parts.append(sol.y[1] / V)
    state = sol.y[:, -1]

t = np.concatenate(t_parts)
c = np.concatenate(c_parts)

# Steady-state values from the last dosing interval
last = t >= (N_DOSES - 1) * TAU
cmax_ss = c[last].max()
cmin_ss = c[last].min()

# Check against the closed-form steady-state equation
tt = np.linspace(0, TAU, 2401)
c_ss_formula = (F * DOSE * KA / (V * (KA - ke))) * (
    np.exp(-ke * tt) / (1 - np.exp(-ke * TAU))
    - np.exp(-KA * tt) / (1 - np.exp(-KA * TAU)))

# Is it really at steady state? Compare the last two intervals
prev = (t >= (N_DOSES - 2) * TAU) & (t < (N_DOSES - 1) * TAU)
ss_gap = abs(c[prev].max() - cmax_ss)

print(f"Steady-state Cmax (ODE)     : {cmax_ss:.2f} mg/L")
print(f"Steady-state Cmin (ODE)     : {cmin_ss:.2f} mg/L")
print(f"Steady-state Cmax (formula) : {c_ss_formula.max():.2f} mg/L")
print(f"Steady-state Cmin (formula) : {c_ss_formula.min():.2f} mg/L")
print(f"Max difference ODE vs formula: {abs(cmax_ss - c_ss_formula.max()):.2e} mg/L")
print(f"Cmax change over last interval: {ss_gap:.2e} mg/L")
print(f"Window: {WINDOW[0]}-{WINDOW[1]} mg/L -> "
      f"{'inside' if cmin_ss >= WINDOW[0] and cmax_ss <= WINDOW[1] else 'OUTSIDE'}")

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
ax.axhspan(*WINDOW, color="tab:green", alpha=0.15, label="Therapeutic window (illustrative)")
ax.plot(t / 24, c, color="tab:blue", lw=2, label="Concentration")
ax.plot(np.arange(N_DOSES) * TAU / 24, np.zeros(N_DOSES), "^", color="k", ms=6, label="Dose")
ax.axhline(cmax_ss, color="tab:red", ls="--", lw=1)
ax.axhline(cmin_ss, color="tab:purple", ls="--", lw=1)
ax.text(0.1, cmax_ss + 0.3, f"Css,max = {cmax_ss:.1f} mg/L", color="tab:red")
ax.text(0.1, cmin_ss - 1.0, f"Css,min = {cmin_ss:.1f} mg/L", color="tab:purple")

ax.set_xlabel("Time (days)")
ax.set_ylabel("Concentration (mg/L)")
ax.set_title(f"{DOSE:.0f} mg oral every {TAU:.0f} h for 7 days")
ax.set_ylim(bottom=-1)
ax.legend(loc="upper right")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("plots/02_multiple_dosing.png", dpi=150)
print("Saved plots/02_multiple_dosing.png")
plt.show()