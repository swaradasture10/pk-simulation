"""One-compartment PK model: IV bolus and first-order oral absorption."""
import numpy as np
from scipy.integrate import solve_ivp

# Illustrative parameters (NOT a real drug)
DOSE = 500.0   # mg
V = 50.0       # L, volume of distribution
CL = 5.0       # L/h, clearance
KA = 1.0       # 1/h, absorption rate constant
F = 1.0        # bioavailability


def elimination_rate(cl=CL, v=V):
    """ke = CL / V  (1/h)"""
    return cl / v


# ---------- Analytical solutions ----------
def iv_analytical(t, dose=DOSE, v=V, cl=CL):
    ke = elimination_rate(cl, v)
    return (dose / v) * np.exp(-ke * t)


def oral_analytical(t, dose=DOSE, v=V, cl=CL, ka=KA, f=F):
    ke = elimination_rate(cl, v)
    return f * dose * ka / (v * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t))


# ---------- ODE solutions ----------
def iv_ode(t_eval, dose=DOSE, v=V, cl=CL):
    ke = elimination_rate(cl, v)

    def rhs(t, y):
        return [-ke * y[0]]            # y[0] = amount in central compartment (mg)

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), [dose],
                    t_eval=t_eval, rtol=1e-10, atol=1e-12)
    return sol.y[0] / v                # concentration (mg/L)


def oral_ode(t_eval, dose=DOSE, v=V, cl=CL, ka=KA, f=F):
    ke = elimination_rate(cl, v)

    def rhs(t, y):
        gut, central = y
        return [-ka * gut, ka * gut - ke * central]

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), [f * dose, 0.0],
                    t_eval=t_eval, rtol=1e-10, atol=1e-12)
    return sol.y[1] / v


if __name__ == "__main__":
    t = np.linspace(0, 48, 481)

    iv_err = np.max(np.abs(iv_ode(t) - iv_analytical(t)))
    oral_err = np.max(np.abs(oral_ode(t) - oral_analytical(t)))

    print(f"ke = {elimination_rate():.3f} 1/h")
    print(f"Half-life = {np.log(2) / elimination_rate():.2f} h")
    print(f"Max |ODE - analytical| IV   : {iv_err:.2e} mg/L")
    print(f"Max |ODE - analytical| oral : {oral_err:.2e} mg/L")