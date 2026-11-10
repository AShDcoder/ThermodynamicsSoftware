# Gas Property Lab: Thermodynamic Properties of Pure Gases using Cubic Equations of State

**Gas Property Lab** is an interactive desktop GUI application developed in Python (`tkinter` and `ttk`) designed to evaluate fundamental real-gas residual thermodynamic properties, total enthalpy and entropy departures, pure-component fugacity, and vapor-liquid saturation pressures ($P^{\text{sat}}$) across four classical cubic equations of state (EOS).

---

## Key Features

* **Four Unified Cubic Equations of State:** Supports van der Waals (vdW), Redlich-Kwong (RK), Soave-Redlich-Kwong (SRK), and Peng-Robinson (PR) formulations.
* **Built-in Chemical Database:** Pre-configured critical constants ($T_c, P_c$), acentric factors ($\omega$), and ideal-gas heat capacity polynomial coefficients ($A, B, C, D$) for **14 pure gases**:
  * *Hydrocarbons:* Methane, Ethane, Propane, n-Butane, Ethylene, Propylene
  * *Inorganics & Industrials:* Nitrogen, Oxygen, Carbon dioxide ($\text{CO}_2$), Carbon monoxide ($\text{CO}$), Hydrogen, Ammonia, Water, Argon
* **Analytical Cubic Root Solver:** Computes exact real compressibility roots ($Z$) using trigonometric and Cardano methods, automatically filtering for physical domain limits ($Z > \beta$).
* **Automated Phase Stability Analysis:** Evaluates Gibbs energy departures via fugacity coefficients ($\phi$) to isolate thermodynamically stable vapor and liquid roots when inside the two-phase envelope.
* **Five Interactive Calculation Modes:**
  1. **Properties at
