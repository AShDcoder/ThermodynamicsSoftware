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
  1. **Properties at $T, P$:** Evaluates $Z, V, H^R, S^R, \ln\phi, \phi, f$, and total $H$ and $S$ (including ideal-gas contributions relative to $298.15\text{ K}, 1\text{ bar}$). Displays both vapor-like and liquid-like roots when 3 real roots exist.
  2. **Change in $H$ and $S$:** Computes net real-gas and ideal-gas enthalpy ($\Delta H$) and entropy ($\Delta S$) changes between State 1 $(T_1, P_1)$ and State 2 $(T_2, P_2)$.
  3. **Table vs Pressure:** Generates parametric pressure sweeps over $N$ points ($2 \le N \le 60$) at constant temperature with one-click **CSV Export**.
  4. **Compare the 4 EOS:** Side-by-side benchmark comparison of vdW, RK, SRK, and PR predictions at a given $(T, P)$.
  5. **Saturation ($P^{\text{sat}}$):** Calculates vapor-liquid equilibrium saturation pressure ($P^{\text{sat}}$) and enthalpy of vaporization ($\Delta H_{\text{vap}}$) at subcritical temperatures ($T < T_c$) using Wilson's initial guess and equal-fugacity successive substitution.

---

## Theoretical Framework

### 1. Generalized Cubic Equation of State
Following the unified formulation by Smith, Van Ness, and Abbott, all four cubic equations of state are expressed as:

$$P = \frac{RT}{V - b} - \frac{a(T)}{(V + \epsilon b)(V + \sigma b)}$$

where the pure-component parameters are defined by critical constants $(T_c, P_c)$ and universal model constants $(\sigma, \epsilon, \Omega, \Psi)$:

$$b = \Omega \frac{RT_c}{P_c}, \quad a(T) = \Psi \frac{\alpha(T_r) R^2 T_c^2}{P_c}$$

| Equation of State | $\sigma$ | $\epsilon$ | $\Omega$ | $\Psi$ | $\alpha(T_r)$ and $\kappa(\omega)$ |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **van der Waals (vdW)** | $0$ | $0$ | $1/8 \ (0.1250)$ | $27/64 \ (0.42188)$ | $\alpha = 1$ |
| **Redlich-Kwong (RK)** | $1$ | $0$ | $0.08664$ | $0.42748$ | $\alpha = T_r^{-0.5}$ |
| **Soave-Redlich-Kwong (SRK)** | $1$ | $0$ | $0.08664$ | $0.42748$ | $\alpha = \left[1 + \kappa(1 - \sqrt{T_r})\right]^2$<br>$\kappa = 0.480 + 1.574\omega - 0.176\omega^2$ |
| **Peng-Robinson (PR)** | $1 + \sqrt{2}$ | $1 - \sqrt{2}$ | $0.07780$ | $0.45724$ | $\alpha = \left[1 + \kappa(1 - \sqrt{T_r})\right]^2$<br>$\kappa = 0.37464 + 1.54226\omega - 0.26992\omega^2$ |

Introducing the dimensionless parameters $\beta = \frac{bP}{RT}$ and $q = \frac{a(T)}{bRT}$ transforms the EOS into a cubic polynomial in the compressibility factor $Z = \frac{PV}{RT}$:

$$Z^3 + \left[(\sigma + \epsilon)\beta - (1 + \beta)\right]Z^2 + \left[\sigma\epsilon\beta^2 - (1 + \beta)(\sigma + \epsilon)\beta + q\beta\right]Z - \left[(1 + \beta)\sigma\epsilon\beta^2 + q\beta^2\right] = 0$$

### 2. Residual Departure Functions & Fugacity
Residual enthalpy ($H^R$), residual entropy ($S^R$), and the fugacity coefficient ($\phi$) are calculated analytically via the departure integral $I$:

$$I = \begin{cases} \frac{1}{\sigma - \epsilon} \ln\left(\frac{Z + \sigma\beta}{Z + \epsilon\beta}\right) & \text{for RK, SRK, PR } (\sigma \neq \epsilon) \\ \frac{\beta}{Z} & \text{for vdW } (\sigma = \epsilon = 0) \end{cases}$$

$$\frac{H^R}{RT} = Z - 1 + \left(\frac{d\ln\alpha}{d\ln T_r} - 1\right)qI$$

$$\frac{S^R}{R} = \ln(Z - \beta) + \left(\frac{d\ln\alpha}{d\ln T_r}\right)qI$$

$$\ln\phi = Z - 1 - \ln(Z - \beta) - qI \implies f = \phi P$$

### 3. Ideal-Gas Contributions & State Transitions
Total enthalpy ($H = H^{ig} + H^R$) and entropy ($S = S^{ig} + S^R$) are evaluated relative to the standard reference state ($T_0 = 298.15\text{ K}, P_0 = 1.0\text{ bar}$, where $H_0^{ig} = 0, S_0^{ig} = 0$) using temperature-dependent heat capacity polynomials:

$$\frac{C_p^{ig}}{R} = A + BT + CT^2 + DT^{-2}$$

$$H^{ig} = R\left[A(T - T_0) + \frac{B}{2}(T^2 - T_0^2) + \frac{C}{3}(T^3 - T_0^3) - D\left(\frac{1}{T} - \frac{1}{T_0}\right)\right]$$

$$S^{ig} = R\left[A\ln\left(\frac{T}{T_0}\right) + B(T - T_0) + \frac{C}{2}(T^2 - T_0^2) - \frac{D}{2}\left(\frac{1}{T^2} - \frac{1}{T_0^2}\right)\right] - R\ln\left(\frac{P}{P_0}\right)$$

---

## Requirements & Installation

The software relies solely on the **Python 3 Standard Library** (`math`, `csv`, `dataclasses`, `typing`, and `tkinter`)—no external `pip` packages are required.

### 1. Prerequisites
* **Python 3.8 or higher**
* **`tkinter` GUI library** (bundled by default with standard Windows and macOS Python installers).
  * *Linux (Ubuntu/Debian) users:* If `tkinter` is not installed, run:
    ```bash
    sudo apt-get update && sudo apt-get install python3-tk
    ```

### 2. Running the Application
1. Clone this repository or download `thermo_app2.py` to your local directory.
2. Open a terminal or command prompt in the project folder and run:
   ```bash
   python thermo_app2.py
