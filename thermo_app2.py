
from __future__ import annotations

import csv
import math
import tkinter as tk
from dataclasses import dataclass
from tkinter import filedialog, messagebox, ttk
from typing import Callable, List, Optional, Tuple


R_ENERGY = 8.314        
R_VOLUME = 83.14        
T_REF = 298.15          
P_REF = 1.0


@dataclass(frozen=True)
class Gas:
    name: str
    tc: float           
    pc: float           
    omega: float        
    a: float            
    b: float
    c: float
    d: float


GAS_LIST = [
    Gas("Methane",         190.6,  45.99,  0.012, 1.702,  9.081e-3,  -2.164e-6,  0.0),
    Gas("Ethane",          305.3,  48.72,  0.100, 1.131,  19.225e-3, -5.561e-6,  0.0),
    Gas("Propane",         369.8,  42.48,  0.152, 1.213,  28.785e-3, -8.824e-6,  0.0),
    Gas("n-Butane",        425.1,  37.96,  0.200, 1.935,  36.915e-3, -11.402e-6, 0.0),
    Gas("Ethylene",        282.3,  50.40,  0.087, 1.424,  14.394e-3, -4.392e-6,  0.0),
    Gas("Propylene",       365.6,  46.65,  0.140, 1.637,  22.706e-3, -6.915e-6,  0.0),
    Gas("Nitrogen",        126.2,  34.00,  0.038, 3.280,  0.593e-3,  0.0,        0.040e5),
    Gas("Oxygen",          154.6,  50.43,  0.022, 3.639,  0.506e-3,  0.0,       -0.227e5),
    Gas("Carbon dioxide",  304.2,  73.83,  0.224, 5.457,  1.045e-3,  0.0,       -1.157e5),
    Gas("Carbon monoxide", 132.9,  34.99,  0.048, 3.376,  0.557e-3,  0.0,       -0.031e5),
    Gas("Hydrogen",         33.19, 13.13, -0.216, 3.249,  0.422e-3,  0.0,        0.083e5),
    Gas("Ammonia",         405.7, 112.80,  0.253, 3.578,  3.020e-3,  0.0,       -0.186e5),
    Gas("Water",           647.1, 220.55,  0.345, 3.470,  1.450e-3,  0.0,        0.121e5),
    Gas("Argon",           150.9,  48.98, -0.004, 2.500,  0.0,       0.0,        0.0),
]


@dataclass
class State:
    """Result for one root of the cubic equation."""
    z: float
    v: float            
    h_res: float        
    s_res: float        
    ln_phi: float
    phi: float
    fugacity: float     



def solve_cubic(c2: float, c1: float, c0: float) -> List[float]:
    """Real roots of x^3 + c2*x^2 + c1*x + c0 = 0 (analytical formula)."""
    q = (c2 * c2 - 3.0 * c1) / 9.0
    r = (2.0 * c2 ** 3 - 9.0 * c2 * c1 + 27.0 * c0) / 54.0

    if r * r < q ** 3:                                 
        theta = math.acos(r / math.sqrt(q ** 3))
        m = -2.0 * math.sqrt(q)
        return [m * math.cos((theta + 2.0 * math.pi * k) / 3.0) - c2 / 3.0
                for k in range(3)]

    sign = 1.0 if r >= 0 else -1.0                     
    big_a = -sign * (abs(r) + math.sqrt(r * r - q ** 3)) ** (1.0 / 3.0)
    big_b = 0.0 if big_a == 0 else q / big_a
    return [big_a + big_b - c2 / 3.0]


class CubicEOS:
    """Generalized cubic EOS: P = RT/(V-b) - a/[(V+eps*b)(V+sigma*b)]."""

    NAMES = ("van der Waals", "Redlich-Kwong", "Soave-Redlich-Kwong", "Peng-Robinson")

    def __init__(self, kind: int):
        self.kind = kind
        self.name = self.NAMES[kind]
        self.sigma = {0: 0.0, 3: 1.0 + math.sqrt(2.0)}.get(kind, 1.0)
        self.eps = (1.0 - math.sqrt(2.0)) if kind == 3 else 0.0
        self.omega_b = {0: 1.0 / 8.0, 3: 0.07780}.get(kind, 0.08664)
        self.psi_a = {0: 27.0 / 64.0, 3: 0.45724}.get(kind, 0.42748)


    def _kappa(self, w: float) -> float:
        if self.kind == 2:
            return 0.480 + 1.574 * w - 0.176 * w * w
        return 0.37464 + 1.54226 * w - 0.26992 * w * w

    def alpha(self, tr: float, w: float) -> float:
        if self.kind == 0:
            return 1.0
        if self.kind == 1:
            return tr ** -0.5
        s = 1.0 + self._kappa(w) * (1.0 - math.sqrt(tr))
        return s * s

    def dln_alpha(self, tr: float, w: float) -> float:
        """d ln(alpha) / d ln(Tr)."""
        if self.kind == 0:
            return 0.0
        if self.kind == 1:
            return -0.5
        k = self._kappa(w)
        return -k * math.sqrt(tr) / (1.0 + k * (1.0 - math.sqrt(tr)))


    def beta_q(self, gas: Gas, t: float, p: float) -> Tuple[float, float]:
        tr, pr = t / gas.tc, p / gas.pc
        beta = self.omega_b * pr / tr
        q = self.psi_a * self.alpha(tr, gas.omega) / (self.omega_b * tr)
        return beta, q

    def roots(self, gas: Gas, t: float, p: float) -> List[float]:
        """Physical roots (Z > beta), smallest first."""
        beta, q = self.beta_q(gas, t, p)
        s = self.sigma + self.eps
        pr = self.sigma * self.eps
        c2 = s * beta - (1.0 + beta)
        c1 = pr * beta ** 2 - (1.0 + beta) * s * beta + q * beta
        c0 = -(1.0 + beta) * pr * beta ** 2 - q * beta ** 2
        return sorted(z for z in solve_cubic(c2, c1, c0) if z > beta + 1e-9)


    def state(self, gas: Gas, t: float, p: float, z: float) -> State:
        beta, q = self.beta_q(gas, t, p)
        if abs(self.sigma - self.eps) < 1e-12:                 
            i_int = beta / (z + self.eps * beta)
        else:
            i_int = math.log((z + self.sigma * beta) / (z + self.eps * beta)) / (self.sigma - self.eps)

        d = self.dln_alpha(t / gas.tc, gas.omega)
        ln_phi = z - 1.0 - math.log(z - beta) - q * i_int
        return State(
            z=z,
            v=z * R_VOLUME * t / p,
            h_res=R_ENERGY * t * (z - 1.0 + (d - 1.0) * q * i_int),
            s_res=R_ENERGY * (math.log(z - beta) + d * q * i_int),
            ln_phi=ln_phi,
            phi=math.exp(ln_phi),
            fugacity=math.exp(ln_phi) * p,
        )


    def saturation(self, gas: Gas, t: float) -> Optional[Tuple[float, State, State]]:
        """Return (Psat, vapor, liquid) or None. Equal fugacity iteration."""
        if t >= gas.tc:
            return None
        p = gas.pc * math.exp(5.373 * (1.0 + gas.omega) * (1.0 - gas.tc / t))   
        for _ in range(2000):
            zs = self.roots(gas, t, p)
            if not zs:
                return None
            if len(zs) < 3:                                     
                p *= 0.9 if zs[0] < 0.3 else 1.1
                continue
            vap, liq = self.state(gas, t, p, zs[2]), self.state(gas, t, p, zs[0])
            ratio = math.exp(liq.ln_phi - vap.ln_phi)
            if abs(ratio - 1.0) < 1e-10:
                return p, vap, liq
            p *= ratio
        return None

    def stable(self, gas: Gas, t: float, p: float) -> Tuple[State, str]:
        """Stable state: the root with the lower fugacity. Returns (state, phase)."""
        zs = self.roots(gas, t, p)
        if not zs:
            return self.state(gas, t, p, 1.0), "?"
        if len(zs) == 3:
            vap, liq = self.state(gas, t, p, zs[2]), self.state(gas, t, p, zs[0])
            return (liq, "Liquid") if liq.ln_phi < vap.ln_phi else (vap, "Vapor")

        only = self.state(gas, t, p, zs[0])
        if t >= gas.tc:
            return only, "Gas"
        sat = self.saturation(gas, t)
        if sat is not None:
            return only, "Liquid" if p > sat[0] else "Vapor"
        return only, "Liquid" if only.z < 0.3 else "Vapor"


def ideal_h(gas: Gas, t: float) -> float:
    """Ideal-gas enthalpy relative to T_REF, J/mol."""
    t0 = T_REF
    return R_ENERGY * (gas.a * (t - t0) + gas.b / 2.0 * (t ** 2 - t0 ** 2)
                       + gas.c / 3.0 * (t ** 3 - t0 ** 3) - gas.d * (1.0 / t - 1.0 / t0))


def ideal_s(gas: Gas, t: float, p: float) -> float:
    """Ideal-gas entropy relative to (T_REF, P_REF), J/(mol K)."""
    t0 = T_REF
    s = (gas.a * math.log(t / t0) + gas.b * (t - t0) + gas.c / 2.0 * (t ** 2 - t0 ** 2)
         - gas.d / 2.0 * (1.0 / t ** 2 - 1.0 / t0 ** 2))
    return R_ENERGY * s - R_ENERGY * math.log(p / P_REF)



Cards = List[Tuple[str, str]]
Result = Tuple[Cards, List[str], List[List[str]]]       

STATE_FIELDS = (                                        
    ("Z", 5, "z"),
    ("Molar volume V (cm3/mol)", 3, "v"),
    ("Residual enthalpy H^R (J/mol)", 2, "h_res"),
    ("Residual entropy S^R (J/(mol K))", 4, "s_res"),
    ("ln(phi)", 5, "ln_phi"),
    ("Fugacity coefficient phi", 5, "phi"),
    ("Fugacity f (bar)", 4, "fugacity"),
)


def state_rows(states: List[State]) -> List[List[str]]:
    return [[label] + [f"{getattr(s, attr):.{digits}f}" for s in states]
            for label, digits, attr in STATE_FIELDS]


@dataclass
class Mode:
    fields: List[Tuple[str, str, float, float]]          
    handler: Callable[[Gas, CubicEOS, List[float]], Result]



BG, SIDE, CARD = "#12161f", "#0b0e14", "#1b2230"
ACCENT, TEXT, MUTED, ERROR = "#2ec4b6", "#e8edf5", "#8a97ab", "#ff6b6b"


class ThermoApp:
    def __init__(self) -> None:
        self.root = tk.Tk()
        self.csv_rows: List[List[float]] = []
        self.csv_header = ["P_bar", "Phase", "Z", "V_cm3_mol", "HR_J_mol", "SR_J_molK", "phi", "f_bar"]
        self.entry_vars: List[tk.StringVar] = []
        self.mode_name = ""

        self.modes = {
            "Properties at T, P": Mode([("Temperature T (K)", "300", 1, 3000),
                                        ("Pressure P (bar)", "5", 1e-6, 5000)], self.mode_state),
            "Change in H and S": Mode([("T1 (K)", "300", 1, 3000), ("P1 (bar)", "1", 1e-6, 5000),
                                       ("T2 (K)", "400", 1, 3000), ("P2 (bar)", "20", 1e-6, 5000)],
                                      self.mode_change),
            "Table vs pressure": Mode([("T (K)", "400", 1, 3000), ("P start (bar)", "0.5", 1e-6, 5000),
                                       ("P end (bar)", "40", 1e-6, 5000), ("Points (2-60)", "15", 2, 60)],
                                      self.mode_table),
            "Compare the 4 EOS": Mode([("T (K)", "350", 1, 3000), ("P (bar)", "50", 1e-6, 5000)],
                                      self.mode_compare),
            "Saturation (Psat)": Mode([("T (K), below Tc", "300", 1, 3000)], self.mode_saturation),
        }
        self.build_ui()
        self.select_mode("Properties at T, P")

    
    def gas(self) -> Gas:
        return GAS_LIST[self.gas_box.current()]

    def eos(self) -> CubicEOS:
        return CubicEOS(self.eos_box.current())

    
    def mode_state(self, gas: Gas, eos: CubicEOS, v: List[float]) -> Result:
        t, p = v
        zs = eos.roots(gas, t, p)
        best, phase = eos.stable(gas, t, p)
        if len(zs) == 3:
            states = [eos.state(gas, t, p, zs[2]), eos.state(gas, t, p, zs[0])]
            columns = ["Property", "Vapor-like root", "Liquid-like root"]
        else:
            states = [best]
            columns = ["Property", "Value"]

        rows = state_rows(states)
        rows.append(["Total enthalpy H (J/mol)"] + [f"{ideal_h(gas, t) + s.h_res:.2f}" for s in states])
        rows.append(["Total entropy S (J/(mol K))"] + [f"{ideal_s(gas, t, p) + s.s_res:.4f}" for s in states])
        cards = [("Phase", phase), ("Z", f"{best.z:.4f}"), ("φ", f"{best.phi:.4f}"),
                 ("f (bar)", f"{best.fugacity:.4f}"),
                 ("H (J/mol)", f"{ideal_h(gas, t) + best.h_res:.1f}"),
                 ("S (J/mol K)", f"{ideal_s(gas, t, p) + best.s_res:.3f}")]
        return cards, columns, rows

    def mode_change(self, gas: Gas, eos: CubicEOS, v: List[float]) -> Result:
        t1, p1, t2, p2 = v
        a, phase_a = eos.stable(gas, t1, p1)
        b, phase_b = eos.stable(gas, t2, p2)
        dh_ideal = ideal_h(gas, t2) - ideal_h(gas, t1)
        ds_ideal = ideal_s(gas, t2, p2) - ideal_s(gas, t1, p1)
        dh = dh_ideal + b.h_res - a.h_res
        ds = ds_ideal + b.s_res - a.s_res

        rows = [["Phase", phase_a, phase_b]] + state_rows([a, b])
        cards = [("ΔH real (J/mol)", f"{dh:.2f}"), ("ΔS real (J/mol K)", f"{ds:.4f}"),
                 ("ΔH ideal (J/mol)", f"{dh_ideal:.2f}"), ("ΔS ideal (J/mol K)", f"{ds_ideal:.4f}")]
        return cards, ["Property", "State 1", "State 2"], rows

    def mode_table(self, gas: Gas, eos: CubicEOS, v: List[float]) -> Result:
        t, p_start, p_end, n = v[0], v[1], v[2], int(v[3])
        if p_end < p_start:
            raise ValueError("P end must be larger than P start.")
        self.csv_rows = []
        rows = []
        for i in range(n):
            p = p_start + (p_end - p_start) * i / (n - 1)
            s, phase = eos.stable(gas, t, p)
            self.csv_rows.append([p, phase, s.z, s.v, s.h_res, s.s_res, s.phi, s.fugacity])
            rows.append([f"{p:.3f}", phase, f"{s.z:.4f}", f"{s.v:.2f}", f"{s.h_res:.1f}",
                         f"{s.s_res:.3f}", f"{s.phi:.4f}", f"{s.fugacity:.3f}"])
        columns = ["P (bar)", "Phase", "Z", "V (cm3/mol)", "H^R (J/mol)", "S^R (J/mol K)", "phi", "f (bar)"]
        return [("T (K)", f"{t:.2f}"), ("Points", str(n))], columns, rows

    def mode_compare(self, gas: Gas, eos: CubicEOS, v: List[float]) -> Result:
        t, p = v
        rows = []
        for kind in range(4):
            model = CubicEOS(kind)
            s, phase = model.stable(gas, t, p)
            rows.append([model.name, phase, f"{s.z:.4f}", f"{s.v:.2f}", f"{s.h_res:.1f}",
                         f"{s.s_res:.3f}", f"{s.phi:.4f}"])
        columns = ["EOS", "Phase", "Z", "V (cm3/mol)", "H^R (J/mol)", "S^R (J/mol K)", "phi"]
        return [("T (K)", f"{t:.2f}"), ("P (bar)", f"{p:.2f}")], columns, rows

    def mode_saturation(self, gas: Gas, eos: CubicEOS, v: List[float]) -> Result:
        found = eos.saturation(gas, v[0])
        if found is None:
            raise ValueError(f"No vapor-liquid equilibrium: T must be below Tc = {gas.tc} K and not too low.")
        p_sat, vap, liq = found
        cards = [("Psat (bar)", f"{p_sat:.4f}"), ("ΔH vap (J/mol)", f"{vap.h_res - liq.h_res:.1f}"),
                 ("f vapor (bar)", f"{vap.fugacity:.4f}"), ("f liquid (bar)", f"{liq.fugacity:.4f}")]
        return cards, ["Property", "Saturated vapor", "Saturated liquid"], state_rows([vap, liq])

    
    def build_ui(self) -> None:
        root = self.root
        root.title("Gas Property Lab")
        root.geometry("1020x680")
        root.minsize(860, 560)
        root.configure(bg=BG)

        style = ttk.Style()
        style.theme_use("clam")
        style.configure("Treeview", background=CARD, fieldbackground=CARD, foreground=TEXT,
                        rowheight=26, borderwidth=0)
        style.configure("Treeview.Heading", background=SIDE, foreground=ACCENT,
                        font=("Segoe UI", 10, "bold"))
        style.map("Treeview", background=[("selected", "#2c3a52")])

        # sidebar
        side = tk.Frame(root, bg=SIDE, width=240)
        side.pack(side="left", fill="y")
        side.pack_propagate(False)
        tk.Label(side, text="GAS PROPERTY LAB", bg=SIDE, fg=ACCENT,
                 font=("Segoe UI", 15, "bold")).pack(pady=(24, 2))
        tk.Label(side, text="Cubic EOS · H · S · Fugacity", bg=SIDE, fg=MUTED,
                 font=("Segoe UI", 9)).pack(pady=(0, 18))

        tk.Label(side, text="GAS", bg=SIDE, fg=MUTED, font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=18)
        self.gas_box = ttk.Combobox(side, values=[g.name for g in GAS_LIST], state="readonly")
        self.gas_box.current(2)
        self.gas_box.pack(fill="x", padx=18, pady=(2, 10))

        tk.Label(side, text="EQUATION OF STATE", bg=SIDE, fg=MUTED,
                 font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=18)
        self.eos_box = ttk.Combobox(side, values=list(CubicEOS.NAMES), state="readonly")
        self.eos_box.current(3)
        self.eos_box.pack(fill="x", padx=18, pady=(2, 18))

        self.nav_buttons = {}
        for name in self.modes:
            button = tk.Button(side, text=name, anchor="w", relief="flat", bg=SIDE, fg=TEXT,
                               activebackground=CARD, activeforeground=ACCENT, padx=18, pady=9,
                               font=("Segoe UI", 10), command=lambda n=name: self.select_mode(n))
            button.pack(fill="x")
            self.nav_buttons[name] = button

        # main area
        main = tk.Frame(root, bg=BG)
        main.pack(side="left", fill="both", expand=True, padx=18, pady=14)
        self.title_label = tk.Label(main, text="", bg=BG, fg=TEXT, font=("Segoe UI", 16, "bold"), anchor="w")
        self.title_label.pack(fill="x")
        self.input_frame = tk.Frame(main, bg=BG)
        self.input_frame.pack(fill="x", pady=8)
        self.status = tk.Label(main, text="", bg=BG, fg=MUTED, anchor="w", font=("Segoe UI", 9))
        self.status.pack(fill="x")
        self.card_frame = tk.Frame(main, bg=BG)
        self.card_frame.pack(fill="x", pady=6)

        table_frame = tk.Frame(main, bg=BG)
        table_frame.pack(fill="both", expand=True)
        self.tree = ttk.Treeview(table_frame, show="headings")
        scroll_y = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll_y.set)
        scroll_y.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)

        root.bind("<Return>", lambda event: self.run())

    def select_mode(self, name: str) -> None:
        self.mode_name = name
        for key, button in self.nav_buttons.items():
            button.config(bg=CARD if key == name else SIDE, fg=ACCENT if key == name else TEXT)
        self.title_label.config(text=name)
        self.status.config(text="Enter the values and press Run (or ENTER).", fg=MUTED)
        self.clear_results()

        for widget in self.input_frame.winfo_children():
            widget.destroy()
        self.entry_vars = []
        for col, (label, default, _lo, _hi) in enumerate(self.modes[name].fields):
            tk.Label(self.input_frame, text=label, bg=BG, fg=MUTED,
                     font=("Segoe UI", 9)).grid(row=0, column=col, sticky="w", padx=6)
            var = tk.StringVar(value=default)
            tk.Entry(self.input_frame, textvariable=var, width=13, bg=CARD, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 11)).grid(row=1, column=col, padx=6, pady=3, ipady=4)
            self.entry_vars.append(var)

        last = len(self.modes[name].fields)
        tk.Button(self.input_frame, text="Run", bg=ACCENT, fg="#06201d", relief="flat", padx=22,
                  font=("Segoe UI", 10, "bold"), command=self.run).grid(row=1, column=last, padx=12)
        if name == "Table vs pressure":
            tk.Button(self.input_frame, text="Save CSV", bg=CARD, fg=TEXT, relief="flat", padx=14,
                      font=("Segoe UI", 10), command=self.save_csv).grid(row=1, column=last + 1)

    # ----- running a calculation -----
    def read_values(self) -> List[float]:
        values = []
        for (label, _default, lo, hi), var in zip(self.modes[self.mode_name].fields, self.entry_vars):
            try:
                x = float(var.get())
            except ValueError:
                raise ValueError(f"{label} must be a number.")
            if not lo <= x <= hi:
                raise ValueError(f"{label} must be between {lo:g} and {hi:g}.")
            values.append(x)
        return values

    def run(self) -> None:
        try:
            values = self.read_values()
            result = self.modes[self.mode_name].handler(self.gas(), self.eos(), values)
        except ValueError as err:
            self.status.config(text=f"⚠  {err}", fg=ERROR)
            return
        except ArithmeticError:
            self.status.config(text="⚠  The calculation failed for these values.", fg=ERROR)
            return
        self.status.config(text=f"{self.gas().name}  ·  {self.eos().name}", fg=MUTED)
        self.show_result(*result)

    def clear_results(self) -> None:
        for widget in self.card_frame.winfo_children():
            widget.destroy()
        self.tree.delete(*self.tree.get_children())

    def show_result(self, cards: Cards, columns: List[str], rows: List[List[str]]) -> None:
        self.clear_results()
        for label, value in cards:
            card = tk.Frame(self.card_frame, bg=CARD, padx=12, pady=6)
            card.pack(side="left", padx=(0, 8))
            tk.Label(card, text=label, bg=CARD, fg=MUTED, font=("Segoe UI", 8)).pack(anchor="w")
            tk.Label(card, text=value, bg=CARD, fg=ACCENT, font=("Segoe UI", 13, "bold")).pack(anchor="w")

        self.tree.configure(columns=columns)
        for index, column in enumerate(columns):
            self.tree.heading(column, text=column)
            self.tree.column(column, width=150 if index else 230, anchor="e" if index else "w", stretch=True)
        for row in rows:
            self.tree.insert("", "end", values=row)

    def save_csv(self) -> None:
        if not self.csv_rows:
            messagebox.showinfo("Nothing to save", "Press Run first to make the table.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="results.csv",
                                            filetypes=[("CSV file", "*.csv")])
        if not path:
            return
        with open(path, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(self.csv_header)
            writer.writerows(self.csv_rows)
        messagebox.showinfo("Saved", f"Table saved in {path}")

    def run_app(self) -> None:
        self.root.mainloop()


if __name__ == "__main__":
    ThermoApp().run_app()
