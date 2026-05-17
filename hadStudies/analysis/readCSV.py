import pandas as pd
import sys
import argparse
 
# ---------------------------------------------------------------------------
# Particle catalogue (ONLY physical projectiles used in input CSV)
# ---------------------------------------------------------------------------
PARTICLE_INFO = {
     211: {"probe_type": "pion",    "label": r"$\pi^+$", "cpp_name": "piplus"},
    -211: {"probe_type": "pion",    "label": r"$\pi^-$", "cpp_name": "piminus"},
     111: {"probe_type": "pion",    "label": r"$\pi^0$", "cpp_name": "pizero"},
    2212: {"probe_type": "nucleon", "label": r"$p$",     "cpp_name": "proton"},
    2112: {"probe_type": "nucleon", "label": r"$n$",     "cpp_name": "neutron"},
}
 
NUCLEUS_NAMES = {
    1000060120: ("C12",  r"$^{12}$C",  "c12"),
    1000080160: ("O16",  r"$^{16}$O",  "o16"),
    1000180400: ("Ar40", r"$^{40}$Ar", "ar40"),
    1000260560: ("Fe56", r"$^{56}$Fe", "fe56"),
}
 
MODEL_NAMES = ["hA2018", "hN2018", "INCL", "Geant4"]
MODEL_SHORT = {"hA2018": "hA", "hN2018": "hN", "INCL": "INCL", "Geant4": "G4"}
 
# ---------------------------------------------------------------------------
# Slots
# ---------------------------------------------------------------------------
PION_SLOTS = [
    ("Difference", "Mean", "DifMean", "Diff.", r"$\mu$"),
    ("Difference", "Std",  "DifStd",  "Diff.", r"$\sigma$"),
    ("Sum",        "Mean", "SumMean", "Sum",   r"$\mu$"),
    ("Sum",        "Std",  "SumStd",  "Sum",   r"$\sigma$"),
]
 
NUCLEON_LINEAR_SLOTS = [
    ("Difference", "Mean", "DifMean", "Diff.", r"$\mu$"),
    ("Difference", "Std",  "DifStd",  "Diff.", r"$\sigma$"),
]
 
# ---------------------------------------------------------------------------
# Structs
# ---------------------------------------------------------------------------
SHARED_STRUCTS = """\
#include <map>
#include <variant>
 
struct LinearParams {
    double slope, intercept;
};
 
struct GammaParams {
    double slope, intercept, floor, exp;
};
 
struct ModelParamsLinear {
    LinearParams hA, hN, INCL, G4;
};
 
struct ModelParamsGamma {
    GammaParams hA, hN, INCL, G4;
};
 
enum class ProbeType { Pion, Nucleon };
 
struct NucleonTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsGamma  SumGamma;
};
 
struct PionTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsLinear SumMean;
    ModelParamsLinear SumStd;
};
 
struct ParticleTargetParams {
    ProbeType probeType;
    std::variant<PionTargetParams, NucleonTargetParams> params;
};
 
// Keyed by [PDG code][target PDG code]
std::map<int, std::map<int, ParticleTargetParams>> fitParams;
"""
 
# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------
def _load_csv(path):
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [str(c).strip() for c in df.columns]
 
    df["particle"] = pd.to_numeric(df["particle"], errors="coerce")
    df["target"]   = pd.to_numeric(df["target"],   errors="coerce")
 
    df = df.dropna(subset=["particle", "target"])
 
    df["particle"] = df["particle"].astype(int)
    df["target"]   = df["target"].astype(int)
 
    return df
 
 
def _fmt(val, decimals=4):
    try:
        f = float(val)
        if pd.isna(f):
            return "---"
        return f"{f:.{decimals}f}"
    except Exception:
        return "---"
 
 
def _get(df, target, pdg, category, quantity):
    return df[
        (df["target"]   == target) &
        (df["particle"] == pdg) &
        (df["category"] == category) &
        (df["quantity"] == quantity)
    ]
 
 
def _particles_in(df):
    particles = df["particle"].dropna().astype(int).unique()
    return sorted([p for p in particles if p in PARTICLE_INFO])
 
 
def particle_info(pdg):
    return PARTICLE_INFO[pdg]
 
 
# ---------------------------------------------------------------------------
# C++ helpers
# ---------------------------------------------------------------------------
def _linear_block(df, target, pdg, cat, qty, indent=12):
    pad = " " * indent
    row = _get(df, target, pdg, cat, qty)
    lines = [pad + "{"]
    for i, model in enumerate(MODEL_NAMES):
        ms    = MODEL_SHORT[model]
        comma = "," if i != len(MODEL_NAMES) - 1 else ""
        if row.empty:
            lines.append(f"{pad}    {{0.0, 0.0}}{comma}  // {ms} missing")
        else:
            s = row.iloc[0][f"{model}_slope"]
            b = row.iloc[0][f"{model}_intercept"]
            lines.append(f"{pad}    {{{s:.4f}, {b:.4f}}}{comma}  // {ms}")
    lines.append(pad + "}")
    return lines
 
 
def _gamma_block(df, target, pdg, indent=12):
    pad = " " * indent
    row = _get(df, target, pdg, "Sum", "Gamma")
    lines = [pad + "{"]
    for i, model in enumerate(MODEL_NAMES):
        ms    = MODEL_SHORT[model]
        comma = "," if i != len(MODEL_NAMES) - 1 else ""
        if row.empty:
            lines.append(f"{pad}    {{0.0, 0.0, 0.0, 1.0}}{comma}  // {ms} missing")
            continue
        s  = row.iloc[0][f"{model}_slope"]
        b  = row.iloc[0][f"{model}_intercept"]
        fl = row.iloc[0][f"{model}_floor"] if f"{model}_floor" in row.columns else 0.0
        ex = row.iloc[0][f"{model}_exp"]   if f"{model}_exp"   in row.columns else 1.0
        lines.append(f"{pad}    {{{s:.4f}, {b:.4f}, {fl:.4f}, {ex:.4f}}}{comma}  // {ms}")
    lines.append(pad + "}")
    return lines
 
 
# ---------------------------------------------------------------------------
# Unified C++ output
# ---------------------------------------------------------------------------
def csv_to_cpp(df):
    pdgs = _particles_in(df)
    pion_pdgs    = [p for p in pdgs if PARTICLE_INFO[p]["probe_type"] == "pion"]
    nucleon_pdgs = [p for p in pdgs if PARTICLE_INFO[p]["probe_type"] == "nucleon"]
 
    lines = ["void makeAllParams() {"]
 
    for pdg in nucleon_pdgs + pion_pdgs:
        info       = PARTICLE_INFO[pdg]
        probe_type = info["probe_type"]
        cpp_name   = info["cpp_name"]
        is_pion    = (probe_type == "pion")
 
        variant_type  = "PionTargetParams"    if is_pion else "NucleonTargetParams"
        probe_enum    = "ProbeType::Pion"     if is_pion else "ProbeType::Nucleon"
        slots         = PION_SLOTS            if is_pion else NUCLEON_LINEAR_SLOTS
 
        lines.append(f"")
        lines.append(f"    // PDG {pdg} ({cpp_name})")
 
        for target, (nuc, _, _) in NUCLEUS_NAMES.items():
            lines.append(f"    fitParams[{pdg}][{target}] = {{  // {nuc}")
            lines.append(f"        {probe_enum},")
            lines.append(f"        {variant_type}{{")
 
            last_slot_idx = len(slots) - 1
            for slot_idx, (cat, qty, field, _, _) in enumerate(slots):
                slot_comma = "," if (slot_idx != last_slot_idx or not is_pion) else ""
                # For nucleons there's always the gamma block after linear slots
                if not is_pion and slot_idx == last_slot_idx:
                    slot_comma = ","
 
                lines.append(f"            // {field}")
                block = _linear_block(df, target, pdg, cat, qty, indent=12)
                block[-1] += slot_comma
                lines.extend(block)
 
            if not is_pion:
                lines.append(f"            // SumGamma")
                block = _gamma_block(df, target, pdg, indent=12)
                lines.extend(block)
 
            lines.append(f"        }}")
            lines.append(f"    }};")
 
    lines.append("}")
    return "\n".join(lines)
 
 
# ---------------------------------------------------------------------------
# LaTeX output
# ---------------------------------------------------------------------------
def csv_to_latex(df, pdg, decimals=4):
    info       = PARTICLE_INFO[pdg]
    probe_type = info["probe_type"]
    label      = info["label"]
    cpp_name   = info["cpp_name"]

    slots        = PION_SLOTS if probe_type == "pion" else NUCLEON_LINEAR_SLOTS
    lin_col_spec = "|c|" + "c|" * (2 * len(MODEL_NAMES))

    lines = []

    for target, (_, nuc_tex, nuc_id) in NUCLEUS_NAMES.items():
        L = []

        L.append(r"\begin{table}[h!]")
        L.append(rf"\caption{{Fit parameters for {label} on {nuc_tex}}}")
        L.append(r"\centering")
        L.append(rf"\begin{{tabular}}{{{lin_col_spec}}}")
        L.append(r"\hline")

        header    = [""] + [rf"\multicolumn{{2}}{{c|}}{{{m}}}" for m in MODEL_NAMES]
        subheader = [""] + [x for _ in MODEL_NAMES for x in [r"$Slope$", r"$y-intercept$"]]

        L.append(" & ".join(header) + r" \\")
        L.append(r"\hline")
        L.append(" & ".join(subheader) + r" \\")
        L.append(r"\hline")

        for cat, qty, _, row_label, qty_label in slots:
            row   = _get(df, target, pdg, cat, qty)
            cells = [f"{row_label} {qty_label}"]

            for model in MODEL_NAMES:
                if row.empty:
                    cells += ["---", "---"]
                else:
                    cells += [
                        _fmt(row.iloc[0][f"{model}_slope"],     decimals),
                        _fmt(row.iloc[0][f"{model}_intercept"], decimals),
                    ]

            L.append(" & ".join(cells) + r" \\")
            L.append(r"\hline")

        L.append(r"\end{tabular}")
        L.append(rf"\label{{tab:{cpp_name}_{nuc_id}}}")
        L.append(r"\end{table}")

        lines.append("\n".join(L))

        # Second table: nucleon Sum/Gamma (4 parameters)
        if probe_type == "nucleon":
            G = []
            gamma_col_spec = "|c|" + "c|" * (4 * len(MODEL_NAMES))

            G.append(r"\begin{table}[h!]")
            G.append(rf"\caption{{Sum $\Gamma$ fit parameters for {label} on {nuc_tex}}}")
            G.append(r"\centering")
            G.append(rf"\begin{{tabular}}{{{gamma_col_spec}}}")
            G.append(r"\hline")

            gamma_header    = [""] + [rf"\multicolumn{{4}}{{c|}}{{{m}}}" for m in MODEL_NAMES]
            gamma_subheader = [""] + [x for _ in MODEL_NAMES for x in [r"$Slope$", r"$Initial$", r"$Floor$", r"$Power$"]]

            G.append(" & ".join(gamma_header) + r" \\")
            G.append(r"\hline")
            G.append(" & ".join(gamma_subheader) + r" \\")
            G.append(r"\hline")

            row   = _get(df, target, pdg, "Sum", "Gamma")
            cells = ["Sum $\\Gamma$"]

            for model in MODEL_NAMES:
                if row.empty:
                    cells += ["---", "---", "---", "---"]
                else:
                    s  = _fmt(row.iloc[0][f"{model}_slope"],     decimals)
                    b  = _fmt(row.iloc[0][f"{model}_intercept"], decimals)
                    fl = _fmt(row.iloc[0].get(f"{model}_floor", 0.0), decimals)
                    ex = _fmt(row.iloc[0].get(f"{model}_exp",   1.0), decimals)
                    cells += [s, b, fl, ex]

            G.append(" & ".join(cells) + r" \\")
            G.append(r"\hline")
            G.append(r"\end{tabular}")
            G.append(rf"\label{{tab:{cpp_name}_{nuc_id}_gamma}}")
            G.append(r"\end{table}")

            lines.append("\n".join(G))

    return "\n\n".join(lines)
 
# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csvfiles", nargs="+")
    parser.add_argument("--latex",    action="store_true")
    parser.add_argument("--decimals", type=int, default=3)
    args = parser.parse_args()
 
    df = pd.concat([_load_csv(p) for p in args.csvfiles], ignore_index=True)
 
    pdgs         = _particles_in(df)
    pion_pdgs    = [p for p in pdgs if PARTICLE_INFO[p]["probe_type"] == "pion"]
    nucleon_pdgs = [p for p in pdgs if PARTICLE_INFO[p]["probe_type"] == "nucleon"]
 
    if not args.latex:
        print(SHARED_STRUCTS)
        print(csv_to_cpp(df))
 
    if args.latex:
        print("% LaTeX output")
        for p in nucleon_pdgs + pion_pdgs:
            print(csv_to_latex(df, p, decimals=args.decimals))
 