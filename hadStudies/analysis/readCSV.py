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
MODEL_SHORT = {"hA2018": "hA2018", "hN2018": "hN2018", "INCL": "INCL++", "Geant4": "Geant4"}

# ---------------------------------------------------------------------------
# Slots
# Tuple: (category, quantity, field_key, display_label, param_type)
# param_type: "linear" | "gamma"
# ---------------------------------------------------------------------------
PION_SLOTS = [
    ("Difference", "Mean",  "DifMean", r"Diff.\ $\mu$",  "linear"),
    ("Difference", "Std",   "DifStd",  r"Diff.\ $\sigma$","linear"),
    ("Sum",        "Mean",  "SumMean", r"Sum $\mu$",      "linear"),
    ("Sum",        "Std",   "SumStd",  r"Sum $\sigma$",   "linear"),
]

NUCLEON_SLOTS = [
    ("Difference", "Mean",  "DifMean",  r"Diff.\ $\mu$",  "linear"),
    ("Difference", "Std",   "DifStd",   r"Diff.\ $\sigma$","linear"),
    ("Sum",        "Gamma", "SumGamma", r"Sum $\Gamma$",   "gamma"),
]

# Columns per model for each param type
NCOLS = {"linear": 2, "gamma": 4}

# ---------------------------------------------------------------------------
# Structs (unchanged)
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


# ---------------------------------------------------------------------------
# C++ helpers (unchanged)
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
# Unified C++ output (unchanged)
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
        slots         = PION_SLOTS            if is_pion else NUCLEON_SLOTS

        lines.append(f"")
        lines.append(f"    // PDG {pdg} ({cpp_name})")

        for target, (nuc, _, _) in NUCLEUS_NAMES.items():
            lines.append(f"    fitParams[{pdg}][{target}] = {{  // {nuc}")
            lines.append(f"        {probe_enum},")
            lines.append(f"        {variant_type}{{")

            linear_slots = [(c, q, f, l, t) for c, q, f, l, t in slots if t == "linear"]
            for slot_idx, (cat, qty, field, _, _) in enumerate(linear_slots):
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
# LaTeX: one table per nucleus, all particles as row groups
# ---------------------------------------------------------------------------
# Column layout per model (widest case = gamma: slope, intercept, floor, power).
# Linear rows use only the first two sub-columns; the remaining two are left blank.
# This keeps a single tabular environment per table.
#
# Column spec:  | particle-label | fit-type | [slope | intercept | floor | power] x N_models |
# ---------------------------------------------------------------------------

def _subheader_row(n_models):
    """Sub-header: one group of 4 per model."""
    cells = [r"\multicolumn{2}{c|}{}", ""]   # particle col + fit-type col
    for _ in range(n_models):
        cells += [r"$m$", r"$b$", r"$c$", r"$\nu$"]
    return " & ".join(cells) + r" \\"


def _model_header_row(n_models):
    """Top header row: model names spanning 4 sub-columns each."""
    cells = [r"\multicolumn{2}{c|}{}", ""]   # particle col + fit-type col (blank)
    for model in MODEL_NAMES[:n_models]:
        ms = MODEL_SHORT[model]
        cells.append(rf"\multicolumn{{4}}{{c|}}{{{ms}}}")
    return " & ".join(cells) + r" \\"


def _data_row(df, target, pdg, cat, qty, param_type, row_label, decimals):
    """Build a single data row. Linear rows leave floor/power blank."""
    cells = [row_label]  # fit-type cell; particle cell added by caller via multirow

    if param_type == "linear":
        row = _get(df, target, pdg, cat, qty)
        for model in MODEL_NAMES:
            if row.empty:
                cells += ["---", "---", "", ""]
            else:
                s = _fmt(row.iloc[0][f"{model}_slope"],     decimals)
                b = _fmt(row.iloc[0][f"{model}_intercept"], decimals)
                cells += [s, b, "", ""]

    elif param_type == "gamma":
        row = _get(df, target, pdg, cat, qty)
        for model in MODEL_NAMES:
            if row.empty:
                cells += ["---", "---", "---", "---"]
            else:
                s  = _fmt(row.iloc[0][f"{model}_slope"],     decimals)
                b  = _fmt(row.iloc[0][f"{model}_intercept"], decimals)
                fl = _fmt(row.iloc[0].get(f"{model}_floor", 0.0), decimals)
                ex = _fmt(row.iloc[0].get(f"{model}_exp",   1.0), decimals)
                cells += [s, b, ex, fl]

    return " & ".join(cells) + r" \\"


THRESH_LABELS = {
    0: "without a threshold",
    1: "with a threshold",
    2: "with a threshold and compound nuclei",
}

PDG_ORDER = [2212, 2112, 211, 111, -211]  # p, n, pi+, pi0, pi-


def csv_to_latex_unified(df, decimals=4, thresh=None):
    """
    Produce one table* per nucleus.  Within each table, particles are row groups
    separated by \\hline.  Each particle has one sub-row per slot.
    Nucleon Sum rows use all 4 gamma sub-columns; linear rows leave floor/power blank.

    thresh: None | 0 | 1 | 2
        Appends a phrase to each caption describing the threshold condition.
    """
    pdgs = _particles_in(df)
    ordered_pdgs = [p for p in PDG_ORDER if p in pdgs]

    n_models = len(MODEL_NAMES)
    # total data columns = 4 sub-columns × n_models
    n_data_cols = 4 * n_models

    # tabular column spec:
    #   col 1: particle label  (centered, fixed width)
    #   col 2: fit-type label
    #   cols 3..: data (4 per model, no separator between sub-cols within a model,
    #             vertical rule between models)
    inner_model_spec = "cccc"
    model_specs = ("|" + inner_model_spec) * n_models + "|"
    col_spec = r"|>{\centering\arraybackslash}p{1.2cm}|l" + model_specs

    # Sub-column header labels — printed only when BOTH linear and gamma slots present
    subheader_labels = {
        "linear": [r"$m$", r"$b$", r"", r""],
        "gamma":  [r"$m$", r"$b$", r"$c$", r"$\nu$"],
    }

    tables = []

    for target, (nuc_ascii, nuc_tex, nuc_id) in NUCLEUS_NAMES.items():
        L = []

        L.append(r"\begin{table*}[htbp]")
        L.append(r"\centering")
        thresh_phrase = (
            f" Results are {THRESH_LABELS[thresh]}." if thresh in THRESH_LABELS else ""
        )
        caption = (
            f"\\caption{{Fit parameters for hadronic rescattering models on {nuc_tex}."
            f"{thresh_phrase} "
            "Linear fits use $y = j + iT$. "
            "Nucleon sum distributions use $\\Gamma = j \\times \\exp(i \\times T^{k}) + l$.}"
        )
        L.append(caption)
        L.append(rf"\label{{tab:fitparams_{nuc_id}}}")
        L.append(rf"\begin{{tabular}}{{{col_spec}}}")
        L.append(r"\hline")

        # --- Model header ---
        header_cells = [r"\multicolumn{2}{c|}{}"]
        for model in MODEL_NAMES:
            ms = MODEL_SHORT[model]
            header_cells.append(rf"\multicolumn{{4}}{{c|}}{{\textbf{{{ms}}}}}")
        L.append(" & ".join(header_cells) + r" \\")
        L.append(r"\hline")

        # --- Sub-column header ---
        sub_cells = ["", ""]
        for _ in MODEL_NAMES:
            sub_cells += [r"$i$", r"$j$", r"$k$", r"$l$"]
        L.append(" & ".join(sub_cells) + r" \\")
        L.append(r"\hline\hline")

        # --- Particle row groups ---
        for pdg_idx, pdg in enumerate(ordered_pdgs):
            info       = PARTICLE_INFO[pdg]
            probe_type = info["probe_type"]
            label      = info["label"]
            slots      = PION_SLOTS if probe_type == "pion" else NUCLEON_SLOTS
            n_slots    = len(slots)

            for slot_idx, (cat, qty, field, row_label, param_type) in enumerate(slots):
                # First sub-row gets the particle multirow label
                if slot_idx == 0:
                    particle_cell = rf"\multirow{{{n_slots}}}{{*}}{{{label}}}"
                else:
                    particle_cell = ""

                data_row = _data_row(df, target, pdg, cat, qty,
                                     param_type, row_label, decimals)

                L.append(f"{particle_cell} & {data_row}")

                # Thin cline between sub-rows within a particle group
                if slot_idx < n_slots - 1:
                    L.append(rf"\cline{{2-{2 + n_data_cols}}}")

            # Thick hline between particle groups (and after last)
            L.append(r"\hline")

        L.append(r"\end{tabular}")
        L.append(r"\end{table*}")

        tables.append("\n".join(L))

    return "\n\n".join(tables)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csvfiles", nargs="+")
    parser.add_argument("--latex",    action="store_true")
    parser.add_argument("--cpp",      action="store_true")
    parser.add_argument("--decimals", type=int, default=4)
    parser.add_argument("--thresh", type=int, default=None,
                        choices=[0, 1, 2],
                        help="0: without a threshold, "
                             "1: with a threshold, "
                             "2: with a threshold and compound nuclei")
    args = parser.parse_args()

    df = pd.concat([_load_csv(p) for p in args.csvfiles], ignore_index=True)

    if args.cpp or not args.latex:
        print(SHARED_STRUCTS)
        print(csv_to_cpp(df))

    if args.latex:
        print(csv_to_latex_unified(df, decimals=args.decimals, thresh=args.thresh))