"""
Peptide SPPS Module
Fmoc/Boc solid-phase peptide synthesis route generator
"""
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from modules.input_parser import AA_MAP, PROBLEMATIC_SEQUENCES


def select_resin(sequence: str, c_terminus: str = "amide") -> dict:
    """Select appropriate resin based on desired C-terminus."""
    if c_terminus == "amide":
        return {
            "name": "Rink Amide MBHA Resin",
            "loading": "0.3-0.7 mmol/g typical",
            "cleavage": "TFA cocktail",
            "product": "C-terminal amide",
            "notes": "Most common for peptide amides. Pre-loaded available."
        }
    elif c_terminus == "acid":
        return {
            "name": "Wang Resin",
            "loading": "0.4-1.1 mmol/g typical",
            "cleavage": "TFA cocktail",
            "product": "C-terminal free acid",
            "notes": "Standard resin for C-terminal acids."
        }
    elif c_terminus == "fragment":
        return {
            "name": "2-Chlorotrityl (2-CTC) Resin",
            "loading": "0.5-1.6 mmol/g typical",
            "cleavage": "1% TFA/DCM or AcOH/TFE/DCM 1:1:8",
            "product": "Side-chain protected fragment (for solution coupling)",
            "notes": "Ideal for fragment condensation. Mild cleavage preserves side-chain PGs."
        }
    else:
        return {
            "name": "Rink Amide MBHA Resin",
            "loading": "0.3-0.7 mmol/g typical",
            "cleavage": "TFA cocktail",
            "product": "C-terminal amide",
            "notes": "Default selection."
        }


def select_coupling_protocol(aa_code: str, position: int, previous_aa: str = "") -> dict:
    """Select coupling reagent based on amino acid and context."""
    # Difficult couplings need stronger reagents
    difficult = {"His", "Arg", "Cys", "Pro", "Aib"}
    beta_branched = {"V", "I", "T"}

    aa_info = AA_MAP.get(aa_code, {})
    aa_name = aa_info.get("name", "")

    if aa_code in beta_branched or (previous_aa in beta_branched and aa_code in beta_branched):
        return {
            "reagent": "HATU / DIPEA",
            "ratio": "3.0 eq AA / 2.9 eq HATU / 6.0 eq DIPEA",
            "solvent": "DMF",
            "time": "45 min, RT",
            "notes": "Beta-branched AA — use HATU for reliable coupling. Consider double coupling."
        }
    elif aa_code == "P":
        return {
            "reagent": "HATU / DIPEA",
            "ratio": "3.0 eq AA / 2.9 eq HATU / 6.0 eq DIPEA",
            "solvent": "DMF",
            "time": "60 min, RT — double coupling recommended",
            "notes": "Proline is a secondary amine — slow coupling. Always double couple."
        }
    elif aa_code == "R":
        return {
            "reagent": "HATU / DIPEA",
            "ratio": "3.0 eq Fmoc-Arg(Pbf)-OH / 2.9 eq HATU / 6.0 eq DIPEA",
            "solvent": "DMF",
            "time": "45 min, RT",
            "notes": "Arg(Pbf) can be slow to dissolve — pre-activate for 2 min before adding resin."
        }
    else:
        return {
            "reagent": "DIC / Oxyma Pure",
            "ratio": "3.0 eq AA / 3.0 eq DIC / 3.0 eq Oxyma",
            "solvent": "DMF",
            "time": "30 min, RT",
            "notes": "Standard green coupling protocol. Low epimerization risk."
        }


def generate_cleavage_cocktail(sequence: str) -> dict:
    """Generate TFA cleavage cocktail based on amino acid composition."""
    has_trp = "W" in sequence
    has_met = "M" in sequence
    has_cys = "C" in sequence
    has_arg = "R" in sequence

    if has_cys and has_met:
        return {
            "cocktail": "TFA / TIS / H2O / EDT",
            "ratio": "92.5 : 2.5 : 2.5 : 2.5",
            "time": "2-3 h, RT",
            "notes": "EDT added as additional scavenger for Cys and Met."
        }
    elif has_trp:
        return {
            "cocktail": "TFA / TIS / H2O / phenol",
            "ratio": "87.5 : 5 : 5 : 2.5",
            "time": "2-3 h, RT",
            "notes": "Phenol added to protect Trp from alkylation by tBu cations."
        }
    elif has_arg:
        return {
            "cocktail": "TFA / TIS / H2O",
            "ratio": "95 : 2.5 : 2.5",
            "time": "2-3 h, RT",
            "notes": "Pbf removal from Arg can be slow — extend cleavage time to 3h if needed."
        }
    else:
        return {
            "cocktail": "TFA / TIS / H2O",
            "ratio": "95 : 2.5 : 2.5",
            "time": "2-3 h, RT",
            "notes": "Standard cocktail for most peptides."
        }


def generate_spps_protocol(
    residues: list,
    scale_mmol: float = 0.1,
    chemistry: str = "fmoc",
    c_terminus: str = "amide"
) -> dict:
    """
    Generate complete SPPS protocol for a peptide sequence.
    residues: list of dicts from parse_peptide_sequence()
    scale_mmol: synthesis scale in mmol
    chemistry: 'fmoc' or 'boc'
    c_terminus: 'amide', 'acid', or 'fragment'
    """
    sequence = "".join(r["code"] for r in residues)
    resin = select_resin(sequence, c_terminus)

    # Calculate resin amount (assuming 0.5 mmol/g loading)
    loading = 0.5  # mmol/g
    resin_amount_g = scale_mmol / loading

    steps = []

    # Step 0: Resin swelling
    steps.append({
        "step": 0,
        "title": "Resin Swelling",
        "details": f"Weigh {resin_amount_g:.2f} g {resin['name']} into synthesis vessel. "
                   f"Add DMF (10 mL/g resin), let swell 30 min with gentle agitation.",
        "time": "30 min",
        "reagents": [f"{resin['name']}: {resin_amount_g:.2f} g", "DMF: ~5 mL"]
    })

    if chemistry == "fmoc":
        # Step 1: Initial Fmoc deprotection
        steps.append({
            "step": 1,
            "title": "Initial Fmoc Deprotection",
            "details": "Add 20% piperidine in DMF (5 mL). Agitate 5 min, drain. "
                      "Repeat with fresh 20% piperidine/DMF for 15 min. "
                      "Wash with DMF × 5.",
            "time": "20 min total",
            "reagents": ["20% piperidine in DMF: 10 mL total", "DMF: 25 mL (washes)"],
            "monitoring": "Kaiser test should be positive (blue = free amine)"
        })

    # Coupling steps (C→N direction for SPPS)
    for i, aa in enumerate(reversed(residues)):
        aa_code = aa["code"]
        prev_code = residues[len(residues) - i - 2]["code"] if i < len(residues) - 1 else ""
        coupling = select_coupling_protocol(aa_code, i + 1, prev_code)

        fmoc_aa = aa.get("fmoc", f"Fmoc-{aa_code}-OH")
        mw_approx = 400  # approximate Fmoc-AA-OH MW
        mass_g = scale_mmol * 3.0 * mw_approx / 1000  # 3 eq

        step_num = (i * 2) + 2

        # Coupling
        steps.append({
            "step": step_num,
            "title": f"Coupling: {fmoc_aa} (Position {len(residues) - i})",
            "details": f"Dissolve {fmoc_aa} ({mass_g:.2f} g, ~3.0 eq) in DMF. "
                      f"Add {coupling['reagent']} ({coupling['ratio']}). "
                      f"Pre-activate 2 min, add to resin. {coupling['time']}. "
                      f"Drain, wash DMF × 3.",
            "time": coupling["time"],
            "reagents": [f"{fmoc_aa}: ~{mass_g:.2f} g (3.0 eq)", coupling["reagent"], "DMF"],
            "monitoring": "Kaiser test should be negative (yellow/clear = complete coupling)",
            "notes": coupling.get("notes", "")
        })

        # Fmoc deprotection (skip after last AA if Boc used for N-terminus)
        if i < len(residues) - 1:
            steps.append({
                "step": step_num + 1,
                "title": f"Fmoc Deprotection (after position {len(residues) - i})",
                "details": "20% piperidine/DMF, 5 min + 15 min. Wash DMF × 5.",
                "time": "20 min",
                "reagents": ["20% piperidine in DMF: 10 mL", "DMF: 25 mL"],
                "monitoring": "Kaiser test positive"
            })

    # Final deprotection and cleavage
    cocktail = generate_cleavage_cocktail(sequence)
    steps.append({
        "step": len(steps) + 1,
        "title": "Cleavage & Global Deprotection",
        "details": f"Add cleavage cocktail ({cocktail['cocktail']}, {cocktail['ratio']}) "
                  f"to dry resin (~1 mL/100 mg resin). "
                  f"Agitate at RT for {cocktail['time']}. "
                  f"Filter resin, wash with neat TFA × 2. "
                  f"Precipitate peptide in cold diethyl ether (45 mL per 5 mL TFA). "
                  f"Centrifuge 4000 rpm × 10 min. Decant ether, dry pellet.",
        "time": cocktail["time"],
        "reagents": [
            f"TFA: main component",
            f"TIS: scavenger",
            f"H2O: scavenger",
            "Cold Et2O: ~50 mL for precipitation"
        ],
        "notes": cocktail.get("notes", "")
    })

    # Purification recommendation
    steps.append({
        "step": len(steps) + 1,
        "title": "Purification — Reverse Phase HPLC",
        "details": "Dissolve crude peptide in 0.1% TFA / H2O (+ MeCN if needed). "
                  "Inject on C18 column. Gradient: 5→95% MeCN in 0.1% TFA / H2O, 30 min. "
                  "Collect main peak, lyophilize.",
        "time": "1-2 days (HPLC + lyophilization)",
        "reagents": ["0.1% TFA in H2O", "0.1% TFA in MeCN", "C18 column"],
        "notes": "Confirm identity by LCMS. Expected purity >95% after RP-HPLC."
    })

    # Check for problematic sequences
    warnings = []
    for motif, warning in PROBLEMATIC_SEQUENCES.items():
        if motif in sequence:
            warnings.append({"motif": motif, "warning": warning})

    return {
        "method": "Fmoc SPPS" if chemistry == "fmoc" else "Boc SPPS",
        "scale": f"{scale_mmol} mmol",
        "sequence": sequence,
        "length": len(residues),
        "resin": resin,
        "cleavage_cocktail": cocktail,
        "steps": steps,
        "warnings": warnings,
        "estimated_crude_yield": f"{scale_mmol * len(residues) * 111.1:.0f} mg theoretical",
        "notes": f"Total {len(steps)} steps. Estimated synthesis time: {len(residues) * 1.5:.0f}h hands-on."
    }
