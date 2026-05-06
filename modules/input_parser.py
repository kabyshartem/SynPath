"""
Input Parser — normalizes SMILES, detects molecule type (organic / peptide),
fetches basic info from PubChem.
"""
import re
import httpx
from typing import Optional

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

AA_MAP = {
    "G": {"name": "Glycine",       "fmoc": "Fmoc-Gly-OH",          "smiles": "NCC(=O)O"},
    "A": {"name": "Alanine",       "fmoc": "Fmoc-Ala-OH",          "smiles": "NC(C)C(=O)O"},
    "V": {"name": "Valine",        "fmoc": "Fmoc-Val-OH",          "smiles": "NC(C(C)C)C(=O)O"},
    "L": {"name": "Leucine",       "fmoc": "Fmoc-Leu-OH",          "smiles": "NC(CC(C)C)C(=O)O"},
    "I": {"name": "Isoleucine",    "fmoc": "Fmoc-Ile-OH",          "smiles": "NC(C(C)CC)C(=O)O"},
    "P": {"name": "Proline",       "fmoc": "Fmoc-Pro-OH",          "smiles": "N1CCCC1C(=O)O"},
    "F": {"name": "Phenylalanine", "fmoc": "Fmoc-Phe-OH",          "smiles": "NC(Cc1ccccc1)C(=O)O"},
    "W": {"name": "Tryptophan",    "fmoc": "Fmoc-Trp(Boc)-OH",     "smiles": "NC(Cc1c[nH]c2ccccc12)C(=O)O"},
    "M": {"name": "Methionine",    "fmoc": "Fmoc-Met-OH",          "smiles": "NC(CCSC)C(=O)O"},
    "S": {"name": "Serine",        "fmoc": "Fmoc-Ser(tBu)-OH",     "smiles": "NC(CO)C(=O)O"},
    "T": {"name": "Threonine",     "fmoc": "Fmoc-Thr(tBu)-OH",     "smiles": "NC(C(O)C)C(=O)O"},
    "C": {"name": "Cysteine",      "fmoc": "Fmoc-Cys(Trt)-OH",     "smiles": "NC(CS)C(=O)O"},
    "Y": {"name": "Tyrosine",      "fmoc": "Fmoc-Tyr(tBu)-OH",     "smiles": "NC(Cc1ccc(O)cc1)C(=O)O"},
    "N": {"name": "Asparagine",    "fmoc": "Fmoc-Asn(Trt)-OH",     "smiles": "NC(CC(=O)N)C(=O)O"},
    "Q": {"name": "Glutamine",     "fmoc": "Fmoc-Gln(Trt)-OH",     "smiles": "NC(CCC(=O)N)C(=O)O"},
    "D": {"name": "Aspartate",     "fmoc": "Fmoc-Asp(OtBu)-OH",    "smiles": "NC(CC(=O)O)C(=O)O"},
    "E": {"name": "Glutamate",     "fmoc": "Fmoc-Glu(OtBu)-OH",    "smiles": "NC(CCC(=O)O)C(=O)O"},
    "K": {"name": "Lysine",        "fmoc": "Fmoc-Lys(Boc)-OH",     "smiles": "NC(CCCCN)C(=O)O"},
    "R": {"name": "Arginine",      "fmoc": "Fmoc-Arg(Pbf)-OH",     "smiles": "NC(CCCNC(=N)N)C(=O)O"},
    "H": {"name": "Histidine",     "fmoc": "Fmoc-His(Trt)-OH",     "smiles": "NC(Cc1cnc[nH]1)C(=O)O"},
}

PROBLEMATIC_SEQUENCES = {
    "DP": "Asp-Pro — high risk of aspartimide formation. Use Fmoc-Asp(OMpe)-OH.",
    "DG": "Asp-Gly — aspartimide risk. Add 0.1M HOBt in DMF during coupling.",
    "DS": "Asp-Ser — aspartimide risk.",
    "DA": "Asp-Ala — moderate aspartimide risk.",
    "NG": "Asn-Gly — deamidation risk.",
    "NS": "Asn-Ser — deamidation risk.",
    "WW": "Trp-Trp — aggregation prone.",
    "VV": "Val-Val — aggregation prone, consider pseudoproline dipeptide.",
    "II": "Ile-Ile — aggregation prone.",
    "PP": "Pro-Pro — difficult coupling, may need double coupling.",
}

# Common compound names that look like AA sequences — always treat as name
KNOWN_COMPOUNDS = {
    "aspirin","ibuprofen","caffeine","paracetamol","acetaminophen",
    "lidocaine","resveratrol","taxol","paclitaxel","morphine","aniline",
    "codeine","dopamine","serotonin","adrenaline","epinephrine",
    "glucose","fructose","sucrose","cholesterol","testosterone",
    "estrogen","progesterone","penicillin","amoxicillin",
    "benzene","toluene","phenol","naphthalene","anthracene",
    "sildenafil","tamiflu","oseltamivir","metformin","salicylic acid",
    "atorvastatin","omeprazole","acetic acid","ethanol","methanol",
}


def detect_input_type(raw: str) -> dict:
    raw = raw.strip()

    # 1. Known compound names — never treat as peptide
    if raw.lower() in KNOWN_COMPOUNDS:
        return {"type": "name", "value": raw}

    # 2. Contains spaces or mixed case without dashes → name
    if " " in raw:
        return {"type": "name", "value": raw}

    # 3. Three-letter peptide: Ala-Gly-Phe or H-Ala-Gly-OH
    three_letter = re.findall(r'[A-Z][a-z]{2}', raw)
    aa_three = {"Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
                "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"}
    if len(three_letter) >= 2 and all(aa in aa_three for aa in three_letter):
        return {"type": "peptide_sequence", "value": raw, "format": "three_letter"}

    # 4. SMILES heuristic
    smiles_chars = set("CNOSPFClBrI()[]=#@+\\/-0123456789")
    if len(raw) > 3 and sum(c in smiles_chars for c in raw) / len(raw) > 0.6:
        return {"type": "smiles", "value": raw}

    # 5. Pure uppercase ALL amino acid letters, 3+ chars → peptide sequence
    if re.match(r'^[ACDEFGHIKLMNPQRSTVWY]{3,}$', raw) and raw.isupper():
        return {"type": "peptide_sequence", "value": raw, "format": "single_letter"}

    # 6. Default → treat as compound name
    return {"type": "name", "value": raw}


def parse_peptide_sequence(seq: str, fmt: str = "single_letter") -> dict:
    if fmt == "single_letter":
        residues = list(seq.upper())
    else:
        three_to_one = {
            "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q",
            "Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K",
            "Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
            "Tyr":"Y","Val":"V"
        }
        found = re.findall(r'[A-Z][a-z]{2}', seq)
        residues = [three_to_one.get(aa, "?") for aa in found]

    result = []
    for aa in residues:
        info = AA_MAP.get(aa, {"name": f"Unknown({aa})", "fmoc": f"Fmoc-{aa}-OH", "smiles": ""})
        result.append({"code": aa, **info})

    warnings = []
    seq_str = "".join(residues)
    for motif, warning in PROBLEMATIC_SEQUENCES.items():
        if motif in seq_str:
            warnings.append(warning)

    return {
        "residues": result,
        "length": len(result),
        "sequence": seq_str,
        "warnings": warnings,
        "mw_estimate": len(result) * 111.1
    }


async def pubchem_lookup(identifier: str, id_type: str = "name") -> Optional[dict]:
    try:
        async with httpx.AsyncClient(timeout=10) as client:
            url = f"{PUBCHEM_URL}/compound/{id_type}/{identifier}/cids/JSON"
            r = await client.get(url)
            if r.status_code != 200:
                return None
            cids = r.json().get("IdentifierList", {}).get("CID", [])
            if not cids:
                return None
            cid = cids[0]

            props_url = f"{PUBCHEM_URL}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,IsomericSMILES,InChIKey/JSON"
            r2 = await client.get(props_url)
            if r2.status_code != 200:
                return None
            props = r2.json().get("PropertyTable", {}).get("Properties", [{}])[0]

            return {
                "cid": cid,
                "smiles": props.get("IsomericSMILES", ""),
                "formula": props.get("MolecularFormula", ""),
                "mw": props.get("MolecularWeight", ""),
                "iupac": props.get("IUPACName", ""),
                "inchikey": props.get("InChIKey", ""),
                "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            }
    except Exception:
        return None


async def smiles_to_info(smiles: str) -> Optional[dict]:
    return await pubchem_lookup(smiles, "smiles")


async def name_to_smiles(name: str) -> Optional[dict]:
    return await pubchem_lookup(name, "name")
