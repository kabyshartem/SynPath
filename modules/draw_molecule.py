"""Draw molecules and reaction schemes using RDKit"""
import base64, io
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False


def smiles_to_png_base64(smiles: str, width=300, height=200, bg_white=True) -> str:
    if not RDKIT_OK or not smiles: return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return ""
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().padding = 0.15
        if bg_white:
            drawer.drawOptions().backgroundColour = (1,1,1,1)
        else:
            drawer.drawOptions().backgroundColour = (1,1,1,0)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("utf-8")
    except Exception:
        return smiles_to_png_base64_pubchem(smiles, width, height)


def smiles_to_png_base64_transparent(smiles: str, width=200, height=150) -> str:
    return smiles_to_png_base64(smiles, width, height, bg_white=False)


def reaction_smiles_to_png_base64(rxn_smiles: str, width=600, height=200) -> str:
    if not RDKIT_OK or not rxn_smiles: return ""
    try:
        from rdkit.Chem import rdChemReactions
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None: return ""
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("utf-8")
    except Exception:
        return ""

def smiles_to_png_base64_pubchem(smiles: str, width: int = 300, height: int = 200) -> str:
    """Fallback: fetch molecule image from PubChem when RDKit fails."""
    try:
        import urllib.request, base64, urllib.parse
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles)}/PNG?image_size={width}x{height}"
        with urllib.request.urlopen(url, timeout=8) as r:
            return base64.b64encode(r.read()).decode("utf-8")
    except Exception:
        return ""
