import anthropic, json, sys, os, re
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config

import httpx as _httpx

async def _fix_smiles(name: str, smiles: str) -> str:
    """Validate SMILES via RDKit; if invalid/empty, fetch from PubChem by name."""
    if smiles:
        try:
            from rdkit import Chem
            if Chem.MolFromSmiles(smiles) is not None:
                return smiles
        except Exception:
            pass
    if not name or len(name) < 2:
        return smiles
    try:
        async with _httpx.AsyncClient(timeout=10) as c:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{_httpx.URL(name)}/property/IsomericSMILES/JSON"
            r = await c.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON")
            if r.status_code == 200:
                data = r.json()
                return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except Exception:
        pass
    return smiles

async def _fix_route_smiles(parsed: dict) -> dict:
    for route in parsed.get("routes", []):
        for step in route.get("steps_with_smiles", []):
            step["reactant_smiles"] = await _fix_smiles(step.get("reactant_name",""), step.get("reactant_smiles",""))
            step["product_smiles"] = await _fix_smiles(step.get("product_name",""), step.get("product_smiles",""))
    return parsed



SYSTEM = "You are SynPath A, expert synthetic chemist. Generate 4-5 retrosynthetic routes for the EXACT target molecule provided. Every route must lead to the EXACT target — not an analog, not a similar compound. If analog-based reasoning is used, the final product must still be the exact target. Use chemistry knowledge when no literature exists. Return ONLY valid JSON object, no markdown, no text outside JSON."

async def analyze_routes(target_name, target_smiles, retrosynthesis_data, literature_papers, availability_data, partial_route=False, start_smiles=""):
    client = anthropic.Anthropic(api_key=config.ANTHROPIC_API_KEY)
    msg = f"TARGET: {target_name}\nSMILES: {target_smiles}\n"
    if partial_route and start_smiles:
        msg += f"USER HAS INTERMEDIATE: {start_smiles}\n"
    valid_lit = [p for p in (literature_papers or [])[:4] if not p.get("error")]
    if valid_lit:
        msg += "LITERATURE:\n"
        for i,p in enumerate(valid_lit):
            c = p.get("extracted_conditions",{})
            msg += f"[{i+1}] {p.get('title','')[:80]}\nDOI: {p.get('doi','')}\nConditions: temp={c.get('temperature','?')}, yield={c.get('yield','?')}\nAbstract: {p.get('abstract_snippet','')[:150]}\n"
    msg += """\nReturn JSON: {"target":"name","smiles":"SMILES","analysis_summary":"text","routes":[{"rank":1,"name":"name","num_steps":2,"overall_yield_estimate":"55%","difficulty":"moderate","cost_estimate":"low","scalability":"good","key_reactions":["rxn1"],"starting_materials":["sm1"],"steps_with_smiles":[{"step":1,"label":"A","reaction":"name","reagents":["r1 (1.0 eq)"],"conditions":"cond","reactant_smiles":"SMILES","reactant_name":"name","product_smiles":"SMILES","product_name":"name","yield":"80%","doi":""}],"pros":["p1"],"cons":["c1"],"literature_support":[{"doi":"","notes":""}]}],"recommended_route":1,"recommendation_rationale":"reason","detailed_procedure":{"title":"title","steps":[{"step":1,"title":"t","reaction":"r","reagents":["r1"],"conditions":"c","workup":"w","purification":"p","yield":"y","doi_reference":"","notes":"n"}]}}"""
    try:
        resp = client.messages.create(model="claude-sonnet-4-6", max_tokens=8192, temperature=0, system=SYSTEM, messages=[{"role":"user","content":msg}])
        raw = resp.content[0].text.strip()
        raw = re.sub(r'^```[a-z]*\s*','',raw); raw = re.sub(r'\s*```$','',raw).strip()
        parsed = None
        try: parsed = json.loads(raw)
        except: pass
        if not parsed:
            try: parsed = json.loads(raw[raw.index("{"):raw.rindex("}")+1])
            except: pass
        if not parsed:
            m = re.search(r'\{.*\}',raw,re.DOTALL)
            if m:
                try: parsed = json.loads(m.group())
                except: pass
        if not parsed: return {"error":f"Parse failed: {raw[:300]}","routes":[],"analysis_summary":"Parse error"}
        if not isinstance(parsed.get("routes"),list) or not parsed["routes"]: return {"error":"No routes","routes":[],"analysis_summary":parsed.get("analysis_summary","")}
        parsed = await _fix_route_smiles(parsed)
        return parsed
    except anthropic.APIStatusError as e: return {"error":f"API {e.status_code}: {e.message}","routes":[],"analysis_summary":str(e)}
    except Exception as ex: return {"error":str(ex),"routes":[],"analysis_summary":str(ex)}

async def analyze_peptide(sequence, spps_protocol, literature_papers, availability_data):
    client = anthropic.Anthropic(api_key=config.ANTHROPIC_API_KEY)
    try:
        resp = client.messages.create(model="claude-sonnet-4-6", max_tokens=2048, messages=[{"role":"user","content":f"PEPTIDE: {sequence}\nLENGTH: {spps_protocol.get('length')} residues\nSCALE: {spps_protocol.get('scale')}\nProvide expert Fmoc SPPS analysis."}])
        return {"analysis":resp.content[0].text}
    except Exception as e: return {"error":str(e)}
