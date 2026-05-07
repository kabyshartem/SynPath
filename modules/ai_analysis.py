import anthropic, json, sys, os, re
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config

import httpx as _httpx

import hashlib as _hashlib
_ROUTE_CACHE: dict = {}
def _cache_key(smiles: str, partial: bool, start: str) -> str:
    return _hashlib.md5(f"{smiles}|{partial}|{start}".encode()).hexdigest()


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


async def _verify_doi(doi: str) -> bool:
    """Check if DOI exists via CrossRef API (free, no key needed)."""
    if not doi or len(doi) < 5:
        return False
    try:
        async with _httpx.AsyncClient(timeout=5) as c:
            r = await c.get(f"https://api.crossref.org/works/{doi}", 
                          headers={"User-Agent": "SynPath/1.0 (mailto:synpath@example.com)"})
            return r.status_code == 200
    except Exception:
        return False

async def _clean_dois(parsed: dict) -> dict:
    """Remove fabricated DOIs — only keep verified ones."""
    for route in parsed.get("routes", []):
        for step in route.get("steps_with_smiles", []):
            doi = step.get("doi", "")
            if doi and not await _verify_doi(doi):
                step["doi"] = ""
        for ref in route.get("literature_support", []):
            doi = ref.get("doi", "")
            if doi and not await _verify_doi(doi):
                ref["doi"] = ""
    return parsed

async def _fix_route_smiles(parsed: dict) -> dict:
    for route in parsed.get("routes", []):
        for step in route.get("steps_with_smiles", []):
            step["reactant_smiles"] = await _fix_smiles(step.get("reactant_name",""), step.get("reactant_smiles",""))
            step["product_smiles"] = await _fix_smiles(step.get("product_name",""), step.get("product_smiles",""))
    return parsed



SYSTEM = "You are SynPath A, expert synthetic chemist. Generate 2-4 retrosynthetic routes (fewer if molecule is simple — quality over quantity). Routes must lead to EXACT target. The final product SMILES in every route MUST match the target SMILES character-for-character. No analogs, no regioisomers, no stereoisomers as final product. CRITICAL DOI RULE: Never fabricate DOIs. For each step provide a real DOI where: (1) the paper uses this EXACT compound, OR (2) the paper uses the SAME reaction on a structurally similar substrate (same ring system, same functional group). DOI must be a real published paper you are certain exists. If you are not 100% certain — leave doi empty. Never guess DOIs. Return ONLY valid JSON, no markdown. Be concise — short pros/cons (max 8 words each), short conditions, no long explanations."

async def analyze_routes(target_name, target_smiles, retrosynthesis_data, literature_papers, availability_data, partial_route=False, start_smiles=""):
    _ck = _cache_key(target_smiles, partial_route, start_smiles)
    if _ck in _ROUTE_CACHE:
        return _ROUTE_CACHE[_ck]
    client = anthropic.Anthropic(api_key=config.ANTHROPIC_API_KEY)
    msg = f"TARGET: {target_name}\nSMILES: {target_smiles}\nWARNING: Final product in every route must have SMILES={target_smiles}\n"
    if partial_route and start_smiles:
        msg += f"USER HAS INTERMEDIATE: {start_smiles}\n"
    valid_lit = [p for p in (literature_papers or [])[:3] if not p.get("error")]
    if valid_lit:
        msg += "LITERATURE:\n"
        for i,p in enumerate(valid_lit):
            c = p.get("extracted_conditions",{})
            msg += f"[{i+1}] {p.get('title','')[:50]}\nDOI: {p.get('doi','')}\nConditions: temp={c.get('temperature','?')}, yield={c.get('yield','?')}\nAbstract: {p.get('abstract_snippet','')[:80]}\n"
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
        if not parsed: print(f"PARSE FAILED: {raw[:200] if raw else None}", flush=True); return {"error":f"Parse failed: {raw[:300]}","routes":[],"analysis_summary":"Parse error"}
        if not isinstance(parsed.get("routes"),list) or not parsed["routes"]: return {"error":"No routes","routes":[],"analysis_summary":parsed.get("analysis_summary","")}
        parsed = await _fix_route_smiles(parsed)
        parsed = await _clean_dois(parsed)
        _ROUTE_CACHE[_ck] = parsed
        return parsed
    except anthropic.APIStatusError as e: print(f"API ERROR: {e.status_code} {e.message}", flush=True); return {"error":f"API {e.status_code}: {e.message}","routes":[],"analysis_summary":str(e)}
    except Exception as ex: print(f"EXCEPTION: {str(ex)}", flush=True); return {"error":str(ex),"routes":[],"analysis_summary":str(ex)}

async def analyze_peptide(sequence, spps_protocol, literature_papers, availability_data):
    client = anthropic.Anthropic(api_key=config.ANTHROPIC_API_KEY)
    try:
        resp = client.messages.create(model="claude-haiku-4-5-20251001", max_tokens=1500, messages=[{"role":"user","content":f"PEPTIDE: {sequence}\nLENGTH: {spps_protocol.get('length')} residues\nSCALE: {spps_protocol.get('scale')}\nProvide expert Fmoc SPPS analysis."}])
        return {"analysis":resp.content[0].text}
    except Exception as e: return {"error":str(e)}
