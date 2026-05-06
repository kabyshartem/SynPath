"""
Literature Mining Layer
Layer 3: PubMed / Europe PMC / PubChem for real synthesis procedures
         Returns actual conditions, yields, DOI references
"""
import httpx
import re
from typing import Optional

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
PUBCHEM   = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBMED    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


async def search_synthesis_literature(compound_name: str, smiles: str = "") -> list:
    """
    Search Europe PMC for synthesis papers on target compound.
    Returns list of papers with DOI, title, abstract, conditions if extractable.
    """
    papers = []

    queries = [
        f'"{compound_name}" synthesis',
        f'"{compound_name}" preparation',
    ]
    if smiles:
        queries.append(f'"{compound_name}" total synthesis')

    try:
        async with httpx.AsyncClient(timeout=15) as client:
            for query in queries[:2]:
                url = f"{EUROPEPMC}/search"
                params = {
                    "query": query,
                    "format": "json",
                    "pageSize": 10,
                    "sort": "CITED desc",
                    "resultType": "core"
                }
                r = await client.get(url, params=params)
                if r.status_code != 200:
                    continue

                results = r.json().get("resultList", {}).get("result", [])
                for paper in results:
                    doi = paper.get("doi", "")
                    title = paper.get("title", "")
                    abstract = paper.get("abstractText", "")

                    # Only include papers that likely contain synthesis info
                    synth_keywords = ["synthesized", "synthesis", "prepared", "yield",
                                     "reaction", "procedure", "method"]
                    abstract_lower = abstract.lower()
                    if not any(kw in abstract_lower for kw in synth_keywords):
                        continue

                    # Extract conditions from abstract using patterns
                    conditions = extract_conditions_from_text(abstract)

                    papers.append({
                        "title": title,
                        "doi": doi,
                        "doi_url": f"https://doi.org/{doi}" if doi else "",
                        "journal": paper.get("journalTitle", ""),
                        "year": paper.get("pubYear", ""),
                        "authors": paper.get("authorString", ""),
                        "abstract_snippet": abstract[:400] + "..." if len(abstract) > 400 else abstract,
                        "extracted_conditions": conditions,
                        "pmid": paper.get("pmid", ""),
                        "source": "Europe PMC"
                    })

            # Deduplicate by DOI
            seen = set()
            unique = []
            for p in papers:
                key = p["doi"] or p["title"]
                if key not in seen:
                    seen.add(key)
                    unique.append(p)

            return unique[:8]

    except Exception as e:
        return [{"error": str(e)}]


def extract_conditions_from_text(text: str) -> dict:
    """
    Extract synthesis conditions from text using regex patterns.
    Finds temperature, solvents, time, yield.
    """
    conditions = {}

    # Temperature patterns: "80°C", "reflux", "0 °C", "-78 °C"
    temp_patterns = [
        r'(-?\d+)\s*°C',
        r'(-?\d+)\s*degrees?\s*C',
        r'(reflux)',
        r'(room temperature|rt|RT|r\.t\.)',
    ]
    temps = []
    for pat in temp_patterns:
        found = re.findall(pat, text, re.IGNORECASE)
        temps.extend(found)
    if temps:
        conditions["temperature"] = temps[0] if len(temps) == 1 else temps

    # Yield patterns: "85% yield", "yield of 72%", "obtained in 65%"
    yield_pat = re.findall(r'(\d{1,3})\s*%\s*(?:yield|isolated|obtained)', text, re.IGNORECASE)
    if not yield_pat:
        yield_pat = re.findall(r'yield(?:ed|s)?\s+(?:of\s+)?(\d{1,3})\s*%', text, re.IGNORECASE)
    if yield_pat:
        conditions["yield"] = yield_pat[0] + "%"

    # Time patterns: "12 h", "overnight", "2 hours"
    time_pat = re.findall(r'(\d+(?:\.\d+)?)\s*h(?:ours?|r)?(?:\s|,|\.)', text, re.IGNORECASE)
    if not time_pat:
        time_pat = re.findall(r'(overnight|16\s*h|12\s*h)', text, re.IGNORECASE)
    if time_pat:
        conditions["time"] = time_pat[0]

    # Common solvents
    solvents = ["DCM", "THF", "DMF", "DMSO", "MeOH", "EtOH", "acetone",
                "toluene", "EtOAc", "hexane", "acetonitrile", "dioxane",
                "CH2Cl2", "chloroform", "AcOH", "water"]
    found_solvents = [s for s in solvents if re.search(r'\b' + s + r'\b', text, re.IGNORECASE)]
    if found_solvents:
        conditions["solvents"] = list(set(found_solvents))

    # Common reagents/reactions
    reactions = {
        "Suzuki": r'Suzuki|Pd.*boronic|boronic.*Pd',
        "Buchwald-Hartwig": r'Buchwald|Hartwig|C-N coupling',
        "Grignard": r'Grignard|RMgX|MgBr',
        "Mitsunobu": r'Mitsunobu|DIAD|DEAD.*PPh3',
        "Swern": r'Swern|oxalyl chloride.*DMSO',
        "Wittig": r'Wittig|ylide|Ph3P=',
        "Hydrogenation": r'H2.*Pd|Pd.*H2|hydrogenation|Pd/C',
        "Reductive amination": r'reductive amination|NaBH.*amine|NaBH3CN',
        "Amide coupling": r'amide.*coupling|EDC|HATU|HBTU|peptide coupling',
        "Fmoc SPPS": r'Fmoc|solid.phase|SPPS|resin',
    }
    detected = [name for name, pat in reactions.items() if re.search(pat, text, re.IGNORECASE)]
    if detected:
        conditions["reaction_types"] = detected

    return conditions


async def get_pubchem_synthesis_refs(cid: int) -> list:
    """Get synthesis-related literature from PubChem for a compound CID."""
    refs = []
    try:
        async with httpx.AsyncClient(timeout=10) as client:
            # Literature references
            url = f"{PUBCHEM}/compound/cid/{cid}/xrefs/PMID/JSON"
            r = await client.get(url)
            if r.status_code == 200:
                pmids = r.json().get("InformationList", {}).get("Information", [{}])[0].get("PMID", [])
                for pmid in pmids[:5]:
                    refs.append({
                        "source": "PubChem/PubMed",
                        "pmid": pmid,
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    })
    except Exception:
        pass
    return refs


async def check_commercial_availability(smiles: str, name: str = "") -> dict:
    """
    Check if compound is commercially available via:
    - PubChem vendor data
    - Enamine catalog API
    """
    availability = {
        "available": False,
        "vendors": [],
        "note": ""
    }

    try:
        async with httpx.AsyncClient(timeout=10) as client:
            # Check Enamine REAL database via their API
            if smiles:
                enamine_url = "https://enaminestore.com/api/v1/search"
                params = {"smiles": smiles, "mode": "exact", "currency": "USD"}
                try:
                    r = await client.get(enamine_url, params=params, timeout=8)
                    if r.status_code == 200:
                        data = r.json()
                        if data.get("found"):
                            availability["available"] = True
                            availability["vendors"].append({
                                "name": "Enamine",
                                "catalog_id": data.get("id", ""),
                                "price_usd": data.get("price", ""),
                                "delivery": "2-3 weeks",
                                "url": f"https://enaminestore.com/catalog/{data.get('id', '')}"
                            })
                except Exception:
                    pass

            # Check PubChem for vendor info
            cid_url = f"{PUBCHEM}/compound/smiles/{smiles}/cids/JSON"
            r2 = await client.get(cid_url)
            if r2.status_code == 200:
                cids = r2.json().get("IdentifierList", {}).get("CID", [])
                if cids:
                    vendor_url = f"{PUBCHEM}/compound/cid/{cids[0]}/xrefs/SourceName/JSON"
                    r3 = await client.get(vendor_url)
                    if r3.status_code == 200:
                        sources = r3.json().get("InformationList", {}).get("Information", [{}])[0].get("SourceName", [])
                        chem_vendors = [s for s in sources if any(
                            v in s for v in ["Sigma", "Aldrich", "TCI", "Acros", "Alfa", "Combi", "Matrix"]
                        )]
                        if chem_vendors:
                            availability["available"] = True
                            for v in chem_vendors[:3]:
                                availability["vendors"].append({"name": v, "source": "PubChem"})

    except Exception:
        pass

    return availability
