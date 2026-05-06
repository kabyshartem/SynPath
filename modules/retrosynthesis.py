"""
Retrosynthesis Engine
Layer 2: RXN4Chemistry (IBM) for AI retrosynthesis tree
         USPTO reaction database for literature-backed reactions
"""
import httpx
import asyncio
from typing import Optional
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config


async def rxn4chem_retrosynthesis(smiles: str) -> Optional[dict]:
    """
    Call IBM RXN4Chemistry API for retrosynthetic analysis.
    Returns tree of routes with SMILES for each step.
    Requires free API key from rxn.app.accelerate.science
    """
    if config.RXN4CHEM_API_KEY == "your-rxn4chem-api-key-here":
        return {"error": "RXN4Chemistry API key not configured", "routes": []}

    headers = {
        "Authorization": config.RXN4CHEM_API_KEY,
        "Content-Type": "application/json"
    }

    try:
        async with httpx.AsyncClient(timeout=60) as client:
            # Step 1: Submit retrosynthesis job
            payload = {
                "smiles": smiles,
                "ai_model": "2020-07-01",
                "max_steps": config.MAX_STEPS,
            }
            r = await client.post(
                f"{config.RXN_BASE_URL}/retrosynthesis",
                json=payload,
                headers=headers
            )
            if r.status_code != 200:
                return {"error": f"RXN4Chem error: {r.status_code}", "routes": []}

            job = r.json()
            job_id = job.get("prediction_id") or job.get("id")
            if not job_id:
                return {"error": "No job ID returned", "routes": []}

            # Step 2: Poll for results
            for _ in range(30):
                await asyncio.sleep(3)
                r2 = await client.get(
                    f"{config.RXN_BASE_URL}/retrosynthesis/{job_id}",
                    headers=headers
                )
                result = r2.json()
                status = result.get("status", "")

                if status == "SUCCESS":
                    return parse_rxn4chem_result(result)
                elif status in ("FAILED", "ERROR"):
                    return {"error": "RXN4Chem job failed", "routes": []}

            return {"error": "Timeout waiting for RXN4Chem", "routes": []}

    except Exception as e:
        return {"error": str(e), "routes": []}


def parse_rxn4chem_result(result: dict) -> dict:
    """Parse RXN4Chemistry response into standardized route format."""
    routes = []
    raw_routes = result.get("retrosynthetic_paths", []) or result.get("sequences", [])

    for i, route in enumerate(raw_routes[:config.MAX_ROUTES]):
        steps = []
        for step in route.get("steps", []) or route.get("reactions", []):
            steps.append({
                "step_number": step.get("step", i + 1),
                "reactants": step.get("reactants", []),
                "product": step.get("product", ""),
                "reaction_smiles": step.get("smiles", ""),
                "confidence": step.get("confidence", 0),
                "reaction_class": step.get("reaction_class", "Unknown"),
            })

        routes.append({
            "id": i + 1,
            "source": "RXN4Chemistry",
            "steps": steps,
            "num_steps": len(steps),
            "overall_confidence": route.get("confidence", 0),
        })

    return {"routes": routes, "source": "RXN4Chemistry"}


async def search_uspto(smiles: str, max_results: int = 20) -> list:
    """
    Search USPTO reaction database via PubChem for reactions
    producing the target compound.
    """
    reactions = []
    try:
        async with httpx.AsyncClient(timeout=15) as client:
            # Search by product SMILES using PubChem
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
            r = await client.get(url)
            if r.status_code != 200:
                return []

            cids = r.json().get("IdentifierList", {}).get("CID", [])
            if not cids:
                return []

            cid = cids[0]

            # Get patent/literature data from PubChem
            lit_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/PatentID/JSON"
            r2 = await client.get(lit_url)
            if r2.status_code == 200:
                patents = r2.json().get("InformationList", {}).get("Information", [{}])[0].get("PatentID", [])
                for p in patents[:5]:
                    reactions.append({
                        "source": "USPTO",
                        "patent_id": p,
                        "url": f"https://patents.google.com/patent/{p}"
                    })

    except Exception:
        pass

    return reactions


async def search_ord(smiles: str) -> list:
    """
    Search Open Reaction Database for reactions involving target SMILES.
    ORD has 1M+ reactions from literature.
    """
    reactions = []
    try:
        async with httpx.AsyncClient(timeout=15) as client:
            # ORD search endpoint
            url = "https://client.open-reaction-database.org/api/query"
            payload = {
                "reaction": {
                    "products": [{"identifiers": [{"type": "SMILES", "value": smiles}]}]
                },
                "limit": 10
            }
            r = await client.post(url, json=payload)
            if r.status_code == 200:
                data = r.json()
                for rxn in data.get("reactions", []):
                    reactions.append({
                        "source": "ORD",
                        "reaction_id": rxn.get("reaction_id", ""),
                        "doi": rxn.get("provenance", {}).get("doi", ""),
                        "yield": rxn.get("outcomes", [{}])[0].get("yield", ""),
                        "conditions": extract_ord_conditions(rxn),
                    })
    except Exception:
        pass

    return reactions


def extract_ord_conditions(rxn: dict) -> dict:
    """Extract reaction conditions from ORD reaction object."""
    conditions = {}
    try:
        setup = rxn.get("setup", {})
        conditions["temperature"] = setup.get("temperature", {}).get("value", "")
        conditions["atmosphere"] = setup.get("atmosphere", {}).get("type", "")

        inputs = rxn.get("inputs", {})
        solvents = []
        reagents = []
        for inp in inputs.values():
            for comp in inp.get("components", []):
                role = comp.get("reaction_role", "")
                name = comp.get("identifiers", [{}])[0].get("value", "")
                if role == "SOLVENT":
                    solvents.append(name)
                elif role in ("REAGENT", "CATALYST"):
                    reagents.append(name)

        conditions["solvents"] = solvents
        conditions["reagents"] = reagents

        for outcome in rxn.get("outcomes", []):
            for product in outcome.get("products", []):
                for measure in product.get("measurements", []):
                    if measure.get("type") == "YIELD":
                        conditions["yield"] = measure.get("percentage", {}).get("value", "")
    except Exception:
        pass

    return conditions
