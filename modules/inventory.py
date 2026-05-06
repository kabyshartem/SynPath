"""
Inventory & Availability Layer
Layer 5: Check internal CSV + external catalogs (Enamine, Sigma)
"""
import pandas as pd
import httpx
import re
from typing import Optional
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config


def load_inventory() -> pd.DataFrame:
    """Load internal inventory from CSV."""
    try:
        df = pd.read_csv(config.INVENTORY_CSV)
        df.columns = df.columns.str.lower().str.strip()
        return df
    except Exception:
        return pd.DataFrame(columns=["name", "smiles", "quantity", "unit", "location", "cas"])


def check_internal(compound_name: str, smiles: str = "", cas: str = "") -> Optional[dict]:
    """Check if compound exists in internal inventory CSV."""
    df = load_inventory()
    if df.empty:
        return None

    # Search by name (case-insensitive partial match)
    name_matches = df[df["name"].str.lower().str.contains(compound_name.lower(), na=False)]

    # Search by CAS
    if cas and "cas" in df.columns:
        cas_matches = df[df["cas"].astype(str).str.contains(cas, na=False)]
        name_matches = pd.concat([name_matches, cas_matches]).drop_duplicates()

    if not name_matches.empty:
        row = name_matches.iloc[0]
        return {
            "found": True,
            "source": "Internal Inventory",
            "name": row.get("name", ""),
            "quantity": f"{row.get('quantity', '')} {row.get('unit', '')}".strip(),
            "location": row.get("location", ""),
            "cas": row.get("cas", ""),
        }
    return None


async def check_enamine(smiles: str) -> Optional[dict]:
    """Check Enamine REAL database for compound availability."""
    if not smiles:
        return None
    try:
        async with httpx.AsyncClient(timeout=8) as client:
            r = await client.get(
                "https://enaminestore.com/api/v1/search",
                params={"smiles": smiles, "mode": "exact", "currency": "USD"}
            )
            if r.status_code == 200:
                data = r.json()
                if data.get("found") or data.get("total", 0) > 0:
                    items = data.get("items", data.get("results", []))
                    item = items[0] if items else {}
                    return {
                        "found": True,
                        "source": "Enamine",
                        "catalog_id": item.get("id", item.get("cat_id", "")),
                        "price": item.get("price", item.get("cost", "N/A")),
                        "purity": item.get("purity", ">95%"),
                        "delivery": "2-4 weeks",
                        "url": f"https://enaminestore.com/catalog/{item.get('id', '')}"
                    }
    except Exception:
        pass
    return None


async def check_sigmaaldrich(name: str, cas: str = "") -> Optional[dict]:
    """Check Sigma-Aldrich via their product API."""
    try:
        async with httpx.AsyncClient(timeout=8) as client:
            query = cas if cas else name
            r = await client.get(
                "https://www.sigmaaldrich.com/api/products",
                params={"query": query, "page": 1, "perpage": 3},
                headers={"Accept": "application/json"}
            )
            if r.status_code == 200:
                data = r.json()
                products = data.get("results", [])
                if products:
                    p = products[0]
                    return {
                        "found": True,
                        "source": "Sigma-Aldrich",
                        "product_number": p.get("productNumber", ""),
                        "name": p.get("name", name),
                        "url": f"https://www.sigmaaldrich.com/US/en/product/sial/{p.get('productNumber', '')}"
                    }
    except Exception:
        pass
    return None


async def full_availability_check(compound_name: str, smiles: str = "", cas: str = "") -> dict:
    """
    Full availability check across all sources.
    Returns color-coded result for UI.
    """
    result = {
        "compound": compound_name,
        "smiles": smiles,
        "internal": None,
        "external": [],
        "status": "unknown",  # "in_stock" | "purchasable" | "synthesize"
        "color": "gray"       # "green" | "blue" | "red"
    }

    # 1. Check internal inventory first
    internal = check_internal(compound_name, smiles, cas)
    if internal:
        result["internal"] = internal
        result["status"] = "in_stock"
        result["color"] = "green"
        return result

    # 2. Check external vendors
    enamine = await check_enamine(smiles)
    if enamine:
        result["external"].append(enamine)

    sigma = await check_sigmaaldrich(compound_name, cas)
    if sigma:
        result["external"].append(sigma)

    if result["external"]:
        result["status"] = "purchasable"
        result["color"] = "blue"
    else:
        result["status"] = "synthesize"
        result["color"] = "red"

    return result


async def check_route_availability(route_intermediates: list) -> list:
    """
    Check availability for all intermediates in a synthesis route.
    route_intermediates: list of {name, smiles, step_number}
    Returns list with availability status for each.
    """
    results = []
    for intermediate in route_intermediates:
        avail = await full_availability_check(
            compound_name=intermediate.get("name", ""),
            smiles=intermediate.get("smiles", ""),
            cas=intermediate.get("cas", "")
        )
        avail["step_number"] = intermediate.get("step_number", 0)
        results.append(avail)
    return results
