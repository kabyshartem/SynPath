"""SynPath A — Backend v1.3"""
import asyncio, uvicorn, json
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
import sys, os
sys.path.append(os.path.dirname(__file__))

from modules.input_parser import detect_input_type, name_to_smiles, smiles_to_info, parse_peptide_sequence
from modules.retrosynthesis import rxn4chem_retrosynthesis, search_ord, search_uspto
from modules.literature import search_synthesis_literature
from modules.inventory import full_availability_check, check_route_availability
from modules.ai_analysis import analyze_routes, analyze_peptide
from modules.peptide.spps import generate_spps_protocol
from modules.draw_molecule import smiles_to_png_base64
import config

app = FastAPI(title="SynPath A", version="1.3.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])
app.mount("/static", StaticFiles(directory="static"), name="static")

class AnalyzeRequest(BaseModel):
    input: str
    start_smiles: Optional[str] = ""
    scale_mmol: Optional[float] = 0.1
    c_terminus: Optional[str] = "amide"

@app.get("/")
async def index(): return FileResponse("index.html")

@app.get("/health")
async def health(): return {"status": "ok", "version": "1.3.0"}


@app.get("/mol_img")
async def mol_img(smiles: str, w: int = 98, h: int = 78):
    from fastapi.responses import Response
    import base64
    img_b64 = smiles_to_png_base64(smiles, w, h)
    if not img_b64:
        return Response(status_code=404)
    img_bytes = base64.b64decode(img_b64)
    return Response(content=img_bytes, media_type="image/png")

@app.post("/analyze")
async def analyze(req: AnalyzeRequest):
    raw = req.input.strip()
    if not raw: raise HTTPException(400, "Empty input")

    detected = detect_input_type(raw)
    input_type = detected["type"]
    compound_info, smiles, target_name = {}, "", raw

    # ── Peptide ──────────────────────────────────────────────────────────────
    if input_type == "peptide_sequence":
        fmt = detected.get("format", "single_letter")
        pdata = parse_peptide_sequence(detected["value"], fmt)
        target_name = detected["value"]
        spps = generate_spps_protocol(pdata["residues"], req.scale_mmol, "fmoc", req.c_terminus)
        pdata["spps_protocol"] = spps
        compound_info = pdata
        lit = await search_synthesis_literature(target_name)
        avail = list(await asyncio.gather(*[full_availability_check(aa["fmoc"]) for aa in pdata["residues"]]))
        ai = await analyze_peptide(target_name, spps, lit, avail)
        return {"success": True, "input_type": "peptide", "compound_info": compound_info,
                "mol_image": "", "routes": [spps], "literature": lit, "availability": avail, "ai_analysis": ai}

    # ── Organic ───────────────────────────────────────────────────────────────
    if input_type == "smiles":
        smiles = detected["value"]
        info = await smiles_to_info(smiles)
        if info: compound_info = info; target_name = info.get("iupac", smiles[:30])
    else:
        info = await name_to_smiles(raw)
        if info: compound_info = info; smiles = info.get("smiles", "")

    mol_image = smiles_to_png_base64(smiles, 400, 250) if smiles else ""

    # Run layers in parallel
    retro_result, ord_res, uspto_res, literature = await asyncio.gather(
        rxn4chem_retrosynthesis(req.start_smiles or smiles),
        search_ord(smiles),
        search_uspto(smiles),
        search_synthesis_literature(target_name, smiles),
        return_exceptions=True
    )

    retro_data = retro_result if isinstance(retro_result, dict) else {"routes": []}
    literature = literature if isinstance(literature, list) else []
    for r in (ord_res if isinstance(ord_res, list) else [])[:3]:
        literature.append({"source": "ORD", **r})

    target_avail = await full_availability_check(target_name, smiles)
    avail_results = [target_avail]

    # AI analysis — always runs regardless of RXN4Chem
    ai_result = await analyze_routes(
        target_name=target_name,
        target_smiles=smiles,
        retrosynthesis_data=retro_data,
        literature_papers=literature,
        availability_data=avail_results,
        partial_route=bool(req.start_smiles),
        start_smiles=req.start_smiles or ""
    )
    ai_result["mol_image"] = mol_image

    return {"success": True, "input_type": input_type, "compound_info": compound_info,
            "mol_image": mol_image, "routes": retro_data.get("routes", []),
            "literature": literature, "availability": avail_results, "ai_analysis": ai_result}

@app.get("/inventory")
async def get_inventory():
    from modules.inventory import load_inventory
    return load_inventory().fillna("").to_dict(orient="records")

@app.post("/inventory/check")
async def check_compound(body: dict):
    return await full_availability_check(body.get("name",""), body.get("smiles",""), body.get("cas",""))

if __name__ == "__main__":
    print("\n" + "="*50 + f"\n  SynPath A · http://127.0.0.1:{config.PORT}\n" + "="*50 + "\n")
    uvicorn.run("backend:app", host=config.HOST, port=config.PORT, reload=False)


# ── Password Protection ────────────────────────────────
import secrets as _secrets
_SITE_PASSWORD = __import__('os').environ.get("SITE_PASSWORD", "synpath2024")
_VALID_TOKENS: set = set()

class _LoginReq(__import__('pydantic').BaseModel):
    password: str

@app.post("/api/login")
async def _login(req: _LoginReq):
    if req.password == _SITE_PASSWORD:
        token = _secrets.token_hex(32)
        _VALID_TOKENS.add(token)
        return {"token": token}
    raise HTTPException(status_code=401, detail="Wrong password")

@app.get("/api/verify")  
async def _verify(x_auth_token: str = __import__('fastapi').Header(None)):
    if x_auth_token and x_auth_token in _VALID_TOKENS:
        return {"ok": True}
    raise HTTPException(status_code=401, detail="Unauthorized")
