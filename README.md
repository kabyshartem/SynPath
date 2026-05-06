# SynPath A — AI Retrosynthesis Engine

## Quick Start (3 steps)

### 1. Install dependencies
```bash
cd synpath_a
pip install -r requirements.txt
```

### 2. Set API keys in config.py
```python
ANTHROPIC_API_KEY = "sk-ant-..."          # console.anthropic.com
RXN4CHEM_API_KEY  = "your-key-here"       # rxn.app.accelerate.science (free)
```

### 3. Run
```bash
python backend.py
```
Open browser: **http://127.0.0.1:8000**

---

## Features

### Organic Synthesis Mode
- Draw molecule in Ketcher (same engine as Reaxys)
- RXN4Chemistry retrosynthesis (IBM, free tier)
- Literature mining: Europe PMC, USPTO, ORD database
- Claude AI ranks routes and writes Org.Synth-style procedure
- Every condition referenced with real DOI

### Peptide Mode
- Enter sequence: `ACDEFGHIK` or `Ala-Gly-Phe-OH`
- Full Fmoc SPPS protocol with coupling conditions
- Automatic warning detection (aspartimide, aggregation-prone)
- Cleavage cocktail recommendation based on sequence
- Building block availability check

### Partial Route (Intermediate Available)
- Have intermediate on shelf?
- Click "HAVE INTERMEDIATE?" → paste its SMILES
- System generates route FROM your intermediate TO target

### Availability Check
- Internal inventory (inventory.csv)
- Enamine REAL catalog (real-time API)
- Sigma-Aldrich catalog
- Color coded: 🟢 in stock · 🔵 purchasable · 🔴 synthesize

---

## Customizing Inventory

Edit `inventory.csv`:
```csv
name,smiles,quantity,unit,location,cas
Fmoc-Gly-OH,O=C(OCC1...)NCC(=O)O,50,g,Freezer-A1,29022-11-5
HATU,...,100,g,Cabinet-B1,148893-10-1
```

---

## Getting API Keys

### Claude (Anthropic)
1. Go to https://console.anthropic.com
2. Create account → API Keys → Create Key
3. Paste in config.py

### RXN4Chemistry (IBM) — FREE
1. Go to https://rxn.app.accelerate.science
2. Register with email
3. Profile → API Key
4. Paste in config.py

---

## Architecture

```
User draws molecule (Ketcher)
    ↓
Layer 1: Input Parser
  - SMILES normalization via PubChem
  - Peptide sequence detection
    ↓
Layer 2: Retrosynthesis Engine
  - RXN4Chemistry (IBM ML retrosynthesis)
  - USPTO patent reactions
  - ORD (Open Reaction Database, 1M+ reactions)
    ↓
Layer 3: Literature Mining
  - Europe PMC full-text search
  - Condition extraction (temp, yield, solvent, time)
  - DOI references
    ↓
Layer 4: Claude AI Analysis
  - Ranks routes by feasibility
  - Generates Org.Synth-style procedure
  - References real DOIs
    ↓
Layer 5: Availability Check
  - Internal CSV inventory
  - Enamine REAL catalog API
  - Sigma-Aldrich API
```

---

## Roadmap

- [ ] Reaxys API integration (when available)
- [ ] ASKCOS (MIT) as fallback retrosynthesis
- [ ] Export to PDF / Word
- [ ] Batch analysis (multiple targets)
- [ ] Custom reaction templates
- [ ] Solution-phase peptide fragment condensation
