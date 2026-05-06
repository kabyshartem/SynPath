import os

ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
RXN4CHEM_API_KEY  = os.environ.get("RXN4CHEM_API_KEY", "")

HOST = "0.0.0.0"
PORT = int(os.environ.get("PORT", 8000))

INVENTORY_CSV = "inventory.csv"
RXN_BASE_URL  = "https://rxn.app.accelerate.science/rxn/api/api/v1"
MAX_ROUTES    = 3
MAX_STEPS     = 5
ENAMINE_API_URL = ""
SIGMA_URL       = "https://www.sigmaaldrich.com/US/en/search#"
C_URL           = ""
