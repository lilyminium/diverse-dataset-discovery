import importlib
import json

_SMARTS_PATH = importlib.resources.files("discoset.data") / "smarts.json"

with open(_SMARTS_PATH, "r") as file:
    SMARTS = json.load(file)
