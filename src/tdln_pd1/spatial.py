"""Auditable TDLN zone representation and cell-to-zone mapping."""

from __future__ import annotations
from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True)
class Zone:
    name: str
    role: str
    resident_cells: tuple[str, ...]


ZONES = (
    Zone("afferent_lymph", "tumour antigen and migratory DC entry", ("migratory DC",)),
    Zone("subcapsular_sinus", "lymph-borne antigen capture", ("SSM", "B cell")),
    Zone("paracortex", "DC–T-cell priming and clonal expansion", ("naive T", "Tpex", "DC", "Treg")),
    Zone("follicle", "B-cell and Tfh interactions", ("B cell", "Tfh", "FDC")),
    Zone("medulla", "plasma-cell niche and efferent exit", ("plasma cell", "macrophage")),
    Zone("efferent_lymph", "activated lymphocyte export", ("effector T",)),
)


def zone_table() -> pd.DataFrame:
    return pd.DataFrame({"zone": z.name, "role": z.role, "resident_cells": ";".join(z.resident_cells)} for z in ZONES)


def assign_zone(cell_type: str) -> str:
    value = cell_type.lower()
    rules = (
        (("dendritic", " dc"), "paracortex"), (("treg", "cd8", "t cell", "tpex"), "paracortex"),
        (("b cell", "tfh", "follic"), "follicle"), (("plasma", "macrophage"), "medulla"),
    )
    for terms, zone in rules:
        if any(term in f" {value}" for term in terms):
            return zone
    return "unresolved"
