"""Prespecified focus gene sets for the 10x neuron-subclass analysis.

Fixed on 2026-10-15 before any 10x whole-transcriptome result was inspected. Membership
follows official mouse gene nomenclature for each OXPHOS complex, mirroring the
MitoCarta3.0 OXPHOS subunit lists, plus a curated adenosine pathway list grouped by
function. Genes absent from the 10x reference are reported, not silently dropped.
"""
from __future__ import annotations

import re

import pandas as pd

EXACT = {
    "OXPHOS_CII": ["Sdha", "Sdhb", "Sdhc", "Sdhd"],
    "OXPHOS_CIII": ["Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrb", "Uqcrq", "Uqcrh", "Uqcr10", "Uqcr11", "Cyc1",
                    "mt-Cytb"],
    "OXPHOS_CIV": ["Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c",
                   "Cox7a1", "Cox7a2", "Cox7b", "Cox7c", "Cox8a", "Cox8b", "Ndufa4", "mt-Co1", "mt-Co2",
                   "mt-Co3"],
    "cytochrome_c": ["Cycs"],
    "mito_biogenesis": ["Ppargc1a", "Ppargc1b", "Tfam", "Nrf1", "Esrra"],
    "creatine_buffer": ["Ckb", "Ckmt1", "Slc6a8"],
    "glycolysis_lactate": ["Hk1", "Pfkm", "Pkm", "Ldha", "Ldhb", "Slc16a1", "Slc16a7", "Slc2a3"],
    "adenosine_production": ["Nt5e", "Entpd1", "Entpd2", "Entpd3", "Nt5c", "Nt5c2"],
    "adenosine_clearance": ["Adk", "Ada"],
    "adenosine_transport": ["Slc29a1", "Slc29a2", "Slc28a1", "Slc28a2", "Slc28a3"],
    "adenosine_receptor": ["Adora1", "Adora2a", "Adora2b", "Adora3"],
    "ATP_release": ["Panx1", "Gja1", "Calhm2"],
}
PATTERN = {
    "OXPHOS_CI": re.compile(r"^(Nduf(a(?!4$)(?!f)\d+|b\d+|c\d|s\d|v\d|ab1)|mt-Nd(1|2|3|4|4l|5|6))$"),
    "OXPHOS_CV": re.compile(r"^(Atp5(?!if1$)(?!s$)[a-z0-9]+|mt-Atp6|mt-Atp8)$"),
}
OXPHOS = ["OXPHOS_CI", "OXPHOS_CII", "OXPHOS_CIII", "OXPHOS_CIV", "OXPHOS_CV"]
THEME = {**{k: "mitochondria_energy" for k in OXPHOS + ["cytochrome_c", "mito_biogenesis", "creatine_buffer",
                                                         "glycolysis_lactate"]},
         **{k: "adenosine" for k in ["adenosine_production", "adenosine_clearance", "adenosine_transport",
                                     "adenosine_receptor", "ATP_release"]}}


def build(symbols, pbm_candidates=()):
    """Return a long table gene, set, theme, present for the given reference symbols."""
    symbols = pd.Index(symbols)
    rows = []
    for name, genes in EXACT.items():
        for gene in genes:
            rows.append({"gene": gene, "set": name, "theme": THEME[name], "present": gene in symbols})
    for name, pattern in PATTERN.items():
        for gene in symbols[symbols.str.match(pattern)]:
            rows.append({"gene": gene, "set": name, "theme": THEME[name], "present": True})
    for gene in pbm_candidates:
        rows.append({"gene": gene, "set": "PBM_candidate62", "theme": "PBM", "present": gene in symbols})
    return pd.DataFrame(rows).drop_duplicates(["gene", "set"])
