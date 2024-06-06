import networkx as nx

import bonesis
from bonesis.asp_encoding import clingo_encode

from boeco import *

species = list("ABECPT")

obs = nx.DiGraph([e.split(" -> ") for e in [
    "APECT -> APEC", "APEC -> APC", "APC -> AP",
    "APET -> APE", "APE -> AP",
    "APT -> AP", "APT -> AT",
    "APC -> AP", "APC -> AC",
    "AET -> AE", "AET -> AT",
    "ECT -> EC",
    "AEC -> AC",
    "ACT -> AT",
    "AP -> A", "A -> _",
    "ET -> T",
    "AE -> A",
    "AT -> T", "AT -> A",
    "EC -> C",
    "CT -> C",
    "AC -> C", "AC -> A",
    "APCT -> APT",
    "AECT -> AET",
    "E -> _",
    "C -> _",
    "ABECT -> ABEC",
    "ABPET -> ABPE",
    "ABPEC -> ABPC",
    "ABPCT -> ABPC",
    "ABET -> ABE", "ABET -> ABT",
    "BECT -> BEC",
    "ABEC -> ABC",
    "ABCT -> ABT", "ABCT -> ABC",
    "ABPE -> ABE", "ABPE -> ABP",
    "ABPT -> ABP",
    "ABPC -> ABP",
    "BET -> BE",
    "BET -> BT",
    "ABE -> AB",
    "BEC -> BC",
    "BCT -> BC",
    "ABT -> AB",
    "ABC -> BC", "ABC -> AB",
    "ABP -> AB",
    "BE -> B",
    "BT -> B",
    "BC -> B",
    "AB -> B", "AB -> A",
    "BPECT -> BPEC",
    "BPEC -> BPC", "BPEC -> BPE",
    "BPCT -> BPC",
    "BPET -> BPE",
    "BPE -> BP",
    "BPC -> BP",
    "BPT -> BP",
    "PECT -> PEC",
    "PEC -> PC", "PEC -> PE",
    "PCT -> PC",
    "PET -> PE",
    "PC -> P",
    "PE -> P",
    "PT -> P"
]])
steady = {state for state, out in obs.out_degree() if out == 0}

def cfg_of_present(present):
    return {i: 1 if i in present else 0 for i in species}

def present_of_cfg(cfg):
    return "".join([a for a in species if cfg[a]])

def model_observations(bo, max_changes=1, monotone=True):
    # transitions (force fully asynchronous)
    with bo.scope_reachability(max_changes=max_changes, monotone=monotone):
        for x, y in obs.edges():
            ~bo.obs(cfg_of_present(x)) >= ~bo.obs(cfg_of_present(y))
    # steady states (fixed points)
    for x in steady:
        bo.fixed(~bo.obs(cfg_of_present(x)))

def minimize_deviation(bo):
    n = 0
    for x in obs.nodes():
        box = ~bo.obs(cfg_of_present(x))
        X = clingo_encode(box.name)
        Z = clingo_encode(f"_fit_{box.name}")
        fixed = set(species)
        for y in obs[x]:
            fixed.difference_update(set(x).symmetric_difference(y))
        rules = []
        for a in fixed:
            N = clingo_encode(a)
            n += 1
            rules.append(f":~ cfg({X},{N},V), eval({Z},{N},-V). [1,{n}]")
        if rules:
            rules.append(f"mcfg({Z},N,V) :- cfg({X},N,V)")
            bo.custom("\n".join(rules))

def model_deviation(f):
    return len(list(deviating_fasync_transitions(f, obs, cfg_of_present)))

dom_any = bonesis.InfluenceGraph.complete(species)

tn_english = parse_trophic_network(species, [
        "A .. B",
        "A .. E",
        "A -- P",
        "A -- C",
        "A -- T",
        "B -- B",
        "B .. P",
        "B -- C",
        "E -- C",
        "E -- T",
        "B -- T",
    ])
dom_english = influence_of_trophic_network(species, tn_english)

tn_goer = parse_trophic_network(species, [
        "A .. B",
        "A .. C",
        "A .. E",
        "A .. T",
        "A -- P",
        "B -- T",
        "B -- C",
        "E .. T",
        "P -- C",
        "P -- T",
        "C -- T",
        "T -- C",
    ])
dom_goer = influence_of_trophic_network(species, tn_goer)

tn_goer_fixed = tn_goer + parse_trophic_network(species, [
        "B .. P",
    ])
dom_goer_fixed = influence_of_trophic_network(species, tn_goer_fixed)

f_goer = bonesis.BooleanNetwork({
    "A": "A&P",
    "B": "B& !(A&!P)",
    "E": 0,
    "P": "P&!A",
    "C": 0,
    "T": "T & !((A&!P) | E | B | C | P)"
})
