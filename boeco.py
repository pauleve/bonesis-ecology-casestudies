import networkx as nx

from colomoto.minibn import *
import bonesis
from bonesis.asp_encoding import clingo_encode as v

def declare_insufficient_trophy(dom, a, b):
    dom.rules.append(f":- clause({v(a)},C,{v(b)},1), #count {{ L: clause({v(a)},C,L,1), L != {v(a)} }} 1")

def declare_no_symmetric_trophy(dom):
    dom.rules.append(":- clause(A,_,B,1), clause(B,_,A,1), B != A")
    dom.rules.append(":- clause(A,_,B,-1), clause(B,_,A,-1), B != A")

def declare_acyclic_trophy_network(dom):
    dom.rules.append("predates(A,B) :- clause(A,_,B,1), B != A."
                     "predates(B,A) :- clause(A,_,B,-1), B != A."
                     "predates(A,C) :- predates(A,B), predates(B,C)."
                     ":- predates(A,A)")

def keysort_influence_graph(ig):
    return (len(ig.edges()), list(sorted(list(ig.edges()))))
def keysort_result(r):
    ig, f = r
    return keysort_influence_graph(ig)

def influence_of_trophic_network(species, tn):
    ig = nx.DiGraph()
    ig.add_nodes_from(species)
    rules = []
    # we allow for self-activations
    for a in species:
        ig.add_edge(a, a, sign=1, label="+")
    # predation rules
    for (a, b, _) in tn:
        for (a,b,s,l) in [(b,a,1,"+"), (a,b,-1,"-")]:
            e = ig.edges.get((a,b))
            if e is None:
                ig.add_edge(a, b, sign=s, label=l)
            elif e["sign"] != s:
                ig.edges[(a,b)]["sign"] = 0
                ig.edges[(a,b)]["label"] = "?"
    dom = bonesis.InfluenceGraph(ig)
    for (a, b, sufficient) in tn:
        if not sufficient:
            declare_insufficient_trophy(dom, a, b)
    return dom


def parse_trophic_network(species, pg):
    edges = []
    for p in pg:
        a, t, b = p.split()
        assert a in species
        assert b in species
        assert t in ["..", "--"]
        edges.append((a, b, t == "--"))
    return edges

def trophic_of_influence_graph(ig):
    pg = nx.DiGraph()
    for a, b, d in ig.edges(data=True):
        if a == b: continue
        if d["sign"] == 1: # a feeds b
            pg.add_edge(b, a)
        elif d["sign"] == -1: # a eats b
            pg.add_edge(a,b)
    return pg

def trophic_of_boolean_network(f):
    ig = f.influence_graph()
    G = trophic_of_influence_graph(ig)
    for b in f:
        dnf = struct_of_dnf(f.ba, f[b])
        if isinstance(dnf, bool):
            continue
        for clause in dnf:
            clause = clause - {(b,True)}
            if len(clause) == 0:
                G[b][b]["sufficient"] = True
            elif len(clause) == 1:
                a, sign = list(clause)[0]
                if sign:
                    G[b][a]["sufficient"] = True
    return [(a,b,d.get("sufficient", False)) for (a,b,d) in G.edges(data=True)]

def signed_edges(g):
    return set((a, b, d["sign"]) for (a,b,d) in g.edges(data=True))

def cfg_diff(x, y):
    return {a: v for a, v in y.items() if x[a] != v}

def deviating_fasync_transitions(f, obs, cfg_of_obs):
    for x in obs.nodes():
        cx = cfg_of_obs(x)
        obs_changes = [cfg_diff(cx, cfg_of_obs(y)) for y in obs[x]]
        z = f(cx)
        dz = cfg_diff(cx, z)
        yield from ((x, a) for a, v in dz.items()
                        if {a:v} not in obs_changes)

def model_deviation(f, obs, cfg_of_obs):
    return len(list(deviating_fasync_transitions(f, obs, cfg_of_obs)))

def fasync_transitions(f, label_cfg):
    d = FullyAsynchronousDynamics(f)
    stg = d.dynamics()
    for cfg in stg.nodes:
        stg.nodes[cfg]["label"] = label_cfg(dict(zip(d.nodes, map(int, cfg))))
    return stg
