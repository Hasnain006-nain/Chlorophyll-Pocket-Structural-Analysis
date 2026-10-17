#!/usr/bin/env python3
"""
Structural and Microenvironment Analysis of A-1 Chlorophyll Binding Pockets
in Far-Red Light Adapted Photosystem I (PSI) Complexes (PDB: 9I9L, Replicate: 9EYS).

This script performs quantitative, reproducible structural measurements at the C2
substituent position of the pseudo-symmetric A-1A and A-1B chlorophyll pairs:
  1. C2 Substituent verification (Chl f C2-formyl vs Chl a C2-methyl)
  2. Potential hydrogen-bonding partners to the formyl oxygen (OMB)
  3. Steric clearance via real formyl oxygen (A-1B) and idealized planar formyl probes (A-1A)
  4. Atomic packing density (heavy atoms within 5.0 Å of C2)
  5. Distance to nearest aromatic ring centroid (e.g., Phe459, Phe708) and nearest Cys sulfur (Cys463)
  6. Approximate net formal charge within 8.0 Å
  7. Sequence/structural divergence against white-light reference homologues (e.g., 1JB0)

Usage:
    python Scripts/a1_pocket_comparison.py --cif Structures/9i9l.cif
    python Scripts/a1_pocket_comparison.py --cif Structures/9eys.cif
"""

import argparse
import csv
from contextlib import redirect_stdout
from datetime import datetime, timezone
from pathlib import Path
import re
import shutil
import sys
import traceback
from urllib.error import HTTPError, URLError
from urllib.request import urlopen
import warnings

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import NeighborSearch
from Bio.PDB.MMCIFParser import MMCIFParser
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RING_N = ("NA", "NB", "NC", "ND")
TAIL = re.compile(r"^C\d{1,2}$")  # phytyl-tail carbons
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
CHG_ATOMS = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "LYS": ("NZ",),
    "ARG": ("NE", "NH1", "NH2"),
}
RING_ATOMS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
}
O_ANGLE = np.deg2rad(57.0)  # C(ring)-C(formyl)=O is 123 deg -> O leaves 57 deg off C2->C axis


# ------------------------------------------------------------------ Geometric Helpers
def xyz(atom):
    return np.asarray(atom.coord, float)


def dist(p, q):
    return float(np.linalg.norm(np.asarray(p) - np.asarray(q)))


def lab(atom):
    res = atom.get_parent()
    return f"{res.get_parent().id}:{res.get_resname()}{res.id[1]}:{atom.get_id()}"


def rlab(res):
    return f"{res.get_parent().id}:{res.get_resname()}{res.id[1]}"


def is_chl(res):
    return "MG" in res and sum(a.element == "N" for a in res) >= 4


def axial_ligand(neighbor_search, res, radius=3.5):
    mg = xyz(res["MG"])
    own = set(res.get_atoms())
    best = None
    for atom in neighbor_search.search(mg, radius, "A"):
        if atom in own or atom.element not in ("N", "O", "S"):
            continue
        d = dist(xyz(atom), mg)
        if best is None or d < best[1]:
            best = (atom, d)
    return best


def find_formyl(res):
    """Identify formyl group: aldehyde carbon has exactly 2 heavy neighbours (1 O <= 1.45 A, 1 C <= 1.70 A)."""
    heavy = [a for a in res if a.element != "H"]
    hits = []
    for c in (a for a in heavy if a.element == "C"):
        nb = [b for b in heavy if b is not c and dist(xyz(b), xyz(c)) <= 1.70]
        os_ = [b for b in nb if b.element == "O" and dist(xyz(b), xyz(c)) <= 1.45]
        cs = [b for b in nb if b.element == "C"]
        if len(nb) == 2 and len(os_) == 1 and len(cs) == 1:
            hits.append((os_[0], c))
    return hits[0] if len(hits) == 1 else None


def c2_substituent(res):
    f = find_formyl(res)
    if f:
        return f[1], f[0]
    return (res["CMB"], None) if "CMB" in res else (None, None)


def attach_atom(res, c):
    cand = [a for a in res if a is not c and a.element == "C" and not TAIL.match(a.get_id())
            and dist(xyz(a), xyz(c)) < 1.75]
    return min(cand, key=lambda a: dist(xyz(a), xyz(c))) if cand else None


def ring_normal(res):
    pts = np.array([xyz(res[n]) for n in RING_N if n in res])
    if len(pts) < 4:
        return None
    return np.linalg.svd(pts - pts.mean(0))[2][-1]


def planar_probes(c, att, normal):
    """Two in-plane positions where a formyl O could sit on the C2 substituent carbon."""
    u = xyz(c) - xyz(att)
    u /= np.linalg.norm(u)
    w = np.cross(normal, u)
    w /= np.linalg.norm(w)
    return [xyz(c) + 1.22 * (np.cos(O_ANGLE) * u + s * np.sin(O_ANGLE) * w) for s in (1, -1)]


def clearance(neighbor_search, pt, own, radius=4.0):
    best = (radius, None)
    for a in neighbor_search.search(np.asarray(pt, float), radius, "A"):
        if a in own:
            continue
        d = dist(xyz(a), pt)
        if d < best[0]:
            best = (d, a)
    return best


def contact_type(atom, d):
    if atom is None:
        return "free"
    if atom.element in ("N", "O") and 2.5 <= d <= 3.3:
        return "H-bond"
    if d < 2.5 or (atom.element in ("C", "S") and d < 3.0):
        return "CLASH"
    return "vdW"


def polar_partners(neighbor_search, pt, own, radius=3.3):
    return sorted(((dist(xyz(a), pt), lab(a)) for a in neighbor_search.search(np.asarray(pt, float), radius, "A")
                   if a not in own and a.element in ("N", "O")))


# ------------------------------------------------------------------ Per-Site Metrics
def site_metrics(neighbor_search, res):
    own = set(res.get_atoms())
    c, o = c2_substituent(res)
    ax = axial_ligand(neighbor_search, res)
    m = dict(
        res=res,
        where=f"{res.get_parent().id}:{res.id[1]}",
        name=res.get_resname(),
        formyl=o is not None,
        axial=f"{lab(ax[0])} ({ax[1]:.2f} A)" if ax else "none within 3.5 A",
    )
    if c is None:
        m["error"] = "no C2 substituent atom found (expected CMB or a formyl C)"
        return m

    cc = xyz(c)
    pocket = {}
    for pt in [cc] + ([xyz(o)] if o is not None else []):
        for a in neighbor_search.search(pt, 4.5, "A"):
            r = a.get_parent()
            if a not in own and r.id[0] == " ":
                pocket[r] = (r.get_parent().id, r.id[1], r.get_resname())
    m["pocket"] = pocket
    m["packing"] = sum(1 for a in neighbor_search.search(cc, 5.0, "A") if a not in own)

    att, n = attach_atom(res, c), ring_normal(res)
    m["attach"] = att.get_id() if att is not None else "?"
    m["c2name"] = c.get_id()
    probes = []
    if att is not None and n is not None:
        for p in planar_probes(c, att, n):
            d, a = clearance(neighbor_search, p, own)
            probes.append(dict(d=d, kind=contact_type(a, d), who=lab(a) if a else "-", p=p,
                               partners=polar_partners(neighbor_search, p, own)))
    m["probes"] = probes
    m["best_probe"] = max((q["d"] for q in probes), default=np.nan)

    if o is not None:
        d, a = clearance(neighbor_search, xyz(o), own)
        m["real_o"] = dict(d=d, kind=contact_type(a, d), who=lab(a) if a else "-")
        m["hbond_partners"] = polar_partners(neighbor_search, xyz(o), own)
        m["real_o_name"] = o.get_id()
    else:
        m["hbond_partners"] = None

    near = {a.get_parent() for a in neighbor_search.search(cc, 9.0, "A")}
    ar = None
    for r in near:
        if r.id[0] == " " and r.get_resname() in RING_ATOMS:
            pts = [xyz(r[x]) for x in RING_ATOMS[r.get_resname()] if x in r]
            if len(pts) >= 5:
                d = dist(np.mean(pts, 0), cc)
                if ar is None or d < ar[0]:
                    ar = (d, rlab(r))
    m["aromatic"] = ar

    cys = [(dist(xyz(r["SG"]), cc), rlab(r)) for r in near
           if r.id[0] == " " and r.get_resname() == "CYS" and "SG" in r]
    m["cys"] = min(cys) if cys else None

    chg = {}
    for a in neighbor_search.search(cc, 8.0, "A"):
        r = a.get_parent()
        if r.id[0] == " " and a.get_id() in CHG_ATOMS.get(r.get_resname(), ()):
            chg[r] = CHARGE[r.get_resname()]
    m["charge8"] = sum(chg.values())
    m["charged_list"] = sorted(rlab(r) for r in chg)
    return m


# ------------------------------------------------------------------ Locating A-1 Sites
def find_pairs(model, neighbor_search, max_a0_dist):
    chls = [r for ch in model for r in ch if is_chl(r)]
    info = {}
    for r in chls:
        ax = axial_ligand(neighbor_search, r)
        info[r] = dict(
            mg=xyz(r["MG"]),
            f=find_formyl(r) is not None,
            lig=f"{ax[0].get_parent().get_resname()}{ax[0].get_parent().id[1]}:{ax[0].get_id()}" if ax else "none"
        )
    a0s = [r for r in chls if info[r]["lig"].startswith("MET")]
    cand = {}
    for a0 in a0s:
        others = [x for x in chls if x is not a0 and not info[x]["lig"].startswith("MET")]
        if not others:
            continue
        x = min(others, key=lambda y: dist(info[y]["mg"], info[a0]["mg"]))
        d = dist(info[x]["mg"], info[a0]["mg"])
        if d <= max_a0_dist:
            cand[x] = (a0, d)

    print(f"{len(chls)} chlorophylls, {len(a0s)} Met-ligated (A0 candidates), {len(cand)} A-1 candidates\n")
    print(f"{'A-1 candidate':<16}{'residue':<8}{'formyl?':<9}{'axial':<14}{'nearest A0':<14}{'Mg-Mg':>7}")
    for x, (a0, d) in sorted(cand.items(), key=lambda kv: (kv[0].get_parent().id, kv[0].id[1])):
        print(f"{x.get_parent().id + ':' + str(x.id[1]):<16}{x.get_resname():<8}{str(info[x]['f']):<9}"
              f"{info[x]['lig']:<14}{a0.get_parent().id + ':' + str(a0.id[1]):<14}{d:7.1f}")
    fs = [x for x in cand if info[x]["f"]]
    nonf = [x for x in cand if not info[x]["f"]]
    pairs = []
    for b in fs:
        if nonf:
            a = min(nonf, key=lambda y: dist(info[y]["mg"], info[b]["mg"]))
            pairs.append((a, b))
    return pairs, cand, info


def get_chl(model, spec):
    cid, num = spec.split(":")
    for r in model[cid]:
        if r.id[1] == int(num) and is_chl(r):
            return r
    sys.exit(f"No chlorophyll {spec}")


# ------------------------------------------------------------------ White-Light Comparison
def std_chain(chain):
    rs = [r for r in chain if r.id[0] == " " and r.get_resname() in AA3]
    return rs, "".join(AA3[r.get_resname()] for r in rs)


def align_map(a, b):
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.end_gap_score = -10.0, -0.5, 0.0
    aln = al.align(a, b)[0]
    mp = {}
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0):
            mp[a0 + k] = b0 + k
    ident = sum(a[i] == b[j] for i, j in mp.items()) / max(len(a), 1)
    return mp, ident


def wl_report(model, wl_model, sites, min_len):
    wl = {ch.id: std_chain(ch) for ch in wl_model}
    wl = {k: v for k, v in wl.items() if len(v[0]) >= min_len}
    if not wl:
        print(" (no white-light chain >= --min-len residues found)")
        return {}
    cache, counts = {}, {}
    for title, m in sites:
        diffs, rows = 0, []
        for r, (cid, num, nm) in sorted(m["pocket"].items(), key=lambda kv: (kv[1][0], kv[1][1])):
            chain = model[cid]
            if cid not in cache:
                frres, frseq = std_chain(chain)
                best = None
                for wid, (wres, wseq) in wl.items():
                    mp, ident = align_map(frseq, wseq)
                    if best is None or ident > best[0]:
                        best = (ident, wid, mp, frres, wres)
                cache[cid] = best
            ident, wid, mp, frres, wres = cache[cid]
            if r not in frres:
                continue
            i = frres.index(r)
            if i in mp:
                w = wres[mp[i]]
                same = w.get_resname() == r.get_resname()
                diffs += (not same)
                rows.append(f" {cid}:{nm}{num:<5} -> WL {wid}:{w.get_resname()}{w.id[1]:<5} "
                            f"{'same' if same else 'DIFFERENT'} (chain identity {ident:.0%})")
            else:
                rows.append(f" {cid}:{nm}{num:<5} -> WL: gap")
        print(f" {title}: {diffs} of {len(m['pocket'])} pocket residues differ from the WL homologue")
        print("\n".join(rows))
        counts[title] = diffs
    return counts


# ------------------------------------------------------------------ Reporting & Output
def fmt_partners(p):
    return "n/a (methyl)" if p is None else (", ".join(f"{l} {d:.2f}" for d, l in p) or "none <= 3.3 A")


def report(tag, m):
    print(f"\n--- {tag}: {m['where']} {m['name']} | C2 substituent: {'FORMYL (Chl f)' if m['formyl'] else 'methyl (Chl a)'}")
    if "error" in m:
        print(" ", m["error"])
        return
    print(f" C2 substituent atom : {m['c2name']} bonded to ring atom {m['attach']}")
    print(f" axial ligand to Mg : {m['axial']}")
    print(f" formyl-O H-bond partners : {fmt_partners(m['hbond_partners'])}")
    if "real_o" in m:
        r = m["real_o"]
        print(f" real formyl O environment : nearest atom {r['who']} at {r['d']:.2f} A ({r['kind']})")
    if "real_o" in m and m["probes"]:
        off = min(dist(q["p"], xyz(m["res"][m["real_o_name"]])) for q in m["probes"])
        print(f" probe validation : nearest probe is {off:.2f} A from real formyl O (expected < ~0.5)")
    for i, q in enumerate(m["probes"], 1):
        part = ", ".join(f"{l} {d:.2f}" for d, l in q["partners"]) or "none"
        print(f" planar formyl probe {i} : clearance {q['d']:.2f} A ({q['kind']}, nearest {q['who']}); polar partners: {part}")
    print(f" heavy atoms within 5 A of C2: {m['packing']}")
    print(f" residues within 4.5 A of C2 : {len(m['pocket'])} -> " +
          " ".join(f"{n}{k}" for (_, k, n) in sorted(m['pocket'].values(), key=lambda t: t[1])))
    print(f" nearest aromatic ring : {('%s at %.2f A' % (m['aromatic'][1], m['aromatic'][0])) if m['aromatic'] else 'none <= 9 A'}")
    print(f" nearest Cys S : {('%s at %.2f A' % (m['cys'][1], m['cys'][0])) if m['cys'] else 'none <= 9 A'}")
    print(f" net charge of residues <= 8 A: {m['charge8']:+d} {m['charged_list']}")


warnings.filterwarnings('ignore', message='Mean of empty slice')


class Tee:
    def __init__(self, *targets):
        self.targets = targets

    def write(self, message):
        for target in self.targets:
            target.write(message)
        return len(message)

    def flush(self):
        for target in self.targets:
            target.flush()


def make_folders(root):
    root = Path(root).expanduser()
    folders = {name: root / name for name in
               ("Structures", "Figures", "Reports", "Tables", "Logs", "Scripts")}
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)
    return folders


def ensure_structure(filename, folder):
    requested = Path(filename).expanduser()
    destination = folder / requested.name
    if requested.is_file():
        if requested.resolve() != destination.resolve():
            shutil.copy2(requested, destination)
        return destination
    if destination.is_file():
        return destination

    match = re.fullmatch(r"([A-Za-z0-9]{4})\.cif", requested.name, flags=re.IGNORECASE)
    if match:
        pdb_id = match.group(1).upper()
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        temporary = destination.with_name(destination.name + ".download")
        print(f"Downloading {pdb_id} from RCSB into {destination} ...")
        try:
            with urlopen(url, timeout=90) as source, temporary.open("wb") as target:
                shutil.copyfileobj(source, target)
            if temporary.stat().st_size == 0:
                raise RuntimeError(f"Downloaded CIF for {pdb_id} was empty")
            temporary.replace(destination)
        except (HTTPError, URLError, OSError, RuntimeError) as exc:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(
                f"Could not download {pdb_id}. Put {requested.name} in "
                f"{folder} and retry. Original error: {exc}"
            ) from exc
        return destination
    raise FileNotFoundError(
        f"Structure file '{filename}' not found. Upload it into {folder}, "
        "or use --cif with an existing file path."
    )


def save_table(path, results, pdb_id, wl_counts):
    columns = (
        "structure", "replicate", "site", "residue", "pigment", "c2_formyl_present",
        "axial_ligand", "formyl_hbond_partners", "best_formyl_probe_clearance_A",
        "real_formyl_O_nearest_distance_A", "real_formyl_O_contact_type",
        "heavy_atoms_within_5A", "pocket_residues_within_4_5A",
        "nearest_aromatic_residue", "nearest_aromatic_distance_A",
        "nearest_cysteine_residue", "nearest_cysteine_S_distance_A",
        "charge_within_8A", "nearby_charged_residues",
        "white_light_different_pocket_residues", "error",
    )
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=columns)
        writer.writeheader()
        for i, (site_a, site_b) in enumerate(results, 1):
            wl = wl_counts[i - 1] if i <= len(wl_counts) else {}
            for site_name, m in (("A-1A", site_a), ("A-1B", site_b)):
                aromatic, cys = m.get("aromatic"), m.get("cys")
                real_o = m.get("real_o", {})
                pocket = m.get("pocket", {})
                hb = m.get("hbond_partners")
                row = {
                    "structure": pdb_id,
                    "replicate": i,
                    "site": site_name,
                    "residue": m.get("where", ""),
                    "pigment": m.get("name", ""),
                    "c2_formyl_present": m.get("formyl", False),
                    "axial_ligand": m.get("axial", ""),
                    "formyl_hbond_partners": fmt_partners(hb) if "hbond_partners" in m else "",
                    "best_formyl_probe_clearance_A": m.get("best_probe", ""),
                    "real_formyl_O_nearest_distance_A": real_o.get("d", ""),
                    "real_formyl_O_contact_type": real_o.get("kind", ""),
                    "heavy_atoms_within_5A": m.get("packing", ""),
                    "pocket_residues_within_4_5A": "; ".join(
                        f"{cid}:{nm}{n}" for cid, n, nm in
                        sorted(pocket.values(), key=lambda v: (v[0], v[1]))
                    ),
                    "nearest_aromatic_residue": aromatic[1] if aromatic else "",
                    "nearest_aromatic_distance_A": aromatic[0] if aromatic else "",
                    "nearest_cysteine_residue": cys[1] if cys else "",
                    "nearest_cysteine_S_distance_A": cys[0] if cys else "",
                    "charge_within_8A": m.get("charge8", ""),
                    "nearby_charged_residues": "; ".join(m.get("charged_list", [])),
                    "white_light_different_pocket_residues": wl.get(site_name, ""),
                    "error": m.get("error", ""),
                }
                writer.writerow(row)


def analyze(args, input_cif, wl_file, png_path, pdf_path, csv_path):
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("FR", str(input_cif))[0]
    neighbor_search = NeighborSearch(list(model.get_atoms()))

    if args.site_a and args.site_b:
        pairs = [(get_chl(model, args.site_a), get_chl(model, args.site_b))]
    else:
        pairs, cand, info = find_pairs(model, neighbor_search, args.max_a0_dist)
        if not pairs:
            sys.exit("\nNo formyl-bearing A-1 candidate found.")

    res_all, wl_counts = [], []
    wl_model = parser.get_structure("WL", str(wl_file))[0] if wl_file else None
    for k, (a, b) in enumerate(pairs, 1):
        print(f"\n================ replicate {k} ================")
        ma, mb = site_metrics(neighbor_search, a), site_metrics(neighbor_search, b)
        report("A-1A", ma)
        report("A-1B", mb)
        res_all.append((ma, mb))
        if wl_model is not None:
            print("\n White-light comparison:")
            valid_sites = [(name, metrics) for name, metrics in (("A-1A", ma), ("A-1B", mb))
                           if "pocket" in metrics]
            wl_counts.append(wl_report(model, wl_model, valid_sites, args.min_len))

    keys = [("best_probe", "O clearance: A simulated / B real (Å)"),
            ("packing", "Other heavy atoms within 5 Å of CMB"),
            ("arom", "Nearest aromatic ring (Å)"),
            ("charge8", "Approx. charge within 8 Å")]

    def val(m, k):
        if "error" in m:
            return np.nan
        if k == "arom":
            return m["aromatic"][0] if m["aromatic"] else np.nan
        return m[k]

    print("\n================ SUMMARY (mean over replicates) ================")
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.6), dpi=300)
    for ax, (k, title) in zip(axs, keys):
        va = np.array([val(a, k) for a, _ in res_all], float)
        vb = np.array([val(b, k) for _, b in res_all], float)
        print(f"{title:<42} A-1A {np.nanmean(va):7.2f} | A-1B {np.nanmean(vb):7.2f} | n={len(va)}")
        ax.bar([0, 1], [np.nanmean(va), np.nanmean(vb)], color=["#3b6ea5", "#d1832b"], width=0.6)
        for j, vv in enumerate((va, vb)):
            ax.scatter(np.full(len(vv), j), vv, color="k", s=10, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["A-1A\n(Chl a)", "A-1B\n(Chl f)"])
        ax.set_title(title, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fb = [bool(b["hbond_partners"]) for _, b in res_all if b.get("hbond_partners") is not None]
    print(f"A-1B formyl has an H-bond partner in {sum(fb)} of {len(fb)} replicates (expected 0)")
    fig.suptitle(f"Corrected A-1 chlorophyll pocket comparison — {input_cif.stem.upper()}", fontsize=10)
    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print("saved", png_path)
    print("saved", pdf_path)
    save_table(csv_path, res_all, input_cif.stem.upper(), wl_counts)
    print("saved", csv_path)


def main():
    repo_root = Path(__file__).resolve().parent.parent
    default_root = repo_root if (repo_root / "Structures").is_dir() else Path("/content/drive/MyDrive/Rabia..Australia")

    ap = argparse.ArgumentParser(description="A-1 chlorophyll pocket comparison analysis")
    ap.add_argument("--cif", default="9i9l.cif", help="PDB mmCIF path or filename (default: 9i9l.cif)")
    ap.add_argument("--site-a")
    ap.add_argument("--site-b")
    ap.add_argument("--max-a0-dist", type=float, default=14.0)
    ap.add_argument("--wl", help="optional white-light PSI mmCIF, e.g. 1jb0.cif")
    ap.add_argument("--min-len", type=int, default=600)
    ap.add_argument("--out", help="custom PNG filename")
    ap.add_argument("--output-root", default=str(default_root), help="Root directory for outputs")
    args = ap.parse_args()

    folders = make_folders(Path(args.output_root))
    input_cif = ensure_structure(args.cif, folders["Structures"])
    wl_file = ensure_structure(args.wl, folders["Structures"]) if args.wl else None
    tag = input_cif.stem.lower()
    if wl_file is not None:
        tag += f"_wl_{wl_file.stem.lower()}"
    if args.site_a and args.site_b:
        tag += f"_{args.site_a.replace(':', '-')}_{args.site_b.replace(':', '-')}"
    image_name = Path(args.out).name if args.out else f"a1_pocket_comparison_{tag}.png"
    png_path = folders["Figures"] / image_name
    pdf_path = png_path.with_suffix(".pdf")
    report_path = folders["Reports"] / f"a1_pocket_comparison_{tag}.txt"
    table_path = folders["Tables"] / f"a1_pocket_comparison_{tag}.csv"
    log_path = folders["Logs"] / f"a1_pocket_comparison_{tag}.log"

    start_time = datetime.now(timezone.utc)
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"Run started (UTC): {start_time.isoformat()}\n")
        log.write(f"Input: {input_cif}\n")
        log.write(f"Output root: {Path(args.output_root)}\n")
        log.flush()
        try:
            with report_path.open("w", encoding="utf-8") as report_file:
                with redirect_stdout(Tee(sys.stdout, report_file)):
                    print("A-1 chlorophyll pocket comparison")
                    print("Structure:", input_cif)
                    print("Save folder:", args.output_root)
                    analyze(args, input_cif, wl_file, png_path, pdf_path, table_path)
            log.write(f"Run finished successfully (UTC): {datetime.now(timezone.utc).isoformat()}\n")
        except BaseException:
            log.write("Run failed. Error details:\n")
            log.write(traceback.format_exc())
            log.flush()
            raise


if __name__ == "__main__":
    main()
