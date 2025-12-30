import math
import csv
import matplotlib.pyplot as plt

# --------- User inputs ---------
PDBQT_FILE = "6m0j_cmip.pdbqt"
ASA_TYPED_FILE = "input_with_atom_types.asa"          # ATOM ... + last column = type (A,C,N,OA,...)
INTERFACE_FILE = "interface_residues.txt"             # chain resseq resname
CHAIN_A = "A"
CHAIN_B = "E"

KEEP_ALA = {"N", "CA", "C", "O", "CB"}  # Ala heavy atoms

# Electrostatics
K_COUL = 332.06371  # kcal·Å/(mol·e^2)
def eps_r(r):
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

# Solvation (fsrf, you called it sigma_dict)
fsrf = {
    "A":  0.111, "C":  0.019, "N": -0.124, "P":  0.000, "SA": 0.026,
    "HO": 0.000, "HD": 0.000, "H":  0.000, "OH": -0.043, "OC": -0.069, "OA": -0.031,
}

# vdW (LJ params)
lj_params = {
    "A":  {"epsilon": 0.09, "sigma": 3.40},
    "C":  {"epsilon": 0.09, "sigma": 3.40},
    "N":  {"epsilon": 0.17, "sigma": 3.25},
    "P":  {"epsilon": 0.20, "sigma": 3.74},
    "SA": {"epsilon": 0.25, "sigma": 3.56},
    "HO": {"epsilon": 0.00, "sigma": 0.00},
    "HD": {"epsilon": 0.02, "sigma": 1.07},
    "H":  {"epsilon": 0.02, "sigma": 2.65},
    "OH": {"epsilon": 0.21, "sigma": 2.96},
    "OC": {"epsilon": 0.21, "sigma": 2.96},
    "OA": {"epsilon": 0.21, "sigma": 2.96},
}

R_ON = 9.9
R_OFF = 10.0

# ---------- Helpers ----------
def dist(a, b):
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def switching(r, r_on, r_off):
    if r <= r_on:
        return 1.0
    if r >= r_off:
        return 0.0
    num = (r_off*r_off - r*r)**2 * (r_off*r_off + 2*r*r - 3*r_on*r_on)
    den = (r_off*r_off - r_on*r_on)**3
    return num / den

def is_removed_by_mutation(atom, mut_chain, mut_resseq):
    """Return True if this atom should be removed under X->Ala mutation."""
    if mut_chain is None:
        return False
    if atom["chain"] != mut_chain or atom["resseq"] != mut_resseq:
        return False
    # Skip GLY special handling outside if needed
    return atom["atom"] not in KEEP_ALA

# ---------- Read interface residues ----------
interface = []
with open(INTERFACE_FILE) as f:
    for raw in f:
        p = raw.split()
        if not p:
            continue
        chain = p[0].upper()
        resseq = p[1]
        resname = p[2].upper()
        interface.append((chain, resseq, resname))

# ---------- Parse PDBQT for electrostatics ----------
pdbqt_atoms = []
with open(PDBQT_FILE) as f:
    for raw in f:
        parts = raw.strip().split()
        if len(parts) == 13 and parts[0] in ("ATOM"):
            atom = parts[2].upper()
            resname = parts[3].upper()
            chain = parts[4].upper()
            resseq = parts[5]
            x, y, z = map(float, parts[6:9])
            q = float(parts[11])
            pdbqt_atoms.append({
                "atom": atom, "resname": resname, "chain": chain, "resseq": resseq,
                "x": x, "y": y, "z": z, "q": q
            })

# ---------- Parse typed ASA for solv + vdW ----------
asa_atoms = []
with open(ASA_TYPED_FILE) as f:
    for raw in f:
        parts = raw.split()
        atom = parts[2].upper()
        resname = parts[3].upper()
        chain = parts[4].upper()
        resseq = parts[5]
        x, y, z = map(float, parts[6:9])
        asa = float(parts[9])
        atype = parts[-1]  # appended type
        asa_atoms.append({
            "atom": atom, "resname": resname, "chain": chain, "resseq": resseq,
            "x": x, "y": y, "z": z, "asa": asa, "type": atype
        })

# Split lists by chain for speed
def split_by_chain(atoms, chain_id):
    return [a for a in atoms if a["chain"] == chain_id]

pA = split_by_chain(pdbqt_atoms, CHAIN_A)
pB = split_by_chain(pdbqt_atoms, CHAIN_B)
aA = split_by_chain(asa_atoms, CHAIN_A)
aB = split_by_chain(asa_atoms, CHAIN_B)

# ---------- Energy terms ----------
def electrostatics(mut_chain=None, mut_resseq=None):
    E = 0.0
    for ai in pA:
        if is_removed_by_mutation(ai, mut_chain, mut_resseq):
            continue
        qi = ai["q"]
        for bj in pB:
            if is_removed_by_mutation(bj, mut_chain, mut_resseq):
                continue
            r = dist(ai, bj)
            if r == 0.0:
                continue
            eps = eps_r(r)
            if eps == 0.0:
                continue
            E += K_COUL * (qi * bj["q"]) / (eps * r)
    return E

def solv(mut_chain=None, mut_resseq=None):
    G = 0.0
    for a in asa_atoms:
        if is_removed_by_mutation(a, mut_chain, mut_resseq):
            continue
        if a["asa"] != 0.0:
            G += a["asa"] * fsrf[a["type"]]
    return G

def vdw(mut_chain=None, mut_resseq=None):
    E = 0.0
    r_off2 = R_OFF * R_OFF
    for ai in aA:
        if is_removed_by_mutation(ai, mut_chain, mut_resseq):
            continue
        pi = lj_params.get(ai["type"])
        if not pi:
            continue
        for bj in aB:
            if is_removed_by_mutation(bj, mut_chain, mut_resseq):
                continue
            pj = lj_params.get(bj["type"])
            if not pj:
                continue

            dx = ai["x"] - bj["x"]
            dy = ai["y"] - bj["y"]
            dz = ai["z"] - bj["z"]
            d2 = dx*dx + dy*dy + dz*dz
            if d2 >= r_off2:
                continue

            r = math.sqrt(d2)
            if r == 0.0:
                continue

            sig = 0.5 * (pi["sigma"] + pj["sigma"])
            eps = math.sqrt(pi["epsilon"] * pj["epsilon"])
            if eps == 0.0:
                continue

            sr = sig / r
            t6 = sr**6
            t12 = t6*t6
            e_lj = 4.0 * eps * (t12 - t6)

            if r > R_ON:
                e_lj *= switching(r, R_ON, R_OFF)

            E += e_lj
    return E

def total_deltaG(mut_chain=None, mut_resseq=None):
    return electrostatics(mut_chain, mut_resseq) + vdw(mut_chain, mut_resseq) + solv(mut_chain, mut_resseq)

# ---------- WT ----------
wt = total_deltaG()
print("WT ΔG(A–B):", wt)

# ---------- Ala scan ----------
results = []
for chain, resseq, resname in interface:
    # GLY special case: no CB to "keep". Easiest: skip or treat as same atoms kept (N,CA,C,O)
    # Here: we still remove atoms not in KEEP_ALA; GLY has no CB anyway, so this is fine.
    mut = total_deltaG(chain, resseq)
    ddg = mut - wt
    results.append((chain, resseq, resname, ddg))

# Save CSV
with open("ala_scan_results.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["chain", "resseq", "resname", "ddG_Ala_minus_WT"])
    for row in results:
        w.writerow(row)

# Plot: sort by ddG descending (most destabilizing at top)
results_sorted = sorted(results, key=lambda x: x[3], reverse=True)
labels = [f"{c}{n}{r}" for c, n, r, _ in results_sorted]
vals = [ddg for _, _, _, ddg in results_sorted]

plt.figure(figsize=(max(10, len(vals)*0.35), 5))
plt.bar(range(len(vals)), vals)
plt.axhline(0.0, linewidth=1)
plt.xticks(range(len(vals)), labels, rotation=90)
plt.ylabel("ΔΔG = ΔG(mut) − ΔG(WT)  (kcal/mol)")
plt.title("Ala scan on interface residues (A–E)")
plt.tight_layout()
plt.savefig("ala_scan_results.png", dpi=200)
print("Wrote ala_scan_results.csv and ala_scan_results.png")
