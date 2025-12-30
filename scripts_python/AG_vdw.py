import math

ASA_TYPED_FILE = "input_with_atom_types.asa"
CHAIN_A = "A"
CHAIN_B = "B"

R_CUT = 8.0

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

asa_atoms = []
with open(ASA_TYPED_FILE) as f:
	for raw in f:
		raw = raw.strip()
		parts = raw.split()

		atom = parts[2]
		resname = parts[3]
		chain = parts[4]
		resseq = parts[5]
		x, y, z = map(float, parts[6:9])

		if len(parts) > 9:
			asa = float(parts[9])

		atype = parts[-1]
		asa_atoms.append({
			"atom": atom, "resname": resname, "chain": chain, "resseq": resseq,
			"x": x, "y": y, "z": z, "asa": asa, "type": atype
		})

def split_by_chain(atoms, chain_id):
	chain_id = chain_id.upper()
	return [a for a in atoms if a["chain"] == chain_id]

aA = split_by_chain(asa_atoms, CHAIN_A)
aB = split_by_chain(asa_atoms, CHAIN_B)

E = 0.0
r_cut2 = R_CUT * R_CUT

for ai in aA:
	pi = lj_params.get(ai["type"])
	if not pi:
		continue

	xi, yi, zi = ai["x"], ai["y"], ai["z"]

	for bj in aB:
		pj = lj_params.get(bj["type"])
		if not pj:
			continue

		dx = xi - bj["x"]
		dy = yi - bj["y"]
		dz = zi - bj["z"]
		d2 = dx*dx + dy*dy + dz*dz
		if d2 >= r_cut2:
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
		E += 4.0 * eps * (t12 - t6)

print(E)
