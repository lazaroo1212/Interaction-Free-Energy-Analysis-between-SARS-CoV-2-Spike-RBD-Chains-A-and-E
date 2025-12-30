import math

K = 332.06371   # kcal·Å/(mol·e^2)

atoms_A = []
atoms_E = []

with open('7EKF_cmip.pdbqt', 'r') as f:
	for line in f:
		line = line.strip().split()
		if len(line) == 13:
			x, y, z = map(float, line[6:9])   # coordinates
			q = float(line[11])              # partial charge
			chain = line[4]
			
			if chain == "A":
				atoms_A.append((x, y, z, q))
			elif chain == "B":
				atoms_E.append((x, y, z, q))
			
def distance(a, b):
	return math.sqrt(
		(a[0] - b[0])**2 +
		(a[1] - b[1])**2 +
		(a[2] - b[2])**2
	)

energy = 0.0

for i in atoms_A:
	qi = i[3]
	for j in atoms_E:
		qj = j[3]
		rij = distance(i, j)
		if rij != 0.0:
			eps = 86.9525/(1-7.7839*math.exp(-0.3153*rij))-8.5525
			energy += K * (qi * qj) / (eps * rij)

print("Electrostatic energy (kcal/mol):", energy)
