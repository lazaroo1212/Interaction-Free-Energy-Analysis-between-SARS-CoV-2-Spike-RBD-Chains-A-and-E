aromatic_atoms = {
	"PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
	"TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
	"TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
	"HIE": {"CG", "CD2", "CE1"},  
}


# Residues where sulfur is considered SA
sulfur_residues = {"MET", "CYS"}

input_file = "7EKF_cmip_interaction.asa"
output_file = "input_with_atom_types.asa"

with open(input_file, "r") as fin, open(output_file, "w") as fout:
	for line in fin:
		if not line.startswith("ATOM"):
			fout.write(line)
			continue

		fields = line.split()
		atom_name = fields[2]
		residue = fields[3]

		# Determine atom type
		atom_type = None

		# Carbon
		if atom_name.startswith("C"):
			if residue in aromatic_atoms and atom_name in aromatic_atoms[residue]:
				atom_type = "A"   # Aromatic C
			else:
				atom_type = "C"   # Aliphatic C

		# Nitrogen
		elif atom_name.startswith("N"):
			atom_type = "N"

		# Sulfur
		elif atom_name.startswith("S"):
			if residue in sulfur_residues:
				atom_type = "SA"

		# Phosphorus
		elif atom_name.startswith("P"):
			atom_type = "P"

		# Hydrogens
		elif atom_name.startswith("H"):
			if atom_name.startswith("HO"):
				atom_type = "HO"
			elif atom_name.startswith("HN") or atom_name.startswith("HS"):
				atom_type = "HD"
			else:
				atom_type = "H"

		# Oxygen
		elif atom_name.startswith("O"):
			# ASP/GLU carboxylate oxygens: OD1/OD2, OE1/OE2 -> OC
			if residue in {"ASP", "GLU"} and (atom_name.startswith("OD") or atom_name.startswith("OE")):
				atom_type = "OC"
			
			# hydroxyl oxygens: SER OG, THR OG1, TYR OH -> OH
			elif (residue == "SER" and atom_name == "OG") or (residue == "THR" and atom_name == "OG1") or (residue == "TYR" and atom_name == "OH"):
				atom_type = "OH"
			else:
				atom_type = "OA"         # all other oxygens
			
		# Write line with new column
		fout.write(line.rstrip() + f"  {atom_type}\n")


sigma_dict = {'A': 0.111, 'C': 0.019, 'N': -0.124, 'P': 0.000, 'SA': 0.026, 'HO': 0.000, 'HD': 0.000, 'H': 0.000, 'OH': -0.043, 'OC': -0.069, 'OA': -0.031}

AG_solv = 0

residues = []

with open('input_with_atom_types.asa', 'r') as f:
	for line in f:
		line = line.strip().split()
		if float(line[9]) != 0.0:
			AG_solv+= float(line[9])*sigma_dict[line[11]]
			if float(line[9])*sigma_dict[line[11]] > 0.1:
				if (line[4], line[5], line[3]) not in residues:
					to_add = (line[4], line[5], line[3], float(line[9])*sigma_dict[line[11]])
					residues.append(to_add)
				
	print(AG_solv)
#	print('Interface Residues:')
#	for res in residues:
#		print(f'Chain: {res[0]}, Residue Number: {res[1]} and Residue: {res[2]}, AG_solv: {res[3]}')
#with open('interface_residues.txt', 'w') as f2:
#	for res in residues:
#		f2.write(f'{res[0]}	{res[1]}	{res[2]}\n')
