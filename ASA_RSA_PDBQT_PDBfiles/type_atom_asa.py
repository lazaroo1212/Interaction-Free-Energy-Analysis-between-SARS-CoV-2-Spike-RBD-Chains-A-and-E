aromatic_atoms = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIE": {"CG", "CD2", "CE1"},  
}


# Residues where sulfur is considered SA
sulfur_residues = {"MET", "CYS"}

input_file = "6m0j_cmip.asa"
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
                atom_type = "H"   # HC

        # Oxygen
        elif atom_name.startswith("O"):
            if atom_name.startswith("OH"):
                atom_type = "OH"
            elif atom_name.startswith("COO"):
                atom_type = "OC"
            else:
                atom_type = "OA"   # generic oxygen if needed

        # Fallback (if no rule applies)
        if atom_type is None:
            atom_type = "X"  # unknown

        # Write line with new column
        fout.write(line.rstrip() + f"  {atom_type}\n")

