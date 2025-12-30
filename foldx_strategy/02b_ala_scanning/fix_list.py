# fix_list.py
# Dictionary to convert 3-letter codes to 1-letter codes
aa_map = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
    'HIS':'H', 'HIE':'H', 'HID':'H', 'HIP':'H', # Histidine variants
    'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P',
    'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

with open("interface_residues.txt", "r") as infile, open("individual_list_corrected.txt", "w") as outfile:
    for line in infile:
        parts = line.split()
        if len(parts) >= 3:
            chain = parts[0]
            resid = parts[1]
            aa_3 = parts[2]
            
            if aa_3 in aa_map:
                aa_1 = aa_map[aa_3]
                # FoldX Format: OriginalChainResIDMutant;
                # Example: SA19A; (Serine Chain A Pos 19 to Alanine)
                mutation_string = f"{aa_1}{chain}{resid}A;"
                outfile.write(mutation_string + "\n")
            else:
                print(f"Warning: Unrecognized amino acid {aa_3}")

print("Done! File 'individual_list_corrected.txt' created.")
