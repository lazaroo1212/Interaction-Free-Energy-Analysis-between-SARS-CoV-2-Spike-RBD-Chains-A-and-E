import Bio.PDB
import pandas as pd
import warnings

# Suppress Biopython PDB warnings
warnings.simplefilter('ignore', Bio.PDB.PDBExceptions.PDBConstructionWarning)

def get_atom_type(atom):
    """Returns 'C', 'N', 'O', 'S' based on element"""
    return atom.element

def is_hydrophobic(residue):
    """Checks if a residue is hydrophobic"""
    hydrophobic_aas = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'PRO']
    return residue.resname in hydrophobic_aas

def is_charged(resname):
    """Returns charge type: 'POS', 'NEG', or 'NEUTRAL'"""
    if resname in ['ARG', 'LYS', 'HIS']: return 'POS'
    if resname in ['ASP', 'GLU']: return 'NEG'
    return 'NEUTRAL'

def classify_interaction(atom_a, atom_b, dist):
    """Classifies interaction based on atoms and distance"""
    res_a = atom_a.get_parent()
    res_b = atom_b.get_parent()
    elem_a = atom_a.element
    elem_b = atom_b.element
    
    # 1. SALT BRIDGES (Ionic) - Distance < 4.0 A
    # Check if residues have opposite charges
    charge_a = is_charged(res_a.resname)
    charge_b = is_charged(res_b.resname)
    
    if dist < 4.0:
        if (charge_a == 'POS' and charge_b == 'NEG') or (charge_a == 'NEG' and charge_b == 'POS'):
            # Ensure it's the side chain atoms interacting (usually N or O)
            if elem_a in ['N', 'O'] and elem_b in ['N', 'O']:
                return "Salt Bridge"

    # 2. HYDROGEN BONDS - Distance < 3.5 A
    # Interaction between N and O, or O and O, or N and N
    if dist < 3.5:
        if (elem_a in ['N', 'O'] and elem_b in ['N', 'O']):
            # Exclude Carbon backbone if needed, but keeping simple for now
            return "Hydrogen Bond"

    # 3. HYDROPHOBIC / VDW - Distance < 4.0 A
    # Carbon-Carbon interaction between hydrophobic residues
    if dist < 4.0:
        if elem_a == 'C' and elem_b == 'C':
            if is_hydrophobic(res_a) and is_hydrophobic(res_b):
                return "Hydrophobic/VdW"
            
    return None

def analyze_interface(pdb_file, chain_1_id='A', chain_2_id='E'):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure("Complex", pdb_file)
    model = structure[0]

    # Get atoms for both chains
    atoms_1 = [atom for atom in model[chain_1_id].get_atoms()]
    atoms_2 = [atom for atom in model[chain_2_id].get_atoms()]

    # Use NeighborSearch (KDTree) for fast distance calculation
    ns = Bio.PDB.NeighborSearch(atoms_2)
    
    interactions = []
    seen_pairs = set()

    print(f"Analyzing interactions between Chain {chain_1_id} and Chain {chain_2_id}...")

    for atom_a in atoms_1:
        # Find neighbors in Chain E within 4.0 Angstroms
        neighbors = ns.search(atom_a.coord, 4.0, level='A')
        
        for atom_b in neighbors:
            dist = atom_a - atom_b
            res_a = atom_a.get_parent()
            res_b = atom_b.get_parent()
            
            # Create a unique ID for the pair to avoid duplicates (e.g. A-B vs B-A)
            pair_id = tuple(sorted([f"{res_a.get_id()[1]}{res_a.resname}", f"{res_b.get_id()[1]}{res_b.resname}"]))
            
            # Classify
            itype = classify_interaction(atom_a, atom_b, dist)
            
            if itype:
                # Store Data
                interactions.append({
                    'Chain_1': chain_1_id,
                    'Residue_1': f"{res_a.resname} {res_a.get_id()[1]}",
                    'Atom_1': atom_a.name,
                    'Chain_2': chain_2_id,
                    'Residue_2': f"{res_b.resname} {res_b.get_id()[1]}",
                    'Atom_2': atom_b.name,
                    'Distance': round(dist, 2),
                    'Type': itype
                })

    return pd.DataFrame(interactions)


pdb_filename = "6m0j_Repair.pdb" 

try:
    df = analyze_interface(pdb_filename)
    

    df_clean = df.sort_values('Distance').drop_duplicates(subset=['Residue_1', 'Residue_2', 'Type'])
    

    output_filename = "Interface_Interactions_Detailed.csv"
    df_clean.to_csv(output_filename, index=False)
    
    print("\nSUCCESS!")
    print(f"Detailed interactions saved to: {output_filename}")
    print("\n--- TOP 10 STRONGEST INTERACTIONS (Shortest Distance) ---")
    print(df_clean[['Residue_1', 'Residue_2', 'Type', 'Distance']].head(10).to_string(index=False))

except FileNotFoundError:
    print(f"Error: Could not find file '{pdb_filename}'. Please check the name.")
except Exception as e:
    print(f"An error occurred: {e}")
