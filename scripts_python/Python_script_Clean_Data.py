cd ~/Escriptori/Bioinformatica/Third_year/BIOPHYSICS/Project/

# PASO 0: seleccionar solo cadenas A y E
python3 -m biobb_structure_checking.check_structure \
    -i 6m0j.pdb \
    -o 6m0j_chains.pdb \
    --output_format pdb \
    --non_interactive \
    --force_save \
    chains --select A,E

# PASO 1: quitar ligandos
python3 -m biobb_structure_checking.check_structure \
    -i 6m0j_chains.pdb \
    -o 6m0j_noligs.pdb \
    --output_format pdb \
    --non_interactive \
    --force_save \
    ligands --remove All

# PASO 2: quitar metales/iones
python3 -m biobb_structure_checking.check_structure \
    -i 6m0j_noligs.pdb \
    -o 6m0j_nomets.pdb \
    --output_format pdb \
    --non_interactive \
    --force_save \
    metals --remove All

# PASO 3: quitar aguas
python3 -m biobb_structure_checking.check_structure \
    -i 6m0j_nomets.pdb \
    -o 6m0j_final.pdb \
    --output_format pdb \
    --non_interactive \
    --force_save \
    water --remove Yes
    
# PASO 4: añadir H + cargas CMIP y generar PDBQT
python3 -m biobb_structure_checking.check_structure \
    -i 6m0j_final.pdb \
    -o 6m0j_cmip.pdbqt \
    --output_format pdbqt \
    --non_interactive \
    --force_save \
    add_hydrogen --add_mode auto --add_charges CMIP

