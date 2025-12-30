# render_rotation.pml
# Script with LABELS for each residue
# ---------------------------------------------------------

reinitialize

# --- 1. LOAD AND PREPARE ---
load /home/davidlazaro/Desktop/Bioinformatics3/BIOFISICA/project_biofisica/BioPhysics/Energy_analysis_project/NACCESS/NACCESS/6m0j_cmip.pdbqt, model1

# Basic visual configuration
bg_color white
set ray_opaque_background, 1
hide everything, all
show cartoon, model1

# --- 2. SELECTIONS ---
select mut_site, \
(chain E and resi 486 and resn PHE) or \
(chain E and resi 505 and resn TYR) or \
(chain E and resi 489 and resn TYR) or \
(chain A and resi 417 and resn TYR) or \
(chain A and resi 37  and resn GLU) or \
(chain A and resi 34  and resn HIE) or \
(chain F and resi 501 and resn ASN) or \
(chain A and resi 83  and resn TYR) or \
(chain E and resi 403 and resn ARG) or \
(chain E and resi 500 and resn THR) or \
(chain A and resi 356 and resn GLU) or \
(chain E and resi 496 and resn GLY) or \
(chain E and resi 475 and resn ALA) or \
(chain E and resi 446 and resn GLY) or \
(chain A and resi 38  and resn ASP) or \
(chain E and resi 453 and resn TYR)

# Environment
select environment, byres (mut_site around 5)

# --- 3. VISUAL STYLE ---
# Protein context
color gray80, model1
set cartoon_transparency, 0.40
set cartoon_smooth_loops, on

# Mutation site (Red sticks)
show sticks, mut_site
color red, mut_site
set stick_radius, 0.25

# Environment (Gray sticks)
show sticks, environment and not mut_site
color gray60, environment and not mut_site

# --- 4. LABELS (NEW SECTION) ---
# We label only the Alpha Carbon (CA) to avoid clutter
# Format: "RESN-RESI" (e.g., TYR-505)
label mut_site and name CA, "%s-%s" % (resn, resi)

# Label styling
set label_color, black
set label_size, 24       # Font size
set label_font_id, 7     # Bold font for better readability
set label_position, (0.5, 0.5, 1) # Offset label slightly from the atom

# --- 5. CAMERA ---
orient environment
zoom environment, 14

# --- 6. IMAGE GENERATION LOOP ---
python
import pymol
from pymol import cmd
import os

# Configuration
prefix = "rotation_frame_"
total_frames = 36  
degrees_per_step = 10
output_path = "/home/davidlazaro/Desktop/Bioinformatics3/BIOFISICA/project_biofisica/BioPhysics/Energy_analysis_project/NACCESS/NACCESS/"

try:
    os.chdir(output_path)
    print("Saving images to: " + output_path)
except:
    print("Saving to default folder.")

print("--- STARTING RENDER ---")

for i in range(total_frames):
    filename = "{}{:03d}.png".format(prefix, i)
    cmd.ray(1200, 1200) 
    cmd.png(filename)
    print("Saved: " + filename)
    cmd.turn("y", degrees_per_step)

print("--- FINISHED ---")
python end
