import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. DATA INPUT
# ==========================================

# A. Global Energies (Section 3.1)
global_energies = {
    'Term': ['Electrostatics', 'Van der Waals', 'Solvation', 'TOTAL'],
    'Energy (kcal/mol)': [-17.79, -36.97, 10.69, -44.07]
}

# B. Interface Residue Data (Section 3.2 - Solvation only available)
# NOTE: To calculate the full interaction formula (Elec + VdW + Solv), 
# individual Elec/VdW values are needed. Since only G_solv is available 
# in the text provided, this plot will show the Desolvation Penalty.
interface_data = [
    # Chain, ResID, AA, Solvation(G_solv)
    ('A', 19, 'SER', 1.03), ('A', 24, 'GLN', 1.20), ('A', 28, 'PHE', 0.11),
    ('A', 30, 'ASP', 1.84), ('A', 31, 'LYS', 3.22), ('A', 34, 'HIE', 2.15),
    ('A', 35, 'GLU', 1.15), ('A', 37, 'GLU', 0.43), ('A', 38, 'ASP', 1.73),
    ('A', 41, 'TYR', 0.87), ('A', 42, 'GLN', 4.90), ('A', 82, 'MET', 0.14),
    ('A', 83, 'TYR', 1.09), ('A', 353, 'LYS', 3.36), ('A', 357, 'ARG', 1.54),
    ('E', 403, 'ARG', 0.19), ('E', 417, 'LYS', 1.28), ('E', 449, 'TYR', 1.06),
    ('E', 486, 'PHE', 0.19), ('E', 489, 'TYR', 1.10), ('E', 493, 'GLN', 9.10),
    ('E', 500, 'THR', 1.52), ('E', 501, 'ASN', 0.94), ('E', 505, 'TYR', 1.10)
]

# C. Hotspot Data (Section 3.3 - Delta Delta G)
hotspot_data = [
    ('E', 486, 'PHE', 14.24),
    ('E', 505, 'TYR', 12.39),
    ('E', 489, 'TYR', 9.38),
    ('A', 41, 'TYR', 4.72),
    ('A', 37, 'GLU', 3.27),
    ('A', 34, 'HIE', 2.67),
    ('E', 501, 'ASN', 2.26),
    ('A', 83, 'TYR', 2.15),
    ('E', 403, 'ARG', 1.92),
    ('E', 500, 'THR', 1.71),
    ('A', 35, 'GLU', 1.62)
]

# Create DataFrames
df_glob = pd.DataFrame(global_energies)
df_int = pd.DataFrame(interface_data, columns=['Chain', 'ResID', 'AA', 'G_Solv'])
# Create a readable label (e.g., A-41-TYR)
df_int['Label'] = df_int['Chain'] + '-' + df_int['ResID'].astype(str) + '-' + df_int['AA']

df_hot = pd.DataFrame(hotspot_data, columns=['Chain', 'ResID', 'AA', 'DDG'])
df_hot['Label'] = df_hot['Chain'] + '-' + df_hot['ResID'].astype(str) + '-' + df_hot['AA']

# ==========================================
# 2. PLOT GENERATION
# ==========================================
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# PLOT 1: Global Energy Contributions
# Shows which forces dominate (VdW vs Elec vs Solv)
colors = ['#3498db', '#2ecc71', '#e74c3c', '#9b59b6']
axes[0].bar(df_glob['Term'], df_glob['Energy (kcal/mol)'], color=colors)
axes[0].set_title('Global Energy Contributions (Section 3.1)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Energy (kcal/mol)')
axes[0].axhline(0, color='black', linewidth=0.8)
# Add value labels on top of bars
for i, v in enumerate(df_glob['Energy (kcal/mol)']):
    axes[0].text(i, v + (1 if v > 0 else -3), str(v), ha='center', fontweight='bold')

# PLOT 2: Desolvation Penalty per Residue (From Section 3.2)
# Shows which residues "pay" the highest cost to interact
top_solv = df_int.sort_values('G_Solv', ascending=False).head(10)
sns.barplot(x='G_Solv', y='Label', data=top_solv, ax=axes[1], palette='Reds')
axes[1].set_title('Top 10 Desolvation Penalties (+G_solv)', fontsize=12, fontweight='bold')
axes[1].set_xlabel('Desolvation Energy (kcal/mol)')

# PLOT 3: Critical Hotspots (From Section 3.3)
# Shows residues that actually stabilize the complex
sns.barplot(x='DDG', y='Label', data=df_hot, ax=axes[2], palette='viridis')
axes[2].set_title('CRITICAL HOTSPOTS (Stability Contribution)', fontsize=12, fontweight='bold')
axes[2].set_xlabel('ΔΔG Binding (kcal/mol)')

plt.tight_layout()
plt.savefig('Energy_Analysis_Plots.png', dpi=300)
print("Success! Plots saved as 'Energy_Analysis_Plots.png'.")
