import numpy as np
import matplotlib.pyplot as plt

# --- File paths ---
dbns_file = "forceCoeffs.dat"
exp_file = "reference_data.dat"

output_pdf = "alpha_vs_Cl_simple.pdf"

# --- Load DBNS data ---
# Assumes forceCoeffs.dat has columns: t, Cm, q, alpha, Ma, Cn
# and starts from row 1 (skipping any headers or comments)
t_list, Cm_list, Cd_list, Cl_list = [], [], [], []

with open(dbns_file, 'r') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        t_list.append(float(parts[0]))
        Cd_list.append(float(parts[1]))  # Cd in your Gnuplot ($3)
        Cl_list.append(float(parts[2]))  # Cd in your Gnuplot ($3)
        Cm_list.append(float(parts[3]))  # Cl in your Gnuplot ($2)

t = np.array(t_list)
Cd = np.array(Cd_list)
Cl = np.array(Cl_list)
Cm = np.array(Cm_list)

# --- Compute alpha and corrected Cl ---
alpha_deg = 0.016 - 2.51 * np.sin(39.978 * t)
alpha_rad = np.deg2rad(alpha_deg)
Cn = Cd * np.sin(alpha_rad) + Cl * np.cos(alpha_rad)

# --- Load experimental data ---
alpha_exp, Cn_exp = [], []

with open(exp_file, 'r') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        parts = line.strip().split()
        if len(parts) < 7:
            continue
        alpha_exp.append(float(parts[4]))  # column 5 = alpha
        Cn_exp.append(float(parts[6]))     # column 7 = Cl

# --- Plot ---
plt.figure(figsize=(7, 6))
plt.plot(alpha_deg, Cn, label='DBNS', color='blue', linewidth=2)
plt.scatter(alpha_exp, Cn_exp, label='Experiment', color='red', edgecolors='black', s=50)

plt.title("Alpha vs Cn", fontsize=14)
plt.xlabel("Alpha (degrees)", fontsize=12)
plt.ylabel("Cl", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim([-3, 3])
plt.ylim([-0.4, 0.6])
plt.legend(loc='upper right')
plt.tight_layout()

# Save as PDF
plt.savefig(output_pdf)
print(f"Plot saved to {output_pdf}")
