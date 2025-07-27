import pandas as pd
import matplotlib.pyplot as plt

# file = "out_GW200311_115853.csv"
# file = "out_GW200216_220804.csv"
# file = "out_GW200129_065458_seperate.csv"
# file = "out_GW200105_162426.csv"
file = "out.csv"


# df = pd.read_csv(file)

# methods_to_compare = ['Hoogeveen','Greedy','AntColony','Gurobi','Genetic']
# # methods_to_compare = ['s-Hoogeveen','s-Greedy', 'st-Hoogeveen','st-Greedy']
# # methods_to_compare = ['mstNaive','mstHoogeveen', 'mstLemon', 'mstLemonHoogeveen', 'mstLemonHoogeveen2']

# fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:pink', 'tab:cyan']
# markers = ['o', 's', 'D', 'p', 'x', 'o']

# #SumProb vs Budget
# for method, color, marker in zip(methods_to_compare, colors, markers):
#     subset = df[df['Method'] == method]
#     subset = subset.sort_values('Budget')
#     axes[0].plot(subset['Budget'], subset['SumProb'], label=method, color=color, marker=marker)
    
# axes[0].set_ylabel('SumProb')
# axes[0].set_title('SumProb vs Deadline')
# # axes[0].set_yscale('log')
# axes[0].legend()
# axes[0].grid(True)

# #TimeSec vs Budget
# for method, color, marker in zip(methods_to_compare, colors, markers):
#     subset = df[df['Method'] == method]
#     subset = subset.sort_values('Budget')
#     axes[1].plot(subset['Budget'], subset['TimeSec'], label=method, color=color, marker=marker)
    
# axes[1].set_xlabel('Deadline')
# # axes[1].set_xlim(0, 500)
# axes[1].set_ylabel('TimeSec')
# axes[1].set_yscale('log')
# axes[1].set_title('TimeSec vs Deadline')
# # axes[1].set_title('large dataset (1768 nodes)')
# axes[1].legend()
# axes[1].grid(True)

# plt.tight_layout()
# plt.show()

# import pandas as pd
# import matplotlib.pyplot as plt

# # Load the data
# file = "out.csv"
df = pd.read_csv(file)

# Original method names and corresponding display labels
method_display_names = {
    'Hoogeveen': 'GCP',
    'Greedy': 'Greedy',
    'AntColony': 'ACO',
    'Gurobi': 'ILP',
    'Genetic': 'Genetic'
}

methods_to_compare = list(method_display_names.keys())
display_labels = [method_display_names[m] for m in methods_to_compare]

# Plot settings
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:pink']
markers = ['o', 's', 'D', 'p', 'x']

# --- Plot 1: FoM vs Deadline ---
plt.figure(figsize=(6, 3), dpi=300)

for method, label, color, marker in zip(methods_to_compare, display_labels, colors, markers):
    subset = df[df['Method'] == method].sort_values('Budget')
    plt.plot(subset['Budget'], subset['SumProb'], label=label,
             color=color, marker=marker, markersize=4, linewidth=1)

plt.xlabel('Deadline (s)', fontsize=12)
plt.ylabel('Figure of Merit (FoM)', fontsize=12)
# plt.yscale('log')
# plt.title('FoM vs Deadline', fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("FoM_vs_Deadline.png", dpi=300)
plt.show()
plt.close()

# --- Plot 2: Runtime vs Deadline ---
plt.figure(figsize=(6, 3), dpi=300)

for method, label, color, marker in zip(methods_to_compare, display_labels, colors, markers):
    subset = df[df['Method'] == method].sort_values('Budget')
    plt.plot(subset['Budget'], subset['TimeSec'], label=label,
             color=color, marker=marker, markersize=4, linewidth=1)

plt.xlabel('Deadline (s)', fontsize=12)
plt.ylabel('Computation Time (s)', fontsize=12)
plt.yscale('log')
# plt.title('Runtime vs Deadline', fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("Runtime_vs_Deadline.png", dpi=300)
plt.show()
plt.close()



# Ensure consistent types for Budget
df['Budget'] = df['Budget'].astype(float)

methods_to_compare = list(method_display_names.keys())
display_labels = [method_display_names[m] for m in methods_to_compare]

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:pink']
markers = ['o', 's', 'D', 'p', 'x']

# Extract ILP (baseline) reference curve
df_ref = df[df['Method'] == 'Gurobi'].set_index('Budget')['SumProb']

# Create delta plot
plt.figure(figsize=(6, 3), dpi=300)

for method, label, color, marker in zip(methods_to_compare, display_labels, colors, markers):
    if method == 'Gurobi':  # skip baseline
        continue
    subset = df[df['Method'] == method].sort_values('Budget').set_index('Budget')
    
    # Align indices safely using reindex()
    aligned_ref = df_ref.reindex(subset.index)
    delta_fom = subset['SumProb'] - aligned_ref

    plt.plot(subset.index, delta_fom, label=label,
             color=color, marker=marker, markersize=4, linewidth=1)

plt.xlabel('Deadline (s)', fontsize=12)
plt.ylabel('ΔFoM (vs ILP)', fontsize=12)
plt.axhline(0, linestyle='--', color='gray', linewidth=0.5)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("Delta_FoM_vs_Deadline.png", dpi=300)
plt.show()