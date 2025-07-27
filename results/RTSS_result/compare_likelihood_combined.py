# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Define file lists
# file_lists = {
#     "Single": [
#         ("out_single_3.csv", "S1"), ("out_single_2.csv", "S2"), ("out_single_1.csv", "S3"),
#         ("out_single_6.csv", "S4"), ("out_single_5.csv", "S5"), ("out_single_4.csv", "S6"),
#         ("out_single_9.csv", "S7"), ("out_single_8.csv", "S8"), ("out_single_7.csv", "S9"),
#     ],
#     "Single_Deepslow": [
#         ("out_single_deepslow_3.csv", "S1"), ("out_single_deepslow_2.csv", "S2"),
#         ("out_single_deepslow_1.csv", "S3"), ("out_single_deepslow_6.csv", "S4"),
#         ("out_single_deepslow_5.csv", "S5"), ("out_single_deepslow_4.csv", "S6"),
#         ("out_single_deepslow_9.csv", "S7"), ("out_single_deepslow_8.csv", "S8"),
#         ("out_single_deepslow_7.csv", "S9"),
#     ],
#     "Combined": [
#         ("out_combined_3.csv", "S1"), ("out_combined_2.csv", "S2"), ("out_combined_1.csv", "S3"),
#         ("out_combined_6.csv", "S4"), ("out_combined_5.csv", "S5"), ("out_combined_4.csv", "S6"),
#         ("out_combined_7.csv", "S9"),
#     ],
#     "Combined_Deepslow": [
#         ("out_combined_deepslow_3.csv", "S1"), ("out_combined_deepslow_2.csv", "S2"),
#         ("out_combined_deepslow_1.csv", "S3"), ("out_combined_deepslow_6.csv", "S4"),
#         ("out_combined_deepslow_5.csv", "S5"), ("out_combined_deepslow_4.csv", "S6"),
#         ("out_combined_deepslow_7.csv", "S9"),
#     ]
# }

# methods_to_compare = ["Gurobi", "Greedy", "Genetic", "AntColony", "SimulatedAnnealing"]

# # Function to process data
# def process_data(files):
#     all_data = []
#     for file, setting_name in files:
#         try:
#             df = pd.read_csv(file)
#             df["setting"] = setting_name
#             all_data.append(df)
#         except FileNotFoundError:
#             print(f"Warning: File {file} not found. Skipping.")
#     df_all = pd.concat(all_data)
#     df_all.columns = df_all.columns.str.lower()
#     summary = df_all.groupby(["setting", "method"]).agg(
#         mean_prob=("sumprob", "mean"),
#         std_prob=("sumprob", "std"),
#         count=("sumprob", "count")
#     ).reset_index()
#     summary = summary[summary["method"].isin(methods_to_compare)]
#     setting_order = [s for _, s in files]
#     summary["setting"] = pd.Categorical(summary["setting"], categories=setting_order, ordered=True)
#     return summary, setting_order

# # Set up figure for two-column paper (width ~6.5 inches, height adjusted for 2x2 grid)
# fig = plt.figure(figsize=(6.5, 5.5))  # Compact size for two-column format

# # Create 2x2 subplot grid
# axes = []
# for i in range(4):
#     ax = fig.add_subplot(2, 2, i + 1, projection='3d')
#     axes.append(ax)

# # Global min and max for consistent color scaling
# z_min, z_max = float('inf'), float('-inf')

# # Process and plot data for each subplot
# surfaces = []
# for idx, (title, files) in enumerate(file_lists.items()):
#     summary, setting_order = process_data(files)
#     pivot_table = summary.pivot(index="setting", columns="method", values="mean_prob")
#     Z = pivot_table.values
#     X, Y = np.meshgrid(range(len(methods_to_compare)), range(len(setting_order)))
    
#     # Update global z limits
#     z_min = min(z_min, Z.min())
#     z_max = max(z_max, Z.max())
    
#     # Plot surface
#     surf = axes[idx].plot_surface(X, Y, Z, cmap='viridis', edgecolor='gray', linewidth=0.2, alpha=0.9)
#     surfaces.append(surf)
    
#     # Customize axes
#     # axes[idx].set_title(title, fontsize=10, pad=5)
#     axes[idx].set_xticks(range(len(methods_to_compare)))
#     axes[idx].set_xticklabels(methods_to_compare, rotation=45, ha='right', fontsize=8)
#     axes[idx].set_yticks(range(len(setting_order)))
#     axes[idx].set_yticklabels(setting_order, rotation=-30, ha='left', fontsize=8)
#     axes[idx].view_init(elev=15, azim=-45)
#     axes[idx].tick_params(axis='z', labelsize=7)
#     # axes[idx].set_zlabel('Collected Likelihood', fontsize=8, labelpad=5)

# # Add a single shared color bar
# cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
# cbar = fig.colorbar(surfaces[0], cax=cbar_ax)
# cbar.set_label('Collected Likelihood', fontsize=10)
# cbar.ax.tick_params(labelsize=8)

# # Set consistent z-limits across all subplots
# for ax in axes:
#     ax.set_zlim(z_min, z_max)

# # Adjust layout to prevent overlap and fit color bar
# plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.05, wspace=0.2, hspace=0.3)
# plt.savefig("four_3d_surface_plots_shared_cbar.png", dpi=300, bbox_inches='tight')
# plt.show()



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define file lists (only two datasets)
file_lists = {
    "Single": [
        ("out_single_3.csv", "S1"), ("out_single_2.csv", "S2"), ("out_single_1.csv", "S3"),
        ("out_single_6.csv", "S4"), ("out_single_5.csv", "S5"), ("out_single_4.csv", "S6"),
        ("out_single_9.csv", "S7"), ("out_single_8.csv", "S8"), ("out_single_7.csv", "S9"),
    ],
    "Single-Deepslow": [
        ("out_single_deepslow_3.csv", "S1"), ("out_single_deepslow_2.csv", "S2"),
        ("out_single_deepslow_1.csv", "S3"), ("out_single_deepslow_6.csv", "S4"),
        ("out_single_deepslow_5.csv", "S5"), ("out_single_deepslow_4.csv", "S6"),
        ("out_single_deepslow_9.csv", "S7"), ("out_single_deepslow_8.csv", "S8"),
        ("out_single_deepslow_7.csv", "S9"),
    ]
}

# Shortened method names
methods_to_compare = ["Gurobi", "Greedy", "Genetic", "AntColony", "SimulatedAnnealing"]
short_methods = ["Gurobi", "Greedy", "Genetic", "AntColony", "Annealing"]  # Shortened versions

# Function to process data
def process_data(files):
    all_data = []
    for file, setting_name in files:
        try:
            df = pd.read_csv(file)
            df["setting"] = setting_name
            all_data.append(df)
        except FileNotFoundError:
            print(f"Warning: File {file} not found. Skipping.")
    df_all = pd.concat(all_data)
    df_all.columns = df_all.columns.str.lower()
    summary = df_all.groupby(["setting", "method"]).agg(
        mean_prob=("sumprob", "mean"),
        std_prob=("sumprob", "std"),
        count=("sumprob", "count")
    ).reset_index()
    summary = summary[summary["method"].isin(methods_to_compare)]
    setting_order = [s for _, s in files]
    summary["setting"] = pd.Categorical(summary["setting"], categories=setting_order, ordered=True)
    return summary, setting_order

# Set up figure for two-column paper (width ~6.5 inches, height for 1x2 layout)
fig = plt.figure(figsize=(6.5, 3.0))  # Adjusted height for two side-by-side plots

# Create 1x2 subplot grid (side by side)
axes = []
for i in range(2):
    ax = fig.add_subplot(1, 2, i + 1, projection='3d')
    axes.append(ax)

# Process and plot data for each subplot
for idx, (title, files) in enumerate(file_lists.items()):
    summary, setting_order = process_data(files)
    pivot_table = summary.pivot(index="setting", columns="method", values="mean_prob")
    Z = pivot_table.values
    X, Y = np.meshgrid(range(len(methods_to_compare)), range(len(setting_order)))
    
    # Plot surface without color bar
    axes[idx].plot_surface(X, Y, Z, cmap='viridis', edgecolor='gray', linewidth=0.2, alpha=0.9)
    
    # Customize axes with shortened method names
    # axes[idx].set_title(title, fontsize=11)
    axes[idx].set_xticks(range(len(short_methods)))
    axes[idx].set_xticklabels(short_methods, rotation=90, ha='right', fontsize=10)
    axes[idx].set_yticks(range(len(setting_order)))
    axes[idx].set_yticklabels(setting_order, rotation=-90, ha='left', fontsize=8)
    axes[idx].set_zlabel('Collected Likelihood', fontsize=10)
    axes[idx].view_init(elev=15, azim=-42)
    axes[idx].tick_params(axis='z', labelsize=9)

# Adjust layout without color bar
plt.tight_layout()
# plt.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.15, wspace=0.2)
plt.savefig("two_3d_surface_plots_short_names.png", dpi=300, bbox_inches='tight')
plt.show()



# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Define file lists (only two datasets)
# file_lists = {
#     "Single": [
#         ("out_single_3.csv", "S1"), ("out_single_2.csv", "S2"), ("out_single_1.csv", "S3"),
#         ("out_single_6.csv", "S4"), ("out_single_5.csv", "S5"), ("out_single_4.csv", "S6"),
#         ("out_single_9.csv", "S7"), ("out_single_8.csv", "S8"), ("out_single_7.csv", "S9"),
#     ],
#     "Single-Deepslow": [
#         ("out_single_deepslow_3.csv", "S1"), ("out_single_deepslow_2.csv", "S2"),
#         ("out_single_deepslow_1.csv", "S3"), ("out_single_deepslow_6.csv", "S4"),
#         ("out_single_deepslow_5.csv", "S5"), ("out_single_deepslow_4.csv", "S6"),
#         ("out_single_deepslow_9.csv", "S7"), ("out_single_deepslow_8.csv", "S8"),
#         ("out_single_deepslow_7.csv", "S9"),
#     ]
# }

# # Method names (slightly adjusted for clarity)
# methods_to_compare = ["Gurobi", "Greedy", "Genetic", "AntColony", "SimulatedAnnealing"]
# short_methods = ["Gurobi", "Greedy", "Genetic", "AntColony", "Annealing"]

# # Function to process data
# def process_data(files):
#     all_data = []
#     for file, setting_name in files:
#         try:
#             df = pd.read_csv(file)
#             df["setting"] = setting_name
#             all_data.append(df)
#         except FileNotFoundError:
#             print(f"Warning: File {file} not found. Skipping.")
#     df_all = pd.concat(all_data)
#     df_all.columns = df_all.columns.str.lower()
#     summary = df_all.groupby(["setting", "method"]).agg(
#         mean_prob=("sumprob", "mean"),
#         std_prob=("sumprob", "std"),
#         count=("sumprob", "count")
#     ).reset_index()
#     summary = summary[summary["method"].isin(methods_to_compare)]
#     setting_order = [s for _, s in files]
#     summary["setting"] = pd.Categorical(summary["setting"], categories=setting_order, ordered=True)
#     return summary, setting_order

# # Set up figure for two-column paper (width ~6.5 inches, height for 1x2 layout)
# fig = plt.figure(figsize=(6.5, 3.0))  # Adjusted height for two side-by-side plots

# # Create 1x2 subplot grid (side by side)
# axes = []
# for i in range(2):
#     ax = fig.add_subplot(1, 2, i + 1, projection='3d')
#     axes.append(ax)

# # Process and plot data for each subplot
# for idx, (title, files) in enumerate(file_lists.items()):
#     summary, setting_order = process_data(files)
#     pivot_table = summary.pivot(index="setting", columns="method", values="mean_prob")
#     Z = pivot_table.values
#     X, Y = np.meshgrid(range(len(methods_to_compare)), range(len(setting_order)))
    
#     # Plot surface without color bar
#     axes[idx].plot_surface(X, Y, Z, cmap='viridis', edgecolor='gray', linewidth=0.2, alpha=0.9)
    
#     # Customize axes with method names
#     axes[idx].set_xticks(range(len(short_methods)))
#     axes[idx].set_xticklabels(short_methods, rotation=90, ha='center', fontsize=10)
#     axes[idx].set_yticks(range(len(setting_order)))
#     axes[idx].set_yticklabels(setting_order, rotation=-90, ha='right', fontsize=8)
#     axes[idx].view_init(elev=15, azim=-45)
#     axes[idx].tick_params(axis='z', labelsize=7)

# plt.tight_layout()
# plt.subplots_adjust(top=0.98)
# plt.savefig("two_3d_surface_plots_short_names_tight.png", dpi=300, bbox_inches='tight')
# plt.show()