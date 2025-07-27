import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from scipy.stats import sem
import warnings
import matplotlib as mpl

warnings.filterwarnings("ignore")

# === Config ===
folder_path = "large"
# folder_path = "small_100_80"
# folder_path = "subtract_wcrt2"
method_display_names = {
    'Hoogeveen': 'GCP',
    # 'HoogeveenParallel': 'GCP_Par',
    'Greedy': 'Greedy',
    # 'AntColony': 'ACO',
    'Gurobi': 'ILP',
    'Genetic': 'Genetic',
    # 'Hoogeveen+Genetic': 'GCP+Gen',
}
selected_methods = list(method_display_names.keys())

# Define custom markers for each method
method_markers = {
    'GCP': 'o',
    'GCP_Par': 's',
    'Greedy': 'D',
    'ACO': '^',
    'ILP': 'v',
    'Genetic': 'X',
    'GCP+Gen': 'P'
}

# Set global font sizes
mpl.rcParams.update({
    'font.size': 13,              # base font size
    'axes.titlesize': 13,         # title font size
    'axes.labelsize': 13,         # x/y label font size
    'xtick.labelsize': 12,         # x tick label size
    'ytick.labelsize': 12,         # y tick label size
    'legend.fontsize': 12,         # legend font size
    'figure.titlesize': 13        # suptitle font size
})


def get_csv_files(folder_path):
    return [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.csv')]

def read_and_aggregate(folder_path):
    all_dfs = []
    for filepath in get_csv_files(folder_path):
        try:
            df = pd.read_csv(filepath)
            df["DatasetName"] = os.path.basename(df["Dataset"].iloc[0])
            all_dfs.append(df)
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
    return pd.concat(all_dfs, ignore_index=True)

def read_and_aggregate2(folder_path1, folder_path2):
    files1 = get_csv_files(folder_path1)
    files2 = get_csv_files(folder_path2)
    all_dfs = []
    for filepath1, filepath2 in zip(files1, files2):
        try:
            df1 = pd.read_csv(filepath1)
            df1 = df1[df1["Method"] == "Gurobi"]
            df2 = pd.read_csv(filepath2)
            df2 = df2[df2["Method"] != "Gurobi"]
            df = pd.concat([df1, df2], ignore_index=True)
            df["DatasetName"] = os.path.basename(df["Dataset"].iloc[0])
            all_dfs.append(df)
        except Exception as e:
            print(f"Error reading {filepath1} or {filepath2}: {e}")
    return pd.concat(all_dfs, ignore_index=True)


def normalize_group(group):
    min_val = group["SumProb"].min()
    max_val = group["SumProb"].max()
    if max_val == min_val:
        group["NormSumProb"] = 1.0
    else:
        group["NormSumProb"] = (group["SumProb"] - min_val) / (max_val - min_val)
    return group

# === Data Loading and Filtering ===
df_all = read_and_aggregate(folder_path)
budgets = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
df_all = df_all[df_all['Budget'].isin(budgets)]

# df_all = read_and_aggregate2("small_gurobi", "small_100_80")
df_filtered = df_all[df_all['Method'].isin(selected_methods)].copy()
df_filtered["DisplayMethod"] = df_filtered["Method"].map(method_display_names)
df_filtered = df_filtered[(df_filtered['Budget'] > 0) & (df_filtered['Budget'] <= 1200)].copy()

df_normalized = df_filtered.groupby(["DatasetName", "Budget"]).apply(normalize_group).reset_index(drop=True)
df_normalized = df_normalized[(df_normalized['Budget'] > 0) & (df_normalized['Budget'] <= 1200)].copy()

# === Compute Deviation from GCP Baseline ===
baseline_df = df_filtered[df_filtered['DisplayMethod'] == 'GCP'][['DatasetName', 'Budget', 'SumProb']]
baseline_df = baseline_df.rename(columns={'SumProb': 'GCP_SumProb'})

df_with_baseline = df_filtered.merge(baseline_df, on=['DatasetName', 'Budget'])
df_with_baseline['DeviationFromGCP'] = df_with_baseline['SumProb'] - df_with_baseline['GCP_SumProb']

# === Aggregation for Plotting ===
deviation_agg = df_with_baseline.groupby(['DisplayMethod', 'Budget'])['DeviationFromGCP'].agg(['mean', sem]).reset_index()

# 1) Average SumProb (FOM) vs Deadline
agg = df_filtered.groupby(['DisplayMethod', 'Budget'])['SumProb'].agg(['mean', sem]).reset_index()
# plt.figure(figsize=(6, 4))
# for method in agg['DisplayMethod'].unique():
#     sub = agg[agg['DisplayMethod'] == method]
#     plt.plot(sub['Budget'], sub['mean'], label=method)
#     plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
# plt.xlabel("Deadline")
# # plt.ylabel("FOM")
# plt.title("Average FOM vs Deadline")
# plt.legend(title="Method")
# plt.grid(True)
# plt.tight_layout()
# plt.show()

plt.figure(figsize=(6, 3))
for method in agg['DisplayMethod'].unique():
    sub = agg[agg['DisplayMethod'] == method]
    marker = method_markers.get(method, 'o')
    plt.plot(sub['Budget'], sub['mean'], label=method, marker=marker)
    plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("FOM")
# plt.yscale("log")
# plt.title("Average FOM vs Deadline")
plt.legend(title="Method")
plt.grid(True)
plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
plt.savefig("plots/Average_FOM_vs_Deadline.png", dpi=300)
plt.show()

# # 2) Boxplots for SumProb (FOM) at [100, 500, 1000]
# # budgets = [100, 300, 600]
# budgets = [10, 30, 50]
# fig, axes = plt.subplots(1, 3, figsize=(10, 3), sharey=True)
# for i, budget in enumerate(budgets):
#     ax = axes[i]
#     sns.boxplot(data=df_filtered[df_filtered['Budget'] == budget], x="DisplayMethod", y="SumProb", ax=ax)
#     ax.set_title(f"Deadline = {budget}")
#     ax.set_xlabel("Method")
#     ax.tick_params(axis='x', labelrotation=45)
#     if i == 0:
#         ax.set_ylabel("FOM")
#     else:
#         ax.set_ylabel("")
# # plt.suptitle("FOM per Method for Different Deadlines")
# plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig("plots/Box_Average_FOM_vs_Deadline.png", dpi=300)
# plt.show()

# # 3) Normalized SumProb and CI
# df_normalized = df_filtered.groupby(["DatasetName", "Budget"]).apply(normalize_group).reset_index(drop=True)
# agg_norm = df_normalized.groupby(['DisplayMethod', 'Budget'])['NormSumProb'].agg(['mean', sem]).reset_index()
# plt.figure(figsize=(6, 3))
# # for method in agg_norm['DisplayMethod'].unique():
# #     sub = agg_norm[agg_norm['DisplayMethod'] == method]
# #     plt.plot(sub['Budget'], sub['mean'], label=method)
# #     plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
# # plt.xlabel("Deadline")
# # plt.ylabel("Normalized FOM")
# # plt.title("Average Normalized FOM vs Deadline (± CI)")
# # plt.legend(title="Method")
# # plt.grid(True)
# # plt.tight_layout()
# # plt.show()

# for method in agg_norm['DisplayMethod'].unique():
#     sub = agg_norm[agg_norm['DisplayMethod'] == method]
#     marker = method_markers.get(method, 'o')
#     plt.plot(sub['Budget'], sub['mean'], label=method, marker=marker)
#     plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
# plt.xlabel("Deadline")
# plt.ylabel("Normalized FOM")
# # plt.yscale("log")
# # plt.title("Average Normalized FOM vs Deadline (± CI)")
# plt.legend(title="Method")
# plt.grid(True)
# plt.tight_layout()
# plt.savefig("plots/Normalized_Average_FOM_vs_Deadline.png", dpi=300)
# plt.show()

# # 4) Boxplots for Normalized SumProb (FOM) at [100, 500, 1000]
# fig, axes = plt.subplots(1, 3, figsize=(10, 4), sharey=True)
# for i, budget in enumerate(budgets):
#     ax = axes[i]
#     sns.boxplot(data=df_normalized[df_normalized['Budget'] == budget], x="DisplayMethod", y="NormSumProb", ax=ax)
#     ax.set_title(f"Deadline = {budget}")
#     ax.set_xlabel("Method")
#     ax.tick_params(axis='x', labelrotation=45)
#     if i == 0:
#         ax.set_ylabel("Normalized FOM")
#     else:
#         ax.set_ylabel("")
# plt.suptitle("Normalized FOM per Method for Different Deadlines")
# plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig("plots/Box_Normalized_Average_FOM_vs_Deadline.png", dpi=300)
# plt.show()

# 5) Computation Time
agg_time = df_filtered.groupby(['DisplayMethod', 'Budget'])['TimeSec'].agg(['mean', 'min', 'max']).reset_index()
plt.figure(figsize=(6, 3))
for method in agg_time['DisplayMethod'].unique():
    sub = agg_time[agg_time['DisplayMethod'] == method]
    plt.plot(sub['Budget'], sub['mean'], label=method)
    plt.fill_between(sub['Budget'], sub['min'], sub['max'], alpha=0.2)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("Computation Time")
# plt.title("Computation Time vs Deadline (Mean ± Range)")
plt.yscale("log")
plt.legend(title="Method")
plt.grid(True)
plt.xticks(sorted(agg_time['Budget'].unique()))
plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
plt.savefig("plots/Runtime_vs_Deadline.png", dpi=300)
plt.show()

max_runtime_table = agg_time.pivot(index='Budget', columns='DisplayMethod', values='mean')

# Save to CSV or display
# max_runtime_table.to_csv("worst_case_runtime/mean_runtime_table.csv")
print(max_runtime_table)


# 6) Absolute Deviation from GCP (baseline)
baseline = df_filtered[df_filtered['DisplayMethod'] == 'GCP'][['DatasetName', 'Budget', 'SumProb']]
baseline = baseline.rename(columns={'SumProb': 'BaselineFOM'})

df_dev = pd.merge(df_filtered, baseline, on=['DatasetName', 'Budget'])
df_dev['FOMDeviation'] = df_dev['SumProb'] - df_dev['BaselineFOM']

# Aggregate deviation
agg_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMDeviation'].agg(['mean', sem]).reset_index()
plt.figure(figsize=(6, 3))
for method in agg_dev['DisplayMethod'].unique():
    # if method == 'GCP':
    #     continue
    sub = agg_dev[agg_dev['DisplayMethod'] == method]
    plt.plot(sub['Budget'], sub['mean'], label=method)
    plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
plt.axhline(0, linestyle='--', color='black', linewidth=1)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("FOM Deviation from GCP")
# plt.title("Absolute Deviation from GCP")
plt.legend(title="Method")
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/Absolute_Deviation.png", dpi=300)
plt.show()

# 7) Percentage Deviation from GCP (baseline)
df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']

agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()
plt.figure(figsize=(6, 3))
for method in agg_pct_dev['DisplayMethod'].unique():
    # if method == 'GCP':
    #     continue
    sub = agg_pct_dev[agg_pct_dev['DisplayMethod'] == method]
    plt.plot(sub['Budget'], sub['mean'], label=method)
    plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
plt.axhline(0, linestyle='--', color='black', linewidth=1)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("FOM % Deviation from GCP")
# plt.title("Percentage Deviation from GCP")
plt.legend(title="Method")
plt.grid(True)
plt.xticks(sorted(df_dev['Budget'].unique()))
plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
plt.savefig("plots/Percentage_Deviation.png", dpi=300)
plt.show()

import matplotlib.pyplot as plt
from scipy.stats import sem

# --- Compute deviation ---
df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']
agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()
agg_gcp = df_dev[df_dev['DisplayMethod'] == 'GCP'].groupby('Budget')['BaselineFOM'].mean().reset_index()

# --- Plot main figure ---
fig, ax1 = plt.subplots(figsize=(6, 3))

# Plot % deviation with error bars
for method in agg_pct_dev['DisplayMethod'].unique():
    # if method == 'GCP':
    #     continue
    sub = agg_pct_dev[agg_pct_dev['DisplayMethod'] == method]
    ax1.plot(sub['Budget'], sub['mean'], label=method)
    ax1.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)

# Customize left axis
ax1.axhline(0, linestyle='--', color='black', linewidth=1)
ax1.set_xlabel("Deadline (Seconds)")
ax1.set_ylabel("FOM % Deviation from GCP")
ax1.set_xticks(sorted(df_dev['Budget'].unique())[::2])
ax1.grid(True)
ax1.legend(title="Method")
ax1.set_xlim(left=sorted(df_dev['Budget'].unique())[0] - 20)

# --- Add top x-axis showing GCP average FoM ---
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())  # Align limits

# Set tick positions and labels on top
budgets = agg_gcp['Budget']
avg_foms = agg_gcp['BaselineFOM'].round(3)
ax2.set_xticks(budgets[::2])
ax2.set_xticklabels(avg_foms[::2])
ax2.set_xlabel("Mean FoM of GCP", labelpad=6)

# Final layout and save
plt.tight_layout()
plt.savefig("plots/Percentage_Deviation_with_GCP_FoM.png", dpi=300)
plt.show()



import matplotlib.pyplot as plt
from scipy.stats import sem

# Compute deviation
df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']
agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()

agg_time = df_filtered.groupby(['DisplayMethod', 'Budget'])['TimeSec'].agg(['mean', 'min', 'max']).reset_index()

# Prepare x ticks
x_ticks = sorted(df_filtered['Budget'].unique())

# Create side-by-side plot
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.5), dpi=300)
fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(7.2, 3.5), dpi=300)

# --- Left Plot: Computation Time ---
for method in agg_time['DisplayMethod'].unique():
    # if method == "GCP":
    #     continue
    sub = agg_time[agg_time['DisplayMethod'] == method]
    ax1.plot(sub['Budget'], sub['mean'], label=method)
    ax1.fill_between(sub['Budget'], sub['min'], sub['max'], alpha=0.2)
ax1.set_xlabel("Deadline (Seconds)")
# ax1.set_ylabel("Computation Time (log scale)")
ax1.set_yscale("log")
ax1.set_xticks(x_ticks[::3])
ax1.set_title("Computation Time")
ax1.grid(True)

# --- Right Plot: % Deviation ---
for method in agg_pct_dev['DisplayMethod'].unique():
    if method == "Greedy":
        continue
    sub = agg_pct_dev[agg_pct_dev['DisplayMethod'] == method]
    ax2.plot(sub['Budget'], sub['mean'], label=method)
    ax2.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
ax2.axhline(0, linestyle='--', color='black', linewidth=1)
ax2.set_xlabel("Deadline (Seconds)")
# ax2.set_ylabel("FOM % Deviation from GCP")
ax2.set_xticks(x_ticks[::3])
ax2.set_title("Expectation % Deviation")
ax2.grid(True)

# --- Shared Legend ---
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.75, 0.35), ncol=1, frameon=False)

# Layout and save
plt.tight_layout(rect=[0, 0, 1, 0.93])  # leave space on top for legend
plt.savefig("plots/Combined_Runtime_Deviation.pdf", dpi=300, bbox_inches='tight')
plt.show()

# import matplotlib.pyplot as plt
# from scipy.stats import sem
# import pandas as pd

# # Compute deviation
# df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']

# # Aggregate
# agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()
# agg_gcp = df_dev[df_dev['DisplayMethod'] == 'GCP'].groupby('Budget')['BaselineFOM'].mean().reset_index()

# # Print check
# print(agg_pct_dev.head())
# print("Methods found:", agg_pct_dev['DisplayMethod'].unique())

# # Setup plot
# fig, ax1 = plt.subplots(figsize=(6, 3))

# # Plot deviation lines
# methods = [m for m in agg_pct_dev['DisplayMethod'].unique() if m != 'GCP']
# if not methods:
#     print("No non-GCP methods to plot.")
# else:
#     for method in methods:
#         sub = agg_pct_dev[agg_pct_dev['DisplayMethod'] == method]
#         if not sub.empty:
#             ax1.plot(sub['Budget'], sub['mean'], label=method)
#             ax1.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)

#     ax1.axhline(0, linestyle='--', color='black', linewidth=1)
#     ax1.set_xlabel("Deadline")
#     ax1.set_ylabel("FOM % Deviation from GCP")
#     budgets = sorted(df_dev['Budget'].unique())
#     ax1.set_xticks(budgets[::2])
#     ax1.grid(True)
#     ax1.legend(title="Method")

#     # Add top x-axis for average GCP FoM
#     ax2 = ax1.twiny()
#     ax2.set_xlim(ax1.get_xlim())

#     avg_foms = agg_gcp.set_index('Budget').loc[budgets]['BaselineFOM'].round(4)
#     ax2.set_xticks(budgets[::2])
#     ax2.set_xticklabels(avg_foms[::2])
#     ax2.set_xlabel("Average GCP FoM")

#     fig.subplots_adjust(top=0.82)
#     plt.tight_layout()
#     plt.savefig("plots/Percentage_Deviation_with_GCP_FoM.png", dpi=300)
#     plt.show()


# Switch baseline from GCP to ILP
baseline = df_filtered[df_filtered['DisplayMethod'] == 'ILP'][['DatasetName', 'Budget', 'SumProb']]
baseline = baseline.rename(columns={'SumProb': 'BaselineFOM'})

df_dev = pd.merge(df_filtered, baseline, on=['DatasetName', 'Budget'])
df_dev['FOMDeviation'] = df_dev['SumProb'] - df_dev['BaselineFOM']

# Absolute Deviation Plot
agg_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMDeviation'].agg(['mean', sem]).reset_index()
plt.figure(figsize=(6, 3))
for method in agg_dev['DisplayMethod'].unique():
    if method == 'ILP':
        continue
    sub = agg_dev[agg_dev['DisplayMethod'] == method]
    plt.plot(sub['Budget'], sub['mean'], label=method)
    plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
plt.axhline(0, linestyle='--', color='black', linewidth=1)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("FOM Deviation from ILP")
plt.legend(title="Method")
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/Absolute_Deviation_from_ILP.png", dpi=300)
plt.show()

# Percentage Deviation Plot
df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']
agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()

plt.figure(figsize=(6, 3))
for method in agg_pct_dev['DisplayMethod'].unique():
    if method == 'ILP':
        continue
    sub = agg_pct_dev[agg_pct_dev['DisplayMethod'] == method]
    plt.plot(sub['Budget'], sub['mean'], label=method)
    plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
plt.axhline(0, linestyle='--', color='black', linewidth=1)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("Expectation % Deviation from ILP")
plt.legend(title="Method")
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/Percentage_Deviation_from_ILP.png", dpi=300)
plt.show()

