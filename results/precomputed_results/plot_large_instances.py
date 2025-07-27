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
method_display_names = {
    'Hoogeveen': 'GCP',
    'Greedy': 'Greedy',
    'Genetic': 'Genetic',
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

folder_path = "large"
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

# Data Loading and Filtering
df_all = read_and_aggregate(folder_path)
budgets = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
df_all = df_all[df_all['Budget'].isin(budgets)]

# df_all = read_and_aggregate2("small_gurobi", "small_100_80")
df_filtered = df_all[df_all['Method'].isin(selected_methods)].copy()
df_filtered["DisplayMethod"] = df_filtered["Method"].map(method_display_names)
df_filtered = df_filtered[(df_filtered['Budget'] > 0) & (df_filtered['Budget'] <= 1200)].copy()

# Compute Deviation from GCP Baseline
baseline_df = df_filtered[df_filtered['DisplayMethod'] == 'GCP'][['DatasetName', 'Budget', 'SumProb']]
baseline_df = baseline_df.rename(columns={'SumProb': 'GCP_SumProb'})

df_with_baseline = df_filtered.merge(baseline_df, on=['DatasetName', 'Budget'])
df_with_baseline['DeviationFromGCP'] = df_with_baseline['SumProb'] - df_with_baseline['GCP_SumProb']

#Aggregation for Plotting
deviation_agg = df_with_baseline.groupby(['DisplayMethod', 'Budget'])['DeviationFromGCP'].agg(['mean', sem]).reset_index()

# 1) Average SumProb (FOM) vs Deadline
agg = df_filtered.groupby(['DisplayMethod', 'Budget'])['SumProb'].agg(['mean', sem]).reset_index()
plt.figure(figsize=(6, 3))
for method in agg['DisplayMethod'].unique():
    sub = agg[agg['DisplayMethod'] == method]
    marker = method_markers.get(method, 'o')
    plt.plot(sub['Budget'], sub['mean'], label=method, marker=marker)
    plt.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("FOM")
# plt.yscale("log")
plt.title("Average FOM vs Deadline")
plt.legend(title="Method")
plt.grid(True)
plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
# plt.savefig("plots/Average_FOM_vs_Deadline.png", dpi=300)
plt.show()

# 2) Computation Time
agg_time = df_filtered.groupby(['DisplayMethod', 'Budget'])['TimeSec'].agg(['mean', 'min', 'max']).reset_index()
plt.figure(figsize=(6, 3))
for method in agg_time['DisplayMethod'].unique():
    sub = agg_time[agg_time['DisplayMethod'] == method]
    plt.plot(sub['Budget'], sub['mean'], label=method)
    plt.fill_between(sub['Budget'], sub['min'], sub['max'], alpha=0.2)
plt.xlabel("Deadline (Seconds)")
plt.ylabel("Computation Time")
plt.title("Computation Time vs Deadline (Mean ± Range)")
plt.yscale("log")
plt.legend(title="Method")
plt.grid(True)
plt.xticks(sorted(agg_time['Budget'].unique()))
plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
# plt.savefig("plots/Runtime_vs_Deadline.png", dpi=300)
plt.show()

#Print maximum runtime
max_runtime_table = agg_time.pivot(index='Budget', columns='DisplayMethod', values='max')
# Save to CSV or display
# max_runtime_table.to_csv("worst_case_runtime/mean_runtime_table.csv")
print('=====Maximum Runtime=====')
print(max_runtime_table)


# 3) Percentage Deviation from GCP (baseline) with mean FOM on top
baseline = df_filtered[df_filtered['DisplayMethod'] == 'GCP'][['DatasetName', 'Budget', 'SumProb']]
baseline = baseline.rename(columns={'SumProb': 'BaselineFOM'})

df_dev = pd.merge(df_filtered, baseline, on=['DatasetName', 'Budget'])
df_dev['FOMDeviation'] = df_dev['SumProb'] - df_dev['BaselineFOM']

df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']
agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()
agg_gcp = df_dev[df_dev['DisplayMethod'] == 'GCP'].groupby('Budget')['BaselineFOM'].mean().reset_index()

fig, ax1 = plt.subplots(figsize=(6, 3))

# Plot % deviation with error bars
for method in agg_pct_dev['DisplayMethod'].unique():
    # if method == 'GCP':
    #     continue
    sub = agg_pct_dev[agg_pct_dev['DisplayMethod'] == method]
    ax1.plot(sub['Budget'], sub['mean'], label=method)
    ax1.fill_between(sub['Budget'], sub['mean'] - sub['sem'], sub['mean'] + sub['sem'], alpha=0.2)

# left axis
ax1.axhline(0, linestyle='--', color='black', linewidth=1)
ax1.set_xlabel("Deadline (Seconds)")
ax1.set_ylabel("FOM % Deviation from GCP")
ax1.set_xticks(sorted(df_dev['Budget'].unique())[::2])
ax1.grid(True)
ax1.legend(title="Method")
ax1.set_xlim(left=sorted(df_dev['Budget'].unique())[0] - 20)

# top x-axis showing GCP average FoM
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

budgets = agg_gcp['Budget']
avg_foms = agg_gcp['BaselineFOM'].round(3)
ax2.set_xticks(budgets[::2])
ax2.set_xticklabels(avg_foms[::2])
ax2.set_xlabel("Mean FoM of GCP", labelpad=6)

plt.tight_layout()
# plt.savefig("plots/Percentage_Deviation_with_GCP_FoM.png", dpi=300)
plt.show()

# This plot is used on RTSS Paper(Figure 8)
# 4) combining Percentage Deviation with Runtime
# Compute deviation
df_dev['FOMPctDeviation'] = 100 * (df_dev['SumProb'] - df_dev['BaselineFOM']) / df_dev['BaselineFOM']
agg_pct_dev = df_dev.groupby(['DisplayMethod', 'Budget'])['FOMPctDeviation'].agg(['mean', sem]).reset_index()
agg_time = df_filtered.groupby(['DisplayMethod', 'Budget'])['TimeSec'].agg(['mean', 'min', 'max']).reset_index()

# Prepare x ticks
x_ticks = sorted(df_filtered['Budget'].unique())
fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(7.2, 3.5), dpi=300)

# subplot: Computation Time
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

# subplot: % Deviation
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
ax2.set_title("FOM % Deviation")
ax2.grid(True)

# Shared Legend
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.75, 0.35), ncol=1, frameon=False)

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("plots/Combined_Runtime_Deviation.pdf", dpi=300, bbox_inches='tight')
plt.show()
