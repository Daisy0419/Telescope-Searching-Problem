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
folder_path = "subtract_wcrt2"
method_display_names = {
    'Hoogeveen': 'GCP',
    # 'HoogeveenParallel': 'GCP_Par',
    'Greedy': 'Greedy',
    # 'AntColony': 'ACO',
    # 'Gurobi': 'ILP',
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
# budgets = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
# df_all = df_all[df_all['Budget'].isin(budgets)]
# df_all = read_and_aggregate2("small_gurobi", "small_100_80")
df_filtered = df_all[df_all['Method'].isin(selected_methods)].copy()
df_filtered["DisplayMethod"] = df_filtered["Method"].map(method_display_names)
# df_filtered = df_filtered[(df_filtered['Budget'] > 0) & (df_filtered['Budget'] <= 1200)].copy()

df_normalized = df_filtered.groupby(["DatasetName", "Budget"]).apply(normalize_group).reset_index(drop=True)
# df_normalized = df_normalized[(df_normalized['Budget'] > 0) & (df_normalized['Budget'] <= 1200)].copy()

#Compute Deviation from GCP Baseline
baseline_df = df_filtered[df_filtered['DisplayMethod'] == 'GCP'][['DatasetName', 'Budget', 'SumProb']]
baseline_df = baseline_df.rename(columns={'SumProb': 'GCP_SumProb'})

df_with_baseline = df_filtered.merge(baseline_df, on=['DatasetName', 'Budget'])
df_with_baseline['DeviationFromGCP'] = df_with_baseline['SumProb'] - df_with_baseline['GCP_SumProb']

# Aggregation for Plotting
deviation_agg = df_with_baseline.groupby(['DisplayMethod', 'Budget'])['DeviationFromGCP'].agg(['mean', sem]).reset_index()


#This plot is used in RTSS Paper (Figure 9)
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
# plt.title("Average FOM vs Deadline")
plt.legend(title="Method")
plt.grid(True)
# plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
# plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
plt.savefig("plots/Average_FOM_vs_Deadline.png", dpi=300)
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
# plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
# plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
# plt.savefig("plots/Runtime_vs_Deadline.png", dpi=300)
plt.show()


#print max runtime
max_runtime_table = agg_time.pivot(index='Budget', columns='DisplayMethod', values='mean')
# Save to CSV or display
# max_runtime_table.to_csv("worst_case_runtime/mean_runtime_table.csv")
print('===== Maximum Runtime =====')
print(max_runtime_table)


# 3) Absolute Deviation from GCP (baseline)
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
plt.title("Absolute Deviation from GCP")
plt.legend(title="Method")
plt.grid(True)
plt.tight_layout()
# plt.savefig("plots/Absolute_Deviation.png", dpi=300)
plt.show()

# 4) Percentage Deviation from GCP (baseline)
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
plt.title("Percentage Deviation from GCP")
plt.legend(title="Method")
plt.grid(True)
plt.xticks(sorted(df_dev['Budget'].unique()))
# plt.xlim(left = sorted(agg['Budget'].unique())[0]-20)
x_ticks = sorted(agg['Budget'].unique())
# plt.xticks(ticks=x_ticks[::2])
plt.tight_layout()
# plt.savefig("plots/Percentage_Deviation.png", dpi=300)
plt.show()


# This plot is used in the RTSS paper (Figure 10)
# 4) Percentage Deviation from ILP with mean FOM
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

#left axis
ax1.axhline(0, linestyle='--', color='black', linewidth=1)
ax1.set_xlabel("Deadline (Seconds)")
ax1.set_ylabel("FOM % Deviation from GCP")
ax1.set_xticks(sorted(df_dev['Budget'].unique()))
ax1.grid(True)
# ax1.legend(title="Method")
ax1.legend(title="Method", ncol=2, loc='lower left', bbox_to_anchor=(0.3, 0.05))
# ax1.set_xlim(left=sorted(df_dev['Budget'].unique())[0] - 20)

# top x-axis showing GCP average FoM
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())  # Align limits

budgets = agg_gcp['Budget']
avg_foms = agg_gcp['BaselineFOM'].round(3)
ax2.set_xticks(budgets[::2])
ax2.set_xticklabels(avg_foms[::2])
ax2.set_xlabel("Mean FOM of GCP", labelpad=6)

plt.tight_layout()
plt.savefig("plots/Percentage_Deviation_with_GCP_FoM.pdf", dpi=300)
plt.show()
