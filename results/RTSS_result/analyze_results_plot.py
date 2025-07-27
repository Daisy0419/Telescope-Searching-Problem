import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

csv_file = "out.csv"
df = pd.read_csv(csv_file, quotechar='"', skipinitialspace=True)

df['SumProb'] = pd.to_numeric(df['SumProb'])
df['TimeSec'] = pd.to_numeric(df['TimeSec'])

method_colors = {
    "Genetic": "tab:green",
    "Greedy": "tab:orange",
    "Gurobi": "tab:purple",
    "MST": "tab:blue",
    "Christo_lemon": "tab:brown",
    "SimulatedAnnealing": "tab:pink",
    "op_solver": "tab:cyan"
}

grouped = df.groupby(["Dataset", "Budget", "SlewRate", "DwellTime", "Method"]).agg({
    "SumProb": "mean",
    "TimeSec": "mean"
}).reset_index()

unique_settings = grouped[["Dataset", "Budget", "SlewRate", "DwellTime"]].drop_duplicates()

for idx, (dataset, budget, slew_rate, dwell_time) in enumerate(unique_settings.itertuples(index=False)):
    setting_df = grouped[
        (grouped["Dataset"] == dataset) & 
        (grouped["Budget"] == budget) & 
        (grouped["SlewRate"] == slew_rate) & 
        (grouped["DwellTime"] == dwell_time)
    ]

    methods = setting_df["Method"].tolist()
    x_positions = np.arange(len(methods)) 
    colors = [method_colors[method] for method in methods]  

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

    bars1 = ax1.bar(x_positions, setting_df["SumProb"], color=colors, alpha=0.8)
    ax1.set_ylabel("Sum Probability", color="black")
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(methods, rotation=45, ha='right')
    ax1.set_title("Sum Probability Comparison")

    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, height, f"{height:.5f}", ha='center', va='bottom', fontsize=10)


    bars2 = ax2.bar(x_positions, setting_df["TimeSec"], color=colors, alpha=0.8)
    ax2.set_ylabel("Time (sec)", color="black")
    ax2.set_xticks(x_positions)
    ax2.set_xticklabels(methods, rotation=45, ha='right')
    ax2.set_yscale("log") 
    ax2.set_title("Running Time Comparison")

    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, height, f"{height:.5f}", ha='center', va='bottom', fontsize=10)

    plt.suptitle(f"{dataset.split('/')[-1]} | Budget: {budget}, SlewRate: {slew_rate}, DwellTime: {dwell_time}")

    plt.savefig(f"comparison_{dataset.split('/')[-1]}_b{budget}_sr{slew_rate}_dt{dwell_time}.png", 
                dpi=300,
                bbox_inches="tight")
    plt.show()