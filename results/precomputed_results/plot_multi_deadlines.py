import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# Define custom markers for each method
markers = {
    'SO1': 'o',
    'SO2': 's',
    'SO3': 'D',
    'SO4': '^',
    'SO5': 'v',
    'D6': 'X',
    'D7': 'P'
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

base_betas = [1.0, 0.8, 0.5, 0.4, 0.2, 0.5, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1]

base_betas = [1.0, 1.0, 0.5, 0.4, 0.3, 0.2, 0.1, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0]
selected_indices = [1, 2, 5]
labels = ["SO1", "SO2", "SO3", "SO4", "SO5", "D6"]

# Compute segment ranges
ranges = [(0, selected_indices[0])]
for i in range(len(selected_indices) - 1):
    ranges.append((selected_indices[i] + 1, selected_indices[i + 1]))

# Select corresponding beta values
selected_beta = [base_betas[i] for i in selected_indices]


def plotMultiDeadlines():
    plt.figure(figsize=(3, 3), dpi=300)
    selected_deadlines = set()
    for idx, label in zip(selected_indices, labels):
        row = df.iloc[idx]
        fom = row["FOMs"]
        # Accumulate FoM over defined ranges
        segment_fom = [sum(fom[a:b+1]) for a, b in ranges]
        # Compute discounted cumulative FoM
        efom = []
        acc = 0.0
        for beta, f in zip(selected_beta, segment_fom):
            print(beta)
            acc += beta * f
            efom.append(acc)

        # Plot only selected deadlines
        deadlines = row["Budgets"]
        # print(deadlines)
        selected_deadlines = [deadlines[i] for i in selected_indices]
        # print(selected_deadlines)
        # selected_deadlines = [0.0] + selected_deadlines
        # efom = [0.0] + efom

        plt.plot(
            selected_deadlines,
            efom,
            label=label,
            marker=markers.get(label, 'o'),
            linewidth=1.5,
            markersize=6
        )

    # Final touches
    plt.xlabel("Deadlines (Seconds)")
    plt.ylabel("Expected FoM")
    plt.xticks(selected_deadlines)
    plt.legend()
    plt.ylim(bottom=0.03)
    plt.tight_layout()
    plt.savefig("plots/multi_expected_fom.png", dpi=300)
    plt.show()



# Load data
path = "multi_deadline"
fname = "GW200216_220804"
filepath = f"{path}/multi_filtered_{fname}_7dt.csv"
# fname = "GW200322_091133"
# filepath = f"{path}/multi_filtered_{fname}_7dt_separate.csv"
df = pd.read_csv(filepath)
df["FOMs"] = df["FOMs"].apply(lambda x: list(map(float, str(x).split())))
df["Budgets"] = df["Budgets"].apply(lambda x: list(map(float, str(x).split())))
plotMultiDeadlines()

path = "multi_deadline"
# fname = "GW200216_220804"
# filepath = f"{path}/multi_filtered_{fname}_7dt.csv"
fname = "GW200322_091133"
filepath = f"{path}/multi_filtered_{fname}_7dt_separate.csv"
df = pd.read_csv(filepath)
df["FOMs"] = df["FOMs"].apply(lambda x: list(map(float, str(x).split())))
df["Budgets"] = df["Budgets"].apply(lambda x: list(map(float, str(x).split())))
plotMultiDeadlines()
