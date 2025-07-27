# import pandas as pd
# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt

# # # # File list
# # files = [
# #     ("out_single_3.csv", "S1"),
# #     ("out_single_2.csv", "S2"),
# #     ("out_single_1.csv", "S3"),
# #     ("out_single_6.csv", "S4"),
# #     ("out_single_5.csv", "S5"),
# #     ("out_single_4.csv", "S6"),
# #     ("out_single_9.csv", "S7"),
# #     ("out_single_8.csv", "S8"),
# #     ("out_single_7.csv", "S9"),
# # ]

# # files = [
# #     ("out_single_deepslow_3.csv", "S1"),
# #     ("out_single_deepslow_2.csv", "S2"),
# #     ("out_single_deepslow_1.csv", "S3"),
# #     ("out_single_deepslow_6.csv", "S4"),
# #     ("out_single_deepslow_5.csv", "S5"),
# #     ("out_single_deepslow_4.csv", "S6"),
# #     ("out_single_deepslow_9.csv", "S7"),
# #     ("out_single_deepslow_8.csv", "S8"),
# #     ("out_single_deepslow_7.csv", "S9"),
# # ]

# files = [
#     ("out_combined_3.csv", "S1"),
#     ("out_combined_2.csv", "S2"),
#     ("out_combined_1.csv", "S3"),
#     ("out_combined_6.csv", "S4"),
#     ("out_combined_5.csv", "S5"),
#     ("out_combined_4.csv", "S6"),
#     # ("out_combined_9.csv", "S7"),
#     # ("out_combined_8.csv", "S8"),
#     ("out_combined_7.csv", "S9"),
# ]


# # files = [
# #     ("out_combined_deepslow_3.csv", "S1"),
# #     ("out_combined_deepslow_2.csv", "S2"),
# #     ("out_combined_deepslow_1.csv", "S3"),
# #     ("out_combined_deepslow_6.csv", "S4"),
# #     ("out_combined_deepslow_5.csv", "S5"),
# #     ("out_combined_deepslow_4.csv", "S6"),
# #     # ("out_combined_deepslow_9.csv", "S7"),
# #     # ("out_combined_deepslow_8.csv", "S8"),
# #     ("out_combined_deepslow_7.csv", "S9"),
# # ]



# # Methods to compare
# methods_to_compare = ["Gurobi", "Greedy", "Genetic", "AntColony", "SimulatedAnnealing", "Greedy+Genetic"]

# # Load and tag each setting
# all_data = []
# for file, setting_name in files:
#     try:
#         df = pd.read_csv(file)
#         df["setting"] = setting_name
#         all_data.append(df)
#     except FileNotFoundError:
#         print(f"Warning: File {file} not found. Skipping.")

# df_all = pd.concat(all_data)
# df_all.columns = df_all.columns.str.lower()

# # Group stats for TimeSec
# summary_time = df_all.groupby(["setting", "method"]).agg(
#     mean_time=("timesec", "mean"),  # Mean running time
#     std_time=("timesec", "std"),    # Standard deviation for CI
#     count=("timesec", "count")      # Number of samples for CI
# ).reset_index()

# # Calculate 95% CI
# summary_time["sem"] = summary_time["std_time"] / np.sqrt(summary_time["count"])
# summary_time["ci95"] = 1.96 * summary_time["sem"]

# # Filter summary to include only the specified methods
# summary_time = summary_time[summary_time["method"].isin(methods_to_compare)]

# # Check if any methods are missing
# missing_methods = set(methods_to_compare) - set(summary_time["method"])
# if missing_methods:
#     print(f"Warning: The following methods were not found in the data: {missing_methods}")

# # Sort settings
# setting_order = [s for _, s in files]
# summary_time["setting"] = pd.Categorical(summary_time["setting"], categories=setting_order, ordered=True)
# method_order = summary_time["method"].unique()

# # Line Plot with Ribbons
# plt.figure(figsize=(6.5, 4))  # Single-column journal size (adjust to 8x5 for double-column)
# sns.set_context("paper", font_scale=1.2)

# # Plot lines
# sns.lineplot(
#     data=summary_time,
#     x="setting",
#     y="mean_time",
#     hue="method",
#     hue_order=method_order,
#     marker="o",
#     palette="tab10",
#     errorbar=None,  # Disable default CI
#     linewidth=1.5  # Thicker lines for visibility
# )

# # Add 95% CI ribbons
# for method in method_order:
#     subset = summary_time[summary_time["method"] == method]
#     plt.fill_between(
#         range(len(setting_order)),
#         subset["mean_time"] - subset["ci95"],
#         subset["mean_time"] + subset["ci95"],
#         alpha=0.2
#     )

# # Apply log scale with formatted ticks
# plt.yscale('log')
# # plt.ylim(1e-3, None)  # Lower limit at 0.001 seconds, adjust as needed
# # plt.yticks([1e-3, 1e-2, 1e-1, 1, 10, 100], ['0.001', '0.01', '0.1', '1', '10', '100'], fontsize=10)

# # Customize
# # plt.title("Running Time Across Settings with 95% CI", fontsize=14, pad=15)
# # plt.xlabel("Observation Setting", fontsize=12, labelpad=5)
# plt.ylabel("Running Time (seconds)", fontsize=12, labelpad=10)
# plt.xticks(range(len(setting_order)), setting_order, rotation=45, ha='right', fontsize=11)
# plt.legend(title="Method", title_fontsize=12, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.grid(True, which="both", ls="--", alpha=0.2)  # Subtle gridlines

# # Adjust layout
# plt.tight_layout()

# # Save as high-resolution PDF
# plt.savefig("lineplot_with_ribbons_running_time_log_paper.png", dpi=300, bbox_inches='tight')
# plt.show()


# # # Group stats for TimeSec (min, max, and mean for the line)
# summary_time = df_all.groupby(["setting", "method"]).agg(
#     mean_time=("timesec", "mean"),  # Mean for the line
#     min_time=("timesec", "min"),    # Lower bound
#     max_time=("timesec", "max"),    # Upper bound
#     count=("timesec", "count")      # Number of samples (optional)
# ).reset_index()

# # Filter summary to include only the specified methods
# summary_time = summary_time[summary_time["method"].isin(methods_to_compare)]

# # Check if any methods are missing
# missing_methods = set(methods_to_compare) - set(summary_time["method"])
# if missing_methods:
#     print(f"Warning: The following methods were not found in the data: {missing_methods}")

# # Sort settings
# setting_order = [s for _, s in files]
# summary_time["setting"] = pd.Categorical(summary_time["setting"], categories=setting_order, ordered=True)
# method_order = summary_time["method"].unique()

# # Line Plot with Ribbons (Lower and Upper Bounds)
# plt.figure(figsize=(6.5, 4))  # Single-column journal size
# sns.set_context("paper", font_scale=1.2)

# # Plot lines (mean time)
# sns.lineplot(
#     data=summary_time,
#     x="setting",
#     y="mean_time",
#     hue="method",
#     hue_order=method_order,
#     marker="o",
#     palette="tab10",
#     errorbar=None,  # No default CI
#     linewidth=1.5   # Thicker lines for visibility
# )

# # Add ribbons for lower and upper bounds
# for method in method_order:
#     subset = summary_time[summary_time["method"] == method]
#     plt.fill_between(
#         range(len(setting_order)),
#         subset["min_time"],
#         subset["max_time"],
#         alpha=0.2
#     )

# # Apply log scale with formatted ticks
# plt.yscale('log')
# # plt.ylim(1e-5, None)  # Lower limit at 0.001 seconds, adjust as needed
# # plt.yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ['0.00001', '0.0001', '0.001', '0.01', '0.1', '1', '10', '100'], fontsize=10)

# # Customize
# # plt.title("Running Time Range Across Settings", fontsize=14, pad=15)
# plt.xlabel("Setting", fontsize=12, labelpad=5)
# plt.ylabel("Running Time (seconds)", fontsize=12, labelpad=10)
# plt.xticks(range(len(setting_order)), setting_order, rotation=45, ha='right', fontsize=11)
# plt.legend(title="Method", title_fontsize=12, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.grid(True, which="both", ls="--", alpha=0.2)  # Subtle gridlines

# # Adjust layout
# plt.tight_layout()

# # Save as high-resolution PDF
# plt.savefig("lineplot_with_bounds_running_time.png", dpi=300, bbox_inches='tight')
# plt.show()


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

# # # File list
# files = [
#     ("out_single_3.csv", "S1"),
#     ("out_single_2.csv", "S2"),
#     ("out_single_1.csv", "S3"),
#     ("out_single_6.csv", "S4"),
#     ("out_single_5.csv", "S5"),
#     ("out_single_4.csv", "S6"),
#     ("out_single_9.csv", "S7"),
#     ("out_single_8.csv", "S8"),
#     ("out_single_7.csv", "S9"),
# ]

# files = [
#     ("out_single_deepslow_3.csv", "S1"),
#     ("out_single_deepslow_2.csv", "S2"),
#     ("out_single_deepslow_1.csv", "S3"),
#     ("out_single_deepslow_6.csv", "S4"),
#     ("out_single_deepslow_5.csv", "S5"),
#     ("out_single_deepslow_4.csv", "S6"),
#     ("out_single_deepslow_9.csv", "S7"),
#     ("out_single_deepslow_8.csv", "S8"),
#     ("out_single_deepslow_7.csv", "S9"),
# ]

files = [
    ("out_combined_3.csv", "S1"),
    ("out_combined_2.csv", "S2"),
    ("out_combined_1.csv", "S3"),
    ("out_combined_6.csv", "S4"),
    ("out_combined_5.csv", "S5"),
    ("out_combined_4.csv", "S6"),
    # ("out_combined_9.csv", "S7"),
    # ("out_combined_8.csv", "S8"),
    ("out_combined_7.csv", "S9"),
]


# files = [
#     ("out_combined_deepslow_3.csv", "S1"),
#     ("out_combined_deepslow_2.csv", "S2"),
#     ("out_combined_deepslow_1.csv", "S3"),
#     ("out_combined_deepslow_6.csv", "S4"),
#     ("out_combined_deepslow_5.csv", "S5"),
#     ("out_combined_deepslow_4.csv", "S6"),
#     # ("out_combined_deepslow_9.csv", "S7"),
#     # ("out_combined_deepslow_8.csv", "S8"),
#     ("out_combined_deepslow_7.csv", "S9"),
# ]

# Methods to compare
methods_to_compare = ["Gurobi", "Greedy", "Genetic", "AntColony", "SimulatedAnnealing", "Greedy+Genetic"]

# Load and tag each setting
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

# Group stats for TimeSec (for 95% CI plot)
summary_time = df_all.groupby(["setting", "method"]).agg(
    mean_time=("timesec", "mean"),  # Mean running time
    std_time=("timesec", "std"),    # Standard deviation for CI
    count=("timesec", "count")      # Number of samples for CI
).reset_index()

# Calculate 95% CI
summary_time["sem"] = summary_time["std_time"] / np.sqrt(summary_time["count"])
summary_time["ci95"] = 1.96 * summary_time["sem"]

# Filter summary to include only the specified methods
summary_time = summary_time[summary_time["method"].isin(methods_to_compare)]

# Check if any methods are missing
missing_methods = set(methods_to_compare) - set(summary_time["method"])
if missing_methods:
    print(f"Warning: The following methods were not found in the data: {missing_methods}")

# Sort settings
setting_order = [s for _, s in files]
summary_time["setting"] = pd.Categorical(summary_time["setting"], categories=setting_order, ordered=True)
method_order = summary_time["method"].unique()

# Line Plot with Ribbons (95% CI)
plt.figure(figsize=(6.5, 4))  # Single-column journal size
sns.set_context("paper", font_scale=1.2)

# Plot lines
sns.lineplot(
    data=summary_time,
    x="setting",
    y="mean_time",
    hue="method",
    hue_order=method_order,
    marker="o",
    palette="tab10",
    errorbar=None,
    linewidth=1.5
)

# Add 95% CI ribbons
for method in method_order:
    subset = summary_time[summary_time["method"] == method]
    plt.fill_between(
        range(len(setting_order)),
        subset["mean_time"] - subset["ci95"],
        subset["mean_time"] + subset["ci95"],
        alpha=0.2
    )

# Apply log scale with formatted ticks
plt.yscale('log')
# Determine y-axis limits based on data
y_min = max(1e-5, summary_time["mean_time"].min() * 0.5)  # Avoid zero or negative
y_max = summary_time["mean_time"].max() * 2
plt.ylim(y_min, y_max)
# Set major and minor ticks
# plt.yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ['0.00001', '0.0001', '0.001', '0.01', '0.1', '1', '10', '100'], fontsize=10)
plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))

# Customize
plt.ylabel("Running Time (seconds)", fontsize=12, labelpad=10)
plt.xticks(range(len(setting_order)), setting_order, rotation=45, ha='right', fontsize=11)
plt.legend(title="Method", title_fontsize=12, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, which="both", ls="--", alpha=0.2)  # Show both major and minor gridlines

plt.tight_layout()
plt.savefig("lineplot_with_ribbons_running_time_log_paper.png", dpi=300, bbox_inches='tight')
plt.show()
plt.clf()  # Clear the figure to avoid state interference

# Group stats for TimeSec (min, max, and mean for the line)
summary_time = df_all.groupby(["setting", "method"]).agg(
    mean_time=("timesec", "mean"),  # Mean for the line
    min_time=("timesec", "min"),    # Lower bound
    max_time=("timesec", "max"),    # Upper bound
    count=("timesec", "count")      # Number of samples (optional)
).reset_index()

# Filter summary to include only the specified methods
summary_time = summary_time[summary_time["method"].isin(methods_to_compare)]

# Check if any methods are missing
missing_methods = set(methods_to_compare) - set(summary_time["method"])
if missing_methods:
    print(f"Warning: The following methods were not found in the data: {missing_methods}")

# Sort settings
setting_order = [s for _, s in files]
summary_time["setting"] = pd.Categorical(summary_time["setting"], categories=setting_order, ordered=True)
method_order = summary_time["method"].unique()

# Line Plot with Ribbons (Lower and Upper Bounds)
plt.figure(figsize=(6.5, 4))  # Single-column journal size
sns.set_context("paper", font_scale=1.2)

# Plot lines (mean time)
sns.lineplot(
    data=summary_time,
    x="setting",
    y="mean_time",
    hue="method",
    hue_order=method_order,
    marker="o",
    palette="tab10",
    errorbar=None,
    linewidth=1.5
)

# Add ribbons for lower and upper bounds
for method in method_order:
    subset = summary_time[summary_time["method"] == method]
    plt.fill_between(
        range(len(setting_order)),
        subset["min_time"],
        subset["max_time"],
        alpha=0.2
    )

# Apply log scale with formatted ticks
plt.yscale('log')
# Determine y-axis limits based on data
y_min = max(1e-5, summary_time["min_time"].min() * 0.5)  # Use min_time for bounds
y_max = summary_time["max_time"].max() * 2
plt.ylim(y_min, y_max)
# Set major and minor ticks
# plt.yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ['0.00001', '0.0001', '0.001', '0.01', '0.1', '1', '10', '100'], fontsize=10)
plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))

# Customize
plt.xlabel("Setting", fontsize=12, labelpad=5)
plt.ylabel("Running Time (seconds)", fontsize=12, labelpad=10)
plt.xticks(range(len(setting_order)), setting_order, rotation=45, ha='right', fontsize=11)
plt.legend(title="Method", title_fontsize=12, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, which="both", ls="--", alpha=0.2)  # Show both major and minor gridlines

plt.tight_layout()
plt.savefig("lineplot_with_bounds_running_time.png", dpi=300, bbox_inches='tight')
plt.show()
plt.clf()  # Clear the figure

