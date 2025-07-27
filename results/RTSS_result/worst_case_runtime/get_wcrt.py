import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

# Display name mapping
method_display_names = {
    'Hoogeveen': 'GCP',
    'Greedy': 'Greedy',
    'Gurobi': 'ILP',
    'Genetic': 'Genetic',
}

# Load the CSV file (choose one)
# df = pd.read_csv("wcrt_filtered_GW200208_222617_7dt.csv") 
# df = pd.read_csv("wcrt_filtered_GW200216_220804_7dt.csv") 
df = pd.read_csv("wcrt_filtered_GW200220_124850_7dt.csv") 
df = pd.read_csv("wcrt_filtered_GW191113_071753_7dt.csv") 
df = pd.read_csv("wcrt_filtered_GW200306_093714_7dt.csv") 

# Compute max runtime per (Method, Budget)
max_runtime = df.groupby(['Method', 'Budget'])['TimeSec'].max().reset_index()
runtime_table = max_runtime.pivot(index='Budget', columns='Method', values='TimeSec')


mean_runtime = df.groupby(['Method', 'Budget'])['TimeSec'].mean().reset_index()
meantime_table = mean_runtime.pivot(index='Budget', columns='Method', values='TimeSec')

# Candidate functions for fitting
def f_log(x, a, b):         return a * np.log(x) + b
def f_linear(x, a, b):      return a * x + b
def f_quad(x, a, b, c):     return a * x**2 + b * x + c
def f_cubic(x, a, b, c, d): return a * x**3 + b * x**2 + c * x + d

models = {
    # "log": f_log,
    "linear": f_linear,
    "quadratic": f_quad,
    "cubic": f_cubic
}


# Fit and plot per method
for method in runtime_table.columns:
    x = runtime_table.index.values
    y = runtime_table[method].values

    # Clean data: remove NaNs
    mask = ~np.isnan(y)
    x_clean = x[mask]
    y_clean = y[mask]

    fits = {}
    rss_scores = {}

    for name, func in models.items():
        try:
            popt, _ = curve_fit(func, x_clean, y_clean, maxfev=10000)
            y_pred = func(x_clean, *popt)
            rss = np.sum((y_clean - y_pred) ** 2)
            fits[name] = (popt, y_pred)
            rss_scores[name] = rss
        except Exception as e:
            print(f"Failed to fit {name} for {method}: {e}")

    if rss_scores:
        best_model = min(rss_scores, key=rss_scores.get)
        print(f"[{method}] Best model: {best_model} (RSS = {rss_scores[best_model]:.4g})")
    else:
        print(f"[{method}] No valid fit.")
        continue

    # Plot results
    plt.figure(figsize=(6, 3))
    plt.scatter(x_clean, y_clean, color='black', label='Data')

    for name, (popt, y_pred) in fits.items():
        plt.plot(x_clean, y_pred, label=f"{name} fit")

    display_name = method_display_names.get(method, method)
    plt.xlabel("Budget")
    plt.ylabel("WCRT")
    plt.title(f"{display_name} (Best fit: {best_model})")
    plt.legend()
    plt.tight_layout()
    plt.show()


# Load the CSV file (choose one)
# df = pd.read_csv("wcrt_filtered_GW200208_222617_7dt.csv") 
# df = pd.read_csv("wcrt_filtered_GW200216_220804_7dt.csv") 
# df = pd.read_csv("wcrt_filtered_GW200220_124850_7dt.csv") 
# df = pd.read_csv("wcrt_filtered_GW191113_071753_7dt.csv") 
# df = pd.read_csv("wcrt_filtered_GW200306_093714_7dt.csv") 


# # Fit and plot per method
# for method in meantime_table.columns:
#     x = meantime_table.index.values
#     y = meantime_table[method].values

#     # Clean data: remove NaNs
#     mask = ~np.isnan(y)
#     x_clean = x[mask]
#     y_clean = y[mask]

#     fits = {}
#     rss_scores = {}

#     for name, func in models.items():
#         try:
#             popt, _ = curve_fit(func, x_clean, y_clean, maxfev=10000)
#             y_pred = func(x_clean, *popt)
#             rss = np.sum((y_clean - y_pred) ** 2)
#             fits[name] = (popt, y_pred)
#             rss_scores[name] = rss
#         except Exception as e:
#             print(f"Failed to fit {name} for {method}: {e}")

#     if rss_scores:
#         best_model = min(rss_scores, key=rss_scores.get)
#         print(f"[{method}] Best model: {best_model} (RSS = {rss_scores[best_model]:.4g})")
#     else:
#         print(f"[{method}] No valid fit.")
#         continue

#     # Plot results
#     plt.figure(figsize=(6, 3))
#     plt.scatter(x_clean, y_clean, color='black', label='Data')

#     for name, (popt, y_pred) in fits.items():
#         plt.plot(x_clean, y_pred, label=f"{name} fit")

#     display_name = method_display_names.get(method, method)
#     plt.xlabel("Budget")
#     plt.ylabel("WCRT")
#     plt.title(f"{display_name} (Best fit: {best_model})")
#     plt.legend()
#     plt.tight_layout()
#     plt.show()


maps = ["GW200208_222617_7dt", "GW200216_220804_7dt", "GW200220_124850_7dt",
        "GW191113_071753_7dt", "GW200306_093714_7dt"]

for map in maps:
    file =f"wcrt_test2/wcrt_filtered_{map}.csv"
    df = pd.read_csv(file)
    # Compute max runtime per (Method, Budget)
    max_runtime = df.groupby(['Method', 'Budget'])['TimeSec'].max().reset_index()
    runtime_table = max_runtime.pivot(index='Budget', columns='Method', values='TimeSec')
    max_runtime.to_csv(f"get_wcrt2/{map}.csv")
    # print(max_runtime)

maps = [
    "GW200208_222617_7dt",
    "GW200216_220804_7dt",
    "GW200220_124850_7dt",
    "GW191113_071753_7dt",
    "GW200306_093714_7dt"
]

all_runtimes = []

for map_name in maps:
    file = f"wcrt_test2/wcrt_filtered_{map_name}.csv"
    df = pd.read_csv(file)
    
    # Compute max runtime per (Method, Budget)
    max_runtime = df.groupby(['Method', 'Budget'])['TimeSec'].max().reset_index()
    all_runtimes.append(max_runtime)

# Combine all max_runtime DataFrames
combined_df = pd.concat(all_runtimes)

# Compute max runtime across all maps for each (Method, Budget)
global_max = combined_df.groupby(['Method', 'Budget'])['TimeSec'].max().reset_index()

# Optional: Pivot into table format
max_runtime_table = global_max.pivot(index='Budget', columns='Method', values='TimeSec')

# Save to CSV
max_runtime_table.to_csv("max_runtime_across_maps.csv")

# Display result
print(max_runtime_table)