import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

# Load the max runtime table
df = pd.read_csv("mean_runtime_table.csv", index_col=0)

# Candidate functions
def f_log(x, a, b):       return a * np.log(x) + b
def f_linear(x, a, b):    return a * x + b
def f_quad(x, a, b, c):   return a * x**2 + b * x + c
def f_cubic(x, a, b, c, d): return a * x**3 + b * x**2 + c * x + d

models = {
    # "log": f_log,
    "linear": f_linear,
    "quadratic": f_quad,
    "cubic": f_cubic
}

# Loop over methods
for method in df.columns:
    x = df.index.to_numpy()
    y = df[method].to_numpy()

    # Fit each model
    fits = {}
    rss_scores = {}

    for name, func in models.items():
        try:
            popt, _ = curve_fit(func, x, y)
            y_pred = func(np.array(x), *popt)
            rss = np.sum((y - y_pred) ** 2)
            fits[name] = (popt, y_pred)
            rss_scores[name] = rss
        except Exception as e:
            print(f"Failed to fit {name} for {method}: {e}")

    # Select best model
    if rss_scores:
        best_model = min(rss_scores, key=rss_scores.get)
        print(f"[{method}] Best model: {best_model} (RSS = {rss_scores[best_model]:.4g})")
    else:
        print(f"[{method}] No valid fit.")
        continue

    # Plot
    plt.figure(figsize=(6, 3))
    plt.scatter(x, y, color='black', label='Data')
    for name, (popt, y_pred) in fits.items():
        plt.plot(x, y_pred, label=f"{name} fit")
    plt.xlabel("Budget")
    plt.ylabel("WCRT")
    plt.title(f"{method} (Best: {best_model})")
    # plt.title(f"{method} (Best: {best_model})")
    plt.legend()
    plt.tight_layout()
    plt.show()
