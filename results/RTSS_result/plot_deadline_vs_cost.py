# import pandas as pd
# import matplotlib.pyplot as plt
# from pathlib import Path

# # 1)  Load the result table  ──────────────────────────────────────────────
# #
# #    If you already have the file on disk just point to it; otherwise
# #    paste the block you showed into a CSV called “results.csv”.
# #
# DATA_FILE = Path("large_with_Gurobi2/out_filtered_GW191103_012549_7dt.csv")          # <-- adjust if needed
# cols      = ["Method", "Budget", "TotalCost"]   # only the fields we need

# df = pd.read_csv(DATA_FILE, usecols=cols)

# # 2)  Make a scatter-plot of  Budget  vs  TotalCost  ──────────────────────
# #
# #    • Each method gets its own marker/colour from matplotlib’s default
# #      cycle (no explicit colours set, as requested).
# #    • A dashed y = x reference line shows the “exactly on budget” limit.
# #
# fig, ax = plt.subplots(figsize=(7, 7))

# for method, sub in df.groupby("Method"):
#     ax.scatter(sub["Budget"], sub["TotalCost"],
#                label=method, s=60, alpha=0.8)

# # reference line
# lims = [0, df["Budget"].max()*1.05]      # a little head-room at the top
# ax.plot(lims, lims, ls="--", lw=1)       # y = x
# ax.set_xlim(lims)
# ax.set_ylim(lims)

# ax.set_xlabel("Budget")
# ax.set_ylabel("Total cost actually incurred")
# ax.set_title("Every method remains well under budget")
# ax.legend(title="Method", loc="upper left")
# plt.tight_layout()
# plt.show()



import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────
# 1.  Load *all* result-files in the folder
# ──────────────────────────────────────────────────────────────────────────
FOLDER = Path("large_with_Gurobi2")                      # <-- adjust if needed
COLS   = ["Method", "Budget", "TotalCost"]

df = (
    pd.concat(
        (pd.read_csv(f, usecols=COLS) for f in FOLDER.glob("*.csv")),
        ignore_index=True
    )
)

# ──────────────────────────────────────────────────────────────────────────
# 2.  Summarise:  for each (Budget, Method) compute min / max / mean cost
# ──────────────────────────────────────────────────────────────────────────
g = df.groupby(["Budget", "Method"])["TotalCost"]
summ = g.agg(["min", "max", "mean"]).reset_index()

# ──────────────────────────────────────────────────────────────────────────
# 3.  Plot:          mean  ±  (max–min)/2   as symmetric error-bars
# ──────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

for method, sub in summ.groupby("Method"):
    x    = sub["Budget"]
    y    = sub["mean"]
    yerr = [(y - sub["min"]), (sub["max"] - y)]
    ax.errorbar(x, y, yerr=yerr,
                fmt="o",  lw=1, capsize=4,  # marker & tiny caps
                label=method, alpha=.9)

# y = x reference (exactly-on-budget)
lims = [0, summ["Budget"].max()*1.05]
ax.plot(lims, lims, ls="--", lw=1)
ax.set_xlim(lims)
ax.set_ylim(lims)

ax.set_xlabel("Budget (deadline)")
ax.set_ylabel("Total cost")
ax.set_title("Cost range for each budget – all runs, all files")
ax.legend(title="Method")
plt.tight_layout()
plt.show()


df["Slack"] = df["Budget"] - df["TotalCost"]

# 2️⃣  summarise per (Budget, Method)
slack = (
    df.groupby(["Budget", "Method"])["Slack"]
      .agg(["mean", "min", "max"])
      .reset_index()
)

# 3️⃣  plot: mean line + shaded band
fig, ax = plt.subplots(figsize=(8, 6))

for method, sub in slack.groupby("Method"):
    sub = sub.sort_values("Budget")          # make sure x is ordered
    x   = sub["Budget"]
    ax.plot(x, sub["mean"], lw=2, label=method)                # mean
    ax.fill_between(                                           # band
        x,
        sub["min"],
        sub["max"],
        alpha=0.20,
        linewidth=0,
    )

# 4️⃣  cosmetics
ax.set_xlabel("Budget (deadline)")
ax.set_ylabel("Slack  (Budget – TotalCost)")
ax.set_ylim(-0.1,0.1)
ax.set_title("Slack range per budget")
ax.axhline(0, color="k", lw=1, ls="--")      # zero-slack reference
ax.legend(title="Method")
plt.tight_layout()
plt.show()