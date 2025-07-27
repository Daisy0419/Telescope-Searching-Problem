import subprocess
import os
import pandas as pd

# def run(budgets, dataset):
#     for budget in budgets:

#         print(f"\nRunning ./ts with dataset={dataset}, budget={budget}")
#         cmd = [EXECUTABLE, dataset, str(budget)]
#         try:
#             result = subprocess.run(cmd, capture_output=True, text=True, check=True)
#             output = result.stdout
#             print(output)

#         except subprocess.CalledProcessError as e:
#             print(f"Error running command: {cmd}")
#             print(f"Error message: {e.stderr}")


# import os
# import subprocess
# import pandas as pd

# def run(budgets, dataset):
#     for budget_ in budgets:
#         budget, budget_greedy, budget_genetic, budget_gcp = budget_
#         print(f"\nRunning ./ts with dataset={dataset}, budget={budget_}")
#         cmd = [EXECUTABLE, dataset, str(budget), str(budget_greedy), str(budget_genetic), str(budget_gcp)]
#         try:
#             result = subprocess.run(cmd, capture_output=True, text=True, check=True)
#             print(result.stdout)
#         except subprocess.CalledProcessError as e:
#             print(f"Error running command: {cmd}")
#             print(f"Error message: {e.stderr}")

# if __name__ == "__main__":
#     EXECUTABLE = "/home/research/w.yanwang/Telescope-Searching-Problem/build/ts"
#     data_path = "/home/research/w.yanwang/Telescope-Searching-Problem/data_RTSS_0518"

#     skymaps = [
#         "GW200208_222617_7dt",
#         "GW200216_220804_7dt",
#         "GW200220_124850_7dt",
#         "GW191113_071753_7dt",
#         "GW200306_093714_7dt"
#     ]

#     default_name = "out.csv"
#     budgets = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]

#     for skymap in skymaps:
#         updated_budgets = []
#         df = pd.read_csv(f"{data_path}/wcrt/{skymap}.csv")
        
#         for budget in budgets:
#             updated_budget = []
#             updated_budget.append(budget)
#             for method in ["Greedy", "Genetic", "Hoogeveen"]: 
#                 values = df.loc[(df["Method"] == method) & (df["Budget"] == budget), "TimeSec"].values
#                 value = float(values[0]) if len(values) > 0 else 0.0
#                 value = max(0.0, budget - value)
#                 updated_budget.append(value)
#             updated_budgets.append(updated_budget)

#         dataset = os.path.join(data_path, f"large_wcrt/filtered_{skymap}.csv") 
#         print(f"dataset: {dataset}")

#         run(updated_budgets, dataset)

#         new_name = f"out_{skymap}.csv"
#         os.rename(default_name, new_name)

    

import os
import subprocess
import pandas as pd

def run(budgets, dataset):
    for budget_ in budgets:
        budget, budget_greedy, budget_genetic, budget_gcp = budget_
        print(f"\nRunning ./ts with dataset={dataset}, budget={budget_}")
        cmd = [EXECUTABLE, dataset, str(budget), str(budget_greedy), str(budget_genetic), str(budget_gcp)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {cmd}")
            print(f"Error message: {e.stderr}")

if __name__ == "__main__":
    EXECUTABLE = "/home/research/w.yanwang/Telescope-Searching-Problem/build/ts"
    data_path = "/home/research/w.yanwang/Telescope-Searching-Problem/data_RTSS_0518"

    skymaps = [
        "GW200208_222617_7dt",
        "GW200216_220804_7dt",
        "GW200220_124850_7dt",
        "GW191113_071753_7dt",
        "GW200306_093714_7dt"
    ]

    default_name = "out.csv"
    # budgets = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
    budgets = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    for skymap in skymaps:
        print(f"\nProcessing skymap: {skymap}")
        updated_budgets = []
        df = pd.read_csv(f"{data_path}/wcrt2/{skymap}.csv")

        for budget in budgets:
            updated_budget = [budget]
            for method in ["Greedy", "Genetic", "Hoogeveen"]: 
                values = df.loc[(df["Method"] == method) & (df["Budget"] == budget), "TimeSec"].values
                value = float(values[0]) if len(values) > 0 else 0.0
                adjusted = max(0.0, budget - value)
                updated_budget.append(adjusted)
            updated_budgets.append(updated_budget)

        dataset = os.path.join(data_path, f"large_wcrt1/filtered_{skymap}.csv")
        print(f"Using dataset: {dataset}")

        run(updated_budgets, dataset)

        new_name = f"out_{skymap}.csv"
        if os.path.exists(default_name):
            os.rename(default_name, new_name)
            print(f"Saved result to {new_name}")
        else:
            print(f"Warning: {default_name} not found after running {skymap}")
