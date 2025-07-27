import os
import subprocess
import pandas as pd

def run(dataset, budgets, run_case=3):
    for bset in budgets:
        orig_budget, budget_greedy, budget_genetic, budget_gcp = bset
        print(f"\nRunning ./ts with dataset={dataset}, budget={orig_budget}, run_case={run_case}")
        cmd = [EXECUTABLE, dataset, str(run_case), str(orig_budget), str(budget_greedy), str(budget_genetic), str(budget_gcp)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout
            print(output)

        except subprocess.CalledProcessError as e:
            print(f"Error running command: {cmd}")
            print(f"Error message: {e.stderr}")

if __name__ == "__main__":
    EXECUTABLE = "..build/ts"
    data_path = "../data"

    skymaps = [
        "GW200208_222617_7dt",
        "GW200216_220804_7dt",
        "GW200220_124850_7dt",
        "GW191113_071753_7dt",
        "GW200306_093714_7dt"
    ]

    output_dir = os.path.join("..", "results", "recomputed_results", "small_wcet")
    os.makedirs(output_dir, exist_ok=True)
    default_name = os.path.join(output_dir, "out.csv")


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

        run(dataset, updated_budgets)

        new_name = os.path.join(output_dir, f"out_{skymap}.csv")
        if os.path.exists(default_name):
            os.rename(default_name, new_name)
            print(f"Saved result to {new_name}")
        else:
            print(f"Warning: {default_name} not found after running {skymap}")
