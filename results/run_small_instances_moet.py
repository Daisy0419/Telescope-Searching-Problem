import os
import subprocess
import pandas as pd

def run(dataset, budgets, run, run_times, dataset_idx, num_datasets, run_case=3):
    for i, bset in enumerate(budgets):
        orig_budget, budget_greedy, budget_genetic, budget_gcp = bset
        print(f"\nRunning ./ts with dataset={dataset} ({dataset_idx+1} of {num_datasets}), budget={orig_budget} ({i+1} of {len(budgets)}), (run {run+1} of {run_times}) run_case={run_case}")
        cmd = [EXECUTABLE, dataset, str(run_case), str(orig_budget), str(budget_greedy), str(budget_genetic), str(budget_gcp)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout
            print(output)

        except subprocess.CalledProcessError as e:
            print(f"Error running command: {cmd}")
            print(f"Error message: {e.stderr}")

if __name__ == "__main__":
    EXECUTABLE = "../build/ts"
    data_path = "../data"

    skymaps = [
        "GW200208_222617_7dt",
        "GW200216_220804_7dt",
        "GW200220_124850_7dt",
        "GW191113_071753_7dt",
        "GW200306_093714_7dt"
    ]

    output_dir = os.path.join("..", "results", "recomputed_results", "instances_with_moet")
    os.makedirs(output_dir, exist_ok=True)
    default_name = os.path.join(output_dir, "out.csv")

    budgets = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # run each map each budget case for 5 time to compute moet
    run_times = 5
    output_sub_dir0 = os.path.join(output_dir, f"get_moet")
    os.makedirs(output_sub_dir0, exist_ok=True)
    for i, skymap in enumerate(skymaps):
        print(f"\nProcessing skymap: {skymap}")
        updated_budgets = []
        for budget in budgets:
            updated_budget = [budget]
            for method in ["Greedy", "Genetic", "Hoogeveen"]: 
                updated_budget.append(budget) # not account moet 
            updated_budgets.append(updated_budget)

        dataset = os.path.join(data_path, f"large_moet/filtered_{skymap}.csv")
        print(f"Using dataset: {dataset}")

        for runtimes in range(run_times):
            run(dataset, updated_budgets, run=runtimes, run_times=run_times, dataset_idx=i, num_datasets=len(skymaps))

        new_name = os.path.join(output_sub_dir0, f"out_{skymap}.csv")
        if os.path.exists(default_name):
            os.rename(default_name, new_name)
            print(f"Saved result to {new_name}")
        else:
            print(f"Warning: {default_name} not found after running {skymap}")


    # compute moet for each map each case
    output_sub_dir1 = os.path.join(output_dir, f"moet")
    os.makedirs(output_sub_dir1, exist_ok=True)
    for skymap in skymaps:
        file =f"{output_sub_dir0}/out_{skymap}.csv"
        df = pd.read_csv(file)
        # Compute max runtime per (Method, Budget)
        max_runtime = df.groupby(['Method', 'Budget'])['TimeSec'].max().reset_index()
        runtime_table = max_runtime.pivot(index='Budget', columns='Method', values='TimeSec')
        max_runtime.to_csv(f"{output_sub_dir1}/{skymap}.csv")
    


    # incorporate moet to budget and rerun the cases
    output_sub_dir2 = os.path.join(output_dir, f"result_with_moet")
    os.makedirs(output_sub_dir2, exist_ok=True)
    for i, skymap in enumerate(skymaps):
        print(f"\nProcessing skymap: {skymap}")
        updated_budgets = []
        df = pd.read_csv(f"{output_sub_dir1}/{skymap}.csv")

        for budget in budgets:
            updated_budget = [budget]
            for method in ["Greedy", "Genetic", "Hoogeveen"]: 
                values = df.loc[(df["Method"] == method) & (df["Budget"] == budget), "TimeSec"].values
                value = float(values[0]) if len(values) > 0 else 0.0
                adjusted = max(0.0, budget - value)
                updated_budget.append(adjusted)
            updated_budgets.append(updated_budget)

        dataset = os.path.join(data_path, f"large_moet/filtered_{skymap}.csv")
        print(f"Using dataset: {dataset}")

        run(dataset, updated_budgets, run=1, run_times = 1, dataset_idx = i, num_datasets = len(skymaps))

        new_name = os.path.join(output_sub_dir2, f"out_{skymap}.csv")
        if os.path.exists(default_name):
            os.rename(default_name, new_name)
            print(f"Saved result to {new_name}")
        else:
            print(f"Warning: {default_name} not found after running {skymap}")
