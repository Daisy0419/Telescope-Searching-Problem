import subprocess
import os
import sys

def run(dataset, budgets, dataset_idx, num_datasets, run_case=1):
    for i, budget in enumerate(budgets):
        print(f"\nRunning ./ts with dataset={dataset} ({dataset_idx+1} of {num_datasets}), budget={budget} ({i+1} of {len(budgets)}), run_case={run_case}")
        cmd = [EXECUTABLE, dataset, str(run_case), str(budget)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout
            print(output)

        except subprocess.CalledProcessError as e:
            print(f"Error running command: {cmd}")
            print(f"Error message: {e.stderr}")
            sys.exit(1)


def getFiles(folder_path):
    skymaps = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            skymap = os.path.splitext(filename)[0]
            skymaps.append(skymap)
    return skymaps

if __name__ == "__main__":
    EXECUTABLE = "../build/ts" 
    data_path = "../data/small"
    os.makedirs("recomputed_results/small", exist_ok=True)
    default_name = "recomputed_results/small/out.csv"
    budgets = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    skymaps=getFiles(data_path)
    for i, skymap in enumerate(skymaps):
        dataset = os.path.join(data_path, f"{skymap}.csv")
        print(f"dataset: {dataset}")
        run(dataset, budgets, dataset_idx = i, num_datasets = len(skymaps))
        new_name = f"recomputed_results/small/out_{skymap}.csv"
        if os.path.exists(default_name):
            os.rename(default_name, new_name)
        else:
            print(f"[WARN] Expected output file not found: {default_name}")
    