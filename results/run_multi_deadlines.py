run_all = True

import subprocess
import os
import sys

def run(dataset, dataset_idx, num_datasets, run_case=4):
    print(f"\nRunning ./ts with dataset={dataset} ({dataset_idx+1} of {num_datasets}), run_case={run_case}")
    cmd = [EXECUTABLE, dataset, str(run_case)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error message: {e.stderr}")
        sys.exit(1)

def getFiles(folder_path):
    return [
        os.path.splitext(f)[0]
        for f in os.listdir(folder_path)
        if f.endswith(".csv")
    ]

if __name__ == "__main__":
    EXECUTABLE = "../build/ts"
    data_path = "../data/large"
    output_dir = os.path.join("..", "results", "recomputed_results", "multi_deadlines")

    os.makedirs(output_dir, exist_ok=True)

    default_name1 = os.path.join(output_dir, "out.csv")
    default_name2 = os.path.join(output_dir, "out2.csv")

    skymaps = []
    if run_all:
        for filename in os.listdir(folder_path):
            if filename.endswith(".csv"):
                skymap = os.path.splitext(filename)[0]
                skymaps.append(skymap)
    else:
        skymaps = [
            "filtered_GW200322 091133_7dt.csv",
            "filtered_GW200216 220804_7dt.csv",
        ]

    skymaps = getFiles(data_path)

    for i, skymap in enumerate(skymaps):
        dataset = os.path.join(data_path, f"{skymap}.csv")
        print(f"Dataset: {dataset}")
        run(dataset, dataset_idx = i, num_datasets = len(skymaps))

        renamed_1 = os.path.join(output_dir, f"out_{skymap}.csv")
        renamed_2 = os.path.join(output_dir, f"multi_{skymap}.csv")

        if os.path.exists(default_name1):
            os.rename(default_name1, renamed_1)
        else:
            print(f"[WARN] Missing output: {default_name1}")

        if os.path.exists(default_name2):
            os.rename(default_name2, renamed_2)
        else:
            print(f"[WARN] Missing output: {default_name2}")

    