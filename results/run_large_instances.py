import subprocess
import os

def run(dataset, run_case=2):
    print(f"\nRunning ./ts with dataset={dataset}, run_case={run_case}")
    cmd = [EXECUTABLE, dataset, str(run_case)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output = result.stdout
        print(output)

    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error message: {e.stderr}")


def getFiles(folder_path):
    skymaps = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            skymap = os.path.splitext(filename)[0]
            skymaps.append(skymap)
    return skymaps

if __name__ == "__main__":
    EXECUTABLE = "../build/ts" 
    data_path = "../data/large"

    default_name = "out.csv"
    skymaps=getFiles(data_path)
    for skymap in skymaps:
        dataset = os.path.join(data_path, f"{skymap}.csv")
        print(f"dataset: {dataset}")
        run(dataset)
        new_name = f"out_{skymap}.csv"
        os.rename(default_name, new_name)
    