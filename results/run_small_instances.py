import subprocess
import os

def run(budgets, dataset):
    for budget in budgets:

        print(f"\nRunning ./ts with dataset={dataset}, budget={budget}")
        cmd = [EXECUTABLE, dataset, str(budget)]
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
    EXECUTABLE = "/home/research/w.yanwang/Telescope-Searching-Problem/build_again/ts" 
    data_path = "/home/research/w.yanwang/Telescope-Searching-Problem/data_RTSS_0518/large"

    default_name = "out.csv"
    # default_name2 = "out2.csv"
    skymaps=getFiles(data_path)
    # budgets = [10, 20, 30, 40, 50, 60, 70, 80, 90]

    # budgets = [10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800,
    #            850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600]

    # budgets = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600]
    budgets = [10, 200, 400, 600, 800, 1000, 1200, 1400, 1600]

    # budgets=[0]
    for skymap in skymaps:
        # if skymap in ["filtered_GW191103_012549_7dt", "filtered_GW191109_010717_7dt", "filtered_GW191113_071753_7dt"
        # "filtered_GW191126_115259_7dt", "filtered_GW191127_050227_7dt", "filtered_GW191126_115259_7dt"]:
        #     continue
        dataset = os.path.join(data_path, f"{skymap}.csv")
        print(f"dataset: {dataset}")
        run(budgets, dataset)
        new_name = f"out_{skymap}.csv"
        os.rename(default_name, new_name)
        # new_name2 = f"multi_{skymap}.csv"
        # os.rename(default_name2, new_name2)
    
    