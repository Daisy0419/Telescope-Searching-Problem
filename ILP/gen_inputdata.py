import math
import csv
import os

def angular_distance(ra1, dec1, ra2, dec2):
    return math.degrees(math.acos(
        math.sin(math.radians(dec1)) * math.sin(math.radians(dec2)) +
        math.cos(math.radians(dec1)) * math.cos(math.radians(dec2)) *
        math.cos(math.radians(ra1 - ra2))
    ))

def read_data(filename, telescope_speed, dwell_time):
    ra, dec, probability = [], [], []

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader) 

        for row in reader:
            ra.append(float(row[2]))
            dec.append(float(row[3]))
            probability.append(float(row[4]))

    n = len(ra)
    costs = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                costs[i][j] = 100000.0
            else:
                dist_deg = angular_distance(ra[i], dec[i], ra[j], dec[j])
                cost_ij = dist_deg / telescope_speed + dwell_time
                costs[i][j] = cost_ij
                costs[j][i] = cost_ij  # Symmetric

    return costs, probability


def read_data_deep_slow(filename, telescope_speed, dwell_time):
    ra, dec, probability = [], [], []

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header

        for row in reader:
            ra.append(float(row[2]))
            dec.append(float(row[3]))
            probability.append(float(row[4]))

    n = len(ra)
    costs = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                costs[i][j] = 100000.0
            else:
                travel_time = max(abs(ra[j] - ra[i]), abs(dec[j] - dec[i])) / telescope_speed
                costs[i][j] = travel_time + dwell_time
                costs[j][i] = costs[i][j]

    return costs, probability


def write_ampl_data_as_matrix(output_filename, costs, probability, start=1, budget=15):
    n = len(probability)

    with open(output_filename, 'w') as file:
        file.write(f"# A small set of {n} nodes\n")
        file.write("set NODES := ")
        file.write(" ".join(str(i + 1) for i in range(n)))
        file.write(";\n\n")

        file.write(f"param start := {start};\n\n")

        file.write("# Prize for visiting each node\n")
        file.write("param Prize :=\n")
        for i in range(n):
            file.write(f"{i+1} {probability[i]:.6f}\n")
        file.write(";\n\n")

        file.write("# Travel costs between nodes (matrix form)\n")
        file.write("param Cost : ")
        file.write(" ".join(str(j+1) for j in range(n)))
        file.write(" :=\n")

        for i in range(n):
            row_str = " ".join(f"{costs[i][j]:.6f}" for j in range(n))
            file.write(f"{i+1}  {row_str}\n")
        file.write(";\n\n")

        file.write("# Maximum allowable budget\n")
        file.write(f"param Budget := {budget};\n")

def generateData():
    input_file = "../data/7dt_separate.csv"   
    output_file = "ampl_input_matrix.dat"     

    telescope_speed = 50.0
    dwell_time = 1.0

    costs, probability = read_data(input_file, telescope_speed, dwell_time)

    write_ampl_data_as_matrix(
        output_filename=output_file,
        costs=costs,
        probability=probability,
        start=1,
        budget=19
    )
    print(f"AMPL data (matrix form) written to {output_file}")

def batchGenerateData(is_deepSlow = False):
    for dataset in datasets:
        for budget, slew_rate, dwell_time in zip(budgets, slew_rates, dwell_times):
            if not is_deepSlow:
                costs, probability = read_data(dataset, slew_rate, dwell_time)
            else:
                costs, probability = read_data_deep_slow(dataset, slew_rate, dwell_time)

            dataset_name = os.path.basename(dataset)
            dataset_name = os.path.splitext(dataset_name)[0]

            output_file = f"{dataset_name}_b{budget}_sr{slew_rate}_dt{dwell_time}.dat"
            write_ampl_data_as_matrix(
                output_filename=output_file,
                costs=costs,
                probability=probability,
                start=1,
                budget=budget-dwell_time
            )
            print(f"AMPL data (matrix form) written to {output_file}")


if __name__ == "__main__":
    # # Define the parameter grid
    # datasets = ["../data/7dt_combined.csv", "../data/7dt_separate.csv"] 
    # budgets =     [50,  30, 30, 50, 200, 1000]
    # slew_rates =  [0.5, 1,  2,   3,  5,   10]
    # dwell_times = [0,   0,  0.5, 1,  10,  60]
    # batchGenerateData()

    # Define the parameter grid
    datasets = ["../data/deep_slow.csv"] 
    # budgets =     [200, 500, 1000, 2000, 1000, 2000, 3000, 5000]
    # slew_rates =  [1,   1,   1,    1,    1,    1,    1,    1]
    # dwell_times = [1,   10,  10,   10,   60,   60,   60,   60]
    budgets =     [50, 100, 200, 500, 1000, 1500]
    slew_rates =  [1,   1,   1,   1,   1,    1]  
    dwell_times = [5,   10,   15,  30,  60,   90] 
    batchGenerateData(True)

    # datasets = ["../data/7dt_combined.csv", "../data/7dt_separate.csv"] # , "../data/deep_slow.csv"
    # budgets =     [20, 10, 15,  20, 30, 200, 1000]
    # slew_rates =  [20, 50, 50,  50, 50,  50, 50]
    # dwell_times = [0,  0,  0.5, 1,  1, 10,  60]
    # batchGenerateData()

    # datasets = ["../data/deep_slow.csv"]
    # budgets =     [300, 500, 1000,  1500, 1500, 2000]
    # slew_rates =  [1,   1,   1,     1,    1,    1]
    # dwell_times = [0,   1,   1,     1,    10,   60]
    # batchGenerateData()
    
    # generateData()

