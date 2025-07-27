import pandas as pd
import csv
import math

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
                costs[i][j] = 100000.0  # Large cost for self-loops
            else:
                dist_deg = angular_distance(ra[i], dec[i], ra[j], dec[j])
                cost_ij = dist_deg / telescope_speed + dwell_time
                costs[i][j] = cost_ij
                costs[j][i] = cost_ij 

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


filename = '../../data_RTSS/filtered_GW191219_163120.fits_7dt.csv'
costs, probability = read_data(filename, 50, 1)
print(f"costs[3][4]: {costs[3][4]}")
print(f"costs[17][18]: {costs[18][19]}")

# def fix_results(csv_file):
#     df = pd.read_csv(csv_file, quotechar='"', skipinitialspace=True)

#     for idx, row in df.iterrows():
#         dataset = row["Dataset"]
#         method = row["Method"]
#         budget = row["Budget"]
#         slew_rate = row["SlewRate"]
#         dwell_time = row["DwellTime"]
#         expected_cost = row["TotalCost"]

#         if isinstance(row["Path"], str) and row["Path"].strip():
#             path_nodes = list(map(int, row["Path"].split()))
#         else:
#             print(f"Not Valid: {method} | Dataset: {dataset} | Budget: {budget}, SlewRate: {slew_rate}, DwellTime: {dwell_time} (Empty or Invalid Path)")
#             continue

#         # if min(path_nodes) == 1:
#         #     path_nodes = [x - 1 for x in path_nodes]

#         dataset_path = f"../{dataset}"
#         if dataset == "deep_slow.csv":
#             costs, probability = read_data_deep_slow(dataset_path, slew_rate, dwell_time)
#         else:
#             costs, probability = read_data(dataset_path, slew_rate, dwell_time)

#         # if any(node >= len(costs) or node < 0 for node in path_nodes):
#         #     print(f"Not Valid: {method} | Dataset: {dataset} | Budget: {budget}, SlewRate: {slew_rate}, DwellTime: {dwell_time} (Index out of bounds)")
#         #     continue

#         actual_cost = sum(costs[path_nodes[i]][path_nodes[i+1]] for i in range(len(path_nodes)-1))
#         actual_prize = sum(probability[node] for node in path_nodes)

#         df.at[idx, "TotalCost"] = actual_cost
#         df.at[idx, "SumProb"] = actual_prize

#     df.to_csv("out_fixed.csv", index=False)

# fix_results("out_fixed.csv")


