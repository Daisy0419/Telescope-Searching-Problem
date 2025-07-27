import pandas as pd
import csv
import math
import matplotlib.pyplot as plt

def angular_distance(ra1, dec1, ra2, dec2):
    return math.degrees(math.acos(
        math.sin(math.radians(dec1)) * math.sin(math.radians(dec2)) +
        math.cos(math.radians(dec1)) * math.cos(math.radians(dec2)) *
        math.cos(math.radians(ra1 - ra2))
    ))

def compute_prizes_costs(df, slew_rate, dwell_time):
    ra = df['RA'].values
    dec = df['Dec'].values
    probability = df['Probability'].values

    n = len(ra)
    costs = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                costs[i][j] = dwell_time
            else:
                dist_deg = angular_distance(ra[i], dec[i], ra[j], dec[j])
                cost_ij = dist_deg / slew_rate + dwell_time
                costs[i][j] = cost_ij
                costs[j][i] = cost_ij 

    return costs, probability

def compute_prizes_costs_deepslow(df, slewrate, dwell_time):
    ra = df['RA'].values
    dec = df['Dec'].values
    probability = df['Probability'].values

    n = len(ra)
    costs = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                costs[i][j] = dwell_time
            else:
                ra_diff = abs(ra[j] - ra[i])
                ra_diff = min(ra_diff, 360.0 - ra_diff)
                dec_diff = abs(dec[j] - dec[i])
                travel_time = max(ra_diff, dec_diff) / slewrate
                costs[i][j] = travel_time + dwell_time
                costs[j][i] = costs[i][j]

    return costs, probability


def incorporate_initpos(df, init_pos):
    row_to_copy = df[df['Rank'] == init_pos]

    if not row_to_copy.empty:
        new_row = row_to_copy.copy()
        new_row.loc[:, 'Rank'] = 0 
        new_row.loc[:, 'Probability'] = 0.0
        new_df = pd.concat([new_row, df], ignore_index=True)
    else:
        print(f"Rank {init_pos} not found in the DataFrame.")

    print(df.head(5))
    return new_df



def build_graph(df, slew_rate, dwell_time, is_deepslow):
    if is_deepslow:
        return compute_prizes_costs_deepslow(df, slew_rate, dwell_time)
    else:
        return compute_prizes_costs(df, slew_rate, dwell_time)


def compute_path_cost(costs, path):
    cost = 0.0
    for i in range(1, len(path), 1):
        cost += costs[path[i-1]][path[i]]
    return cost

def compute_cdf(path, costs, prizes, deadlines):
    cdf = []
    cost = 0.0
    prize = prizes[path[0]]
    i = 0  
    
    for j in range(1, len(path)):
        cost += costs[path[j-1]][path[j]]
        
        if i < len(deadlines) and cost > deadlines[i]:
            cdf.append(prize)
            i += 1
        prize += prizes[path[j]]

    while i < len(deadlines):
        cdf.append(prize)
        i += 1

    return cdf

data_file='../../data_RTSS/filtered_GW191219_163120.fits_7dt.csv'
# data_file='../../data_RTSS/filtered_GW191219_163120.fits_slow_deep.csv'
# data_file = "../../data_RTSS/filtered_GW200112_155838.fits_7dt.csv"
# data_file = "../../data_RTSS/filtered_GW200112_155838.fits_slow_deep.csv"
data_df = pd.read_csv(data_file)
data_df = incorporate_initpos(data_df, 100)
costs, prizes = build_graph(data_df, slew_rate=50, dwell_time=1, is_deepslow=False)
print(len(costs))


# result_file='result_50_200_500.csv'
result_file='result_100_300_700.csv'
# result_file='result_deepslow_100_300_700.csv'
# result_file='result_deepslow_200_600_1400.csv'
# result_file='result2_200_400_800.csv'
# result_file='result2_200_500_1200.csv'
# result_file='result2_deepslow_200_500_1200.csv'
# result_file='result2_deepslow_500_1500_3000.csv'
result_df = pd.read_csv(result_file)

paths = {}  
methods = ["d1_d2_d3", "d2_d3", "d3"]

for method in methods:
    method_df = result_df[result_df["Method"] == method]
    if len(method_df) > 0:
        paths[method] = list(map(int, method_df.iloc[0]["Path"].split()))
    else:
        print(f"No row found for method {method} with these parameters.")

print("\nAvailable Paths:", list(paths.keys()))


cdfs = {}
deadlines = [0, 100, 200, 300, 400, 500, 600, 700]
# deadlines = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400]

for method in methods:
    if method in paths: 
        # print(paths[method])
        cost = compute_path_cost(costs,paths[method])
        print(f"{method} cost: {cost}")
        cdfs[method] = compute_cdf(paths[method], costs, prizes, deadlines)

#Plot CDFs
plt.figure(figsize=(8, 6))

for method, cdf in cdfs.items():
    plt.plot(deadlines, cdf, marker='o', label=method)

plt.xlabel("Deadline (TimeSec)")
plt.ylabel("Cumulative Probability (SumProb)")
plt.title("CDF vs Deadline for Different Methods")
plt.grid(True)
plt.legend()
plt.show()


# print(f"costs[3][4]: {costs[3][4]}")
# print(f"costs[17][18]: {costs[18][19]}")