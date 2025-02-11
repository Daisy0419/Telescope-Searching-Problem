
import math
import csv

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

def main():
    input_file = "../data/7dt_separate.csv"   
    output_file = "ampl_input_matrix.dat"     

    telescope_speed = 5.0
    dwell_time = 1.0

    costs, probability = read_data(input_file, telescope_speed, dwell_time)

    write_ampl_data_as_matrix(
        output_filename=output_file,
        costs=costs,
        probability=probability,
        start=1,
        budget=39
    )
    print(f"AMPL data (matrix form) written to {output_file}")

if __name__ == "__main__":
    main()

