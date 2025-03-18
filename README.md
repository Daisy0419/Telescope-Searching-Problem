# Telescope Searching Problem

## Overview
Telescope Searching Problem focuses on optimizing telescope observation schedules to maximize the probability of detecting GRB within a given deadline. 

### Key Components
- **Probability**: The probability for observing GRB in each tile, total probability sum up to 1.
- **Slew Time**: The time required to move between tiles, calculated from the angular distance and telescope slewing speed.
- **Dwell Time**: The observation time for a tile, min time for a meaningful observation. 


### Key Components in modeled TSP problem
The problem is modeled as an extension of the prize-collecting Traveling Salesman Problem (TSP) with a budget constraint, where:
- **Budget**: The deadline, consumed by both slewing time (between tiles) and dwell time (per tile).
- **Prize**: The probability of detecting the GRB, derived from the probability distribution.
- **Cost**:  slewing time (moving between tiles) and dwell time (per tile).


## Input Data
1. **7DT (Combined)**  
   - **Tile Size**: 2.5 x 2.5 degrees  
   - **Probability Distribution**: BAYESTAR Sample  
   - **Slew Rate**: 50 degrees per second  
   - **Coordinates**: Right Ascension (RA), Declination (Dec)

2. **7DT (Separate)**  
   - **Tile Size**: 10 x 10 degrees  
   - **Probability Distribution**: BAYESTAR Sample  
   - **Slew Rate**: 50 degrees per second  

3. **Deep Slow Telescope**  
   - **Tile Size**: 0.5 x 0.5 degrees  
   - **Probability Distribution**: BAYESTAR Sample  
   - **Slew Rate**: 1 degree per second (both azimuth and elevation)  
   - **Slew time calculation**: `max(|RA_1 - RA_0|, |Dec_1 - Dec_0|)`


## Cost Computation
1. **Slew Time**  
`θ = arccos(sin(Dec₁)⋅sin(Dec₂) + cos(Dec₁)⋅cos(Dec₂)⋅cos(RA₁ − RA₂))`
`Slew Time = θ / Slew Rate`

2. **Dwell Time**
    Fix to 1, 10, 60, 600 seconds

## Experiment Results:

## Dataset: `7dt_combined.csv` | Budget: `15`, SlewRate: `50`, DwellTime: `1`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.85243 | 2.9312 | 14.6871 | `1, 2, 3, 8, 9, 4, 6, 11, 7, 5, 14, 13, 12, 10` |
| Genetic | 0.85243 | 3.45362 | 14.8548 | `1, 2, 4, 10, 8, 5, 6, 12, 13, 11, 7, 3, 14, 9` |
| Greedy | 0.85243 | 1e-05 | 14.4207 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14` |

## Dataset: `7dt_combined.csv` | Budget: `150`, SlewRate: `50`, DwellTime: `10`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.85243 | 1.73517 | 140.837 | `1, 5, 8, 10, 6, 4, 2, 3, 13, 7, 9, 12, 11, 14` |
| Genetic | 0.85243 | 0.51098 | 140.938 | `1, 8, 5, 10, 4, 2, 6, 3, 13, 11, 12, 7, 14, 9` |
| Greedy | 0.85243 | 1e-05 | 140.421 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14` |

## Dataset: `7dt_combined.csv` | Budget: `500`, SlewRate: `50`, DwellTime: `60`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.71401 | 0.87048 | 480.379 | `1, 8, 6, 2, 7, 3, 5, 4` |
| Genetic | 0.71401 | 0.5108 | 480.347 | `1, 8, 6, 7, 2, 4, 3, 5` |
| Greedy | 0.71401 | 2e-05 | 480.233 | `1, 2, 3, 4, 5, 6, 7, 8` |

## Dataset: `7dt_combined.csv` | Budget: `1550`, SlewRate: `50`, DwellTime: `100`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.86997 | 4.59908 | 1500.95 | `1, 2, 4, 6, 8, 5, 7, 15, 13, 11, 3, 12, 9, 14, 10` |
| Genetic | 0.86997 | 0.63755 | 1500.93 | `1, 2, 4, 8, 6, 15, 10, 9, 14, 12, 13, 11, 7, 3, 5` |
| Greedy | 0.86997 | 1e-05 | 1500.45 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15` |

## Dataset: `7dt_combined.csv` | Budget: `3000`, SlewRate: `50`, DwellTime: `600`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.43898 | 0.35106 | 2400.09 | `1, 4, 3, 2` |
| Genetic | 0.43898 | 2.74338 | 2400.12 | `1, 3, 4, 2` |
| Greedy | 0.43898 | 1e-05 | 2400.09 | `1, 2, 3, 4` |

## Dataset: `7dt_separate.csv` | Budget: `15`, SlewRate: `50`, DwellTime: `1`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.992 | 0.22476 | 14.9423 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 12` |
| Genetic | 0.992 | 1.40707 | 14.946 | `1, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4, 5, 13, 14` |
| Greedy | 0.992 | 1e-05 | 14.8248 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14` |

## Dataset: `7dt_separate.csv` | Budget: `150`, SlewRate: `50`, DwellTime: `10`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.992 | 0.2952 | 141.611 | `1, 2, 4, 3, 7, 5, 6, 9, 13, 8, 10, 11, 14, 12` |
| Genetic | 0.992 | 5.95519 | 141.454 | `1, 7, 6, 9, 12, 2, 11, 10, 14, 4, 5, 13, 3, 8` |
| Greedy | 0.992 | 1e-05 | 140.825 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14` |

## Dataset: `7dt_separate.csv` | Budget: `500`, SlewRate: `50`, DwellTime: `60`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.9515 | 1.19485 | 480.71 | `1, 6, 4, 2, 3, 5, 8, 7` |
| Genetic | 0.9515 | 0.89566 | 480.594 | `1, 7, 6, 2, 8, 3, 4, 5` |
| Greedy | 0.9515 | 1e-05 | 480.449 | `1, 2, 3, 4, 5, 6, 7, 8` |

## Dataset: `7dt_separate.csv` | Budget: `1550`, SlewRate: `50`, DwellTime: `100`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.99287 | 0.27058 | 1501.95 | `1, 6, 3, 2, 13, 4, 5, 11, 7, 10, 14, 9, 8, 15, 12` |
| Genetic | 0.99287 | 0.9639 | 1502.16 | `1, 10, 15, 12, 6, 11, 8, 4, 5, 9, 2, 14, 7, 3, 13` |
| Greedy | 0.99287 | 1e-05 | 1500.89 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15` |

## Dataset: `7dt_separate.csv` | Budget: `3000`, SlewRate: `50`, DwellTime: `600`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.79507 | 0.07111 | 2400.23 | `1, 4, 2, 3` |
| Genetic | 0.79507 | 0.37634 | 2400.23 | `1, 2, 4, 3` |
| Greedy | 0.79507 | 1e-05 | 2400.18 | `1, 2, 3, 4` |

## Dataset: `deep_slow.csv` | Budget: `500`, SlewRate: `1`, DwellTime: `10`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.54477 | 36.2081 | 494.178 | `1, 8, 21, 2, 22, 9, 23, 10, 24, 3, 25, 46, 11, 26, 12, 27, 4, 13, 14, 5, 15, 32, 16, 6` |
| Genetic | 0.52881 | 5.18154 | 490.566 | `1, 8, 21, 2, 9, 23, 10, 24, 3, 11, 26, 12, 4, 28, 13, 14, 30, 5, 15, 16, 6, 17` |
| Greedy | 0.25322 | 5e-05 | 489.686 | `1, 2, 3, 4, 5, 6, 7, 15` |

## Dataset: `deep_slow.csv` | Budget: `1000`, SlewRate: `1`, DwellTime: `10`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.78654 | 94.761 | 989.584 | `1, 8, 2, 21, 41, 22, 42, 9, 23, 10, 24, 3, 25, 11, 26, 12, 48, 27, 4, 28, 13, 14, 30, 5, 15, 31, 54, 16, 33, 6, 17, 34, 58, 18, 36, 7, 37, 62, 19, 20, 39, 64, 38, 63, 90, 60, 35, 59, 125` |
| Genetic | 0.81092 | 6.99113 | 992.537 | `1, 40, 8, 21, 2, 22, 42, 9, 44, 23, 10, 24, 3, 25, 46, 11, 12, 27, 4, 28, 13, 14, 30, 5, 31, 54, 15, 55, 32, 56, 16, 33, 57, 6, 34, 58, 17, 35, 59, 89, 90, 60, 18, 36, 61, 62, 63, 38, 64, 39, 20, 19, 37, 7, 91` |
| Greedy | 0.44838 | 5e-05 | 999.448 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16` |

## Dataset: `deep_slow.csv` | Budget: `2000`, SlewRate: `1`, DwellTime: `10`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.95334 | 306.899 | 1989.47 | `1, 8, 21, 40, 65, 2, 9, 22, 42, 41, 66, 24, 3, 11, 25, 46, 12, 27, 4, 13, 28, 50, 29, 51, 14, 5, 31, 54, 32, 55, 15, 16, 33, 6, 17, 34, 58, 7, 19, 20, 39, 64, 38, 63, 37, 62, 18, 36, 61, 91, 35, 59, 88, 57, 86, 56, 85, 84, 83, 82, 117, 116, 81, 30, 53, 52, 80, 79, 113, 78, 112, 77, 49, 76, 48, 75, 26, 47, 74, 107, 72, 73, 106, 45, 44, 10, 23, 43, 69, 101, 71, 103, 70, 102, 141, 67, 68, 100, 139, 99` |
| Genetic | 0.97174 | 4.20423 | 1998.8 | `1, 8, 65, 40, 21, 2, 41, 66, 67, 22, 42, 9, 43, 23, 69, 101, 68, 100, 99, 98, 97, 96, 102, 10, 71, 44, 70, 45, 104, 103, 24, 3, 72, 46, 25, 11, 26, 47, 105, 106, 73, 74, 75, 48, 27, 76, 77, 78, 79, 80, 81, 31, 54, 83, 84, 85, 87, 88, 89, 90, 60, 91, 61, 92, 62, 94, 95, 64, 39, 20, 38, 63, 19, 93, 37, 7, 36, 18, 35, 59, 17, 58, 34, 6, 57, 86, 33, 16, 56, 32, 55, 15, 82, 5, 53, 30, 14, 52, 29, 51, 13, 28, 50, 4, 49, 12, 107, 108, 109, 110, 111, 112` |
| Greedy | 0.75065 | 6e-05 | 1999.93 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38` |

## Dataset: `deep_slow.csv` | Budget: `1000`, SlewRate: `1`, DwellTime: `60`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.36167 | 21.058 | 993.264 | `1, 8, 21, 2, 9, 10, 3, 11, 4, 13, 14, 5, 15` |
| Genetic | 0.36781 | 0.89742 | 993.264 | `1, 8, 2, 9, 10, 3, 11, 12, 4, 13, 14, 5, 15` |
| Greedy | 0.25475 | 6e-05 | 975.638 | `1, 2, 3, 4, 5, 6, 7, 12` |

## Dataset: `deep_slow.csv` | Budget: `2000`, SlewRate: `1`, DwellTime: `60`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.58631 | 49.3124 | 1995.68 | `1, 2, 8, 21, 22, 42, 9, 10, 3, 11, 12, 27, 26, 4, 13, 14, 5, 15, 31, 6, 17, 7, 18, 36, 19, 37` |
| Genetic | 0.62448 | 1.62929 | 1998.26 | `1, 21, 8, 2, 22, 9, 23, 10, 24, 3, 25, 11, 12, 4, 13, 14, 5, 15, 32, 16, 6, 17, 18, 7, 19, 20, 38` |
| Greedy | 0.48835 | 6e-05 | 1976.75 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18` |

## Dataset: `deep_slow.csv` | Budget: `3000`, SlewRate: `1`, DwellTime: `60`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.74266 | 89.6049 | 2996.45 | `1, 2, 8, 21, 9, 3, 11, 4, 13, 5, 6, 17, 7, 19, 20, 39, 38, 37, 18, 36, 35, 34, 16, 33, 32, 15, 31, 14, 30, 29, 28, 12, 27, 26, 25, 10, 24, 23` |
| Genetic | 0.74036 | 3.08112 | 2999.51 | `1, 21, 2, 8, 9, 10, 3, 11, 4, 12, 13, 5, 14, 15, 6, 34, 35, 36, 7, 20, 19, 18, 17, 33, 16, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 40` |
| Greedy | 0.61977 | 5e-05 | 2950.7 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26` |

## Dataset: `deep_slow.csv` | Budget: `5000`, SlewRate: `1`, DwellTime: `60`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.89879 | 204.234 | 4978.33 | `1, 2, 8, 21, 9, 3, 11, 4, 13, 5, 14, 30, 29, 12, 27, 26, 10, 24, 23, 22, 42, 41, 40, 65, 43, 44, 45, 25, 46, 47, 48, 49, 28, 50, 51, 52, 31, 54, 55, 56, 57, 6, 17, 7, 19, 20, 39, 18, 36, 35, 37, 62, 63, 38, 61, 60, 59, 34, 58, 16, 33, 32, 15, 53, 81` |
| Genetic | 0.90873 | 6.73349 | 4999.41 | `1, 8, 2, 41, 21, 40, 42, 22, 23, 24, 9, 43, 44, 10, 3, 45, 50, 28, 29, 30, 53, 54, 55, 15, 31, 32, 56, 33, 16, 57, 6, 34, 35, 36, 18, 37, 7, 19, 63, 38, 39, 20, 64, 62, 61, 60, 59, 17, 58, 5, 14, 52, 51, 13, 4, 49, 48, 12, 27, 26, 11, 47, 25, 46, 69, 68, 66, 65` |
| Greedy | 0.82454 | 5e-05 | 4985.75 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48` |

## Dataset: `deep_slow.csv` | Budget: `10000`, SlewRate: `1`, DwellTime: `600`

| Method     | SumProb | TimeSec | Path_Cost | Path |
|-----------|---------|---------|-----------|------|
| AntColony | 0.42707 | 30.9023 | 9778.06 | `1, 2, 8, 10, 9, 3, 5, 15, 7, 6, 12, 14, 4, 11, 13` |
| Genetic | 0.44838 | 0.75042 | 9981.03 | `1, 8, 2, 9, 10, 3, 11, 12, 4, 13, 14, 5, 15, 6, 7, 16` |
| Greedy | 0.42707 | 5e-05 | 9810.8 | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15` |

