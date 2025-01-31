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

## Initial Test Result:

### 1. **7DT (Combined)**  
| Algorithm          | Dwell Time (S) | Budget (S) | Sum Probability | Total Cost | Num Tiles in Path | Computation Time (S) |
|------------------- |----------------|------------|-----------------|------------|-------------------|----------------------|
| Greedy Strategy    | 1              | 30         | 0.852429        | 29.9094    | 14                | 2.9e-05              |
| Ant Colony         | 1              | 30         | 0.852429       | 29.9094    | 14                | 4.7569               |
| Genetic Algorithm  | 1              | 30         | 0.939657        | 29.9756    | 22                | 0.2823               |
| Greedy Strategy    | 10             | 300        | 0.970533       | 298.332    | 27                | 7.1e-06              |
| Ant Colony         | 10             | 300        | 0.970533        | 298.714    | 27                | 10.306               |
| Genetic Algorithm  | 10             | 300        | 0.97527        | 297.902    | 28                | 0.3762               |
| Greedy Strategy    | 60             | 1800       | 0.978318       | 1769.1    | 29                | 6.8e-06              |
| Ant Colony         | 60             | 1800       | 0.978318        | 1798.67   | 29                | 11.4599              |
| Genetic Algorithm  | 60             | 1800       | 0.978318        | 1786.08   | 29                | 0.3874               |
| Greedy Strategy    | 600            | 18000      | 0.978318        | 17422.8    | 30                | 6.2e-06              |
| Ant Colony         | 600            | 18000      | 0.978318        | 17455.5    | 30                | 11.23.6              |
| Genetic Algorithm  | 600            | 18000      | 0.978318        | 17478.5    | 30                | 0.3785               |

---

### 2. **7DT (Separate)**  
| Algorithm          | Dwell Time (S) | Budget (S) | Sum Probability | Total Cost | Num Tiles in Path | Computation Time (S) |
|------------------- |----------------|------------|-----------------|------------|-------------------|----------------------|
| Greedy Strategy    | 1              | 50         | 0.99691        | 49.0781    | 22                | 5.5e-06              |
| Ant Colony         | 1              | 50         | 0.997138        | 48.7641    | 24                | 0.8194               |
| Genetic Algorithm  | 1              | 50         | 0.99985        | 49.9598    | 39                | 0.4688                |
| Greedy Strategy    | 10             | 200        | 0.994884        | 198.849    | 18                | 1.6e-06              |
| Ant Colony         | 10             | 200        | 0.994884        | 198.849    | 18                | 0.6935               |
| Genetic Algorithm  | 10             | 200        | 0.994884        | 196.755    | 18                | 0.2241               |
| Greedy Strategy    | 60             | 1000       | 0.993745        | 977.701    | 16                | 5.9e-06              |
| Ant Colony         | 60             | 1000       | 0.993745        | 979.998    | 16                | 0.6411              |
| Genetic Algorithm  | 60             | 1000       | 0.993745        | 996.068    | 16                | 0.1958               |
| Greedy Strategy    | 600            | 10000      | 0.993745        | 9617.7     | 16                | 6.7e-06              |
| Ant Colony         | 600            | 10000      | 0.993745        | 9617.7     | 16                |  0.6690               |
| Genetic Algorithm  | 600            | 10000      | 0.993745        | 9635.5     | 16                | 0.2005               |

---

### 3. **Deep Slow Telescope**  
| Algorithm          | Budget (S) | Sum Probability | Total Cost | Num Tiles in Path | Computation Time (S) |
|------------------- |------------|-----------------|------------|---------------|----------------------|
| Greedy Strategy    | 500        | 0.272592        | 496.815    | 10                | 3.7e-05              |
| Ant Colony         | 500        | 0.709269        | 498.767    | 29                | 175.177              |
| Genetic Algorithm  | 500        | 0.849264        | 499.51     | 68                | 1.19056              |
| Greedy Strategy    | 1000       | 0.53719         | 999.779    | 23                | 1.6e-06              |
| Ant Colony         | 1000       | 0.912846        | 999.967    | 55                | 377.565               |
| Genetic Algorithm  | 1000       | 0.956065        | 999.199     | 116                | 3.4668               |
