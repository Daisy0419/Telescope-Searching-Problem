# Telescope Searching Problem

## Overview
Telescope Searching Problem focuses on optimizing telescope observation schedules to maximize the probability of detecting GRB within a given deadline. 

### Key Components
- **Probability**: The probability for observing GRB in each tile, total probability sum up to 1.
- **Slew Time**: The time required to move between tiles, calculated from the angular distance and telescope slewing speed.
- **Dwell Time**: The observation time for a tile, determining the probability of detecting the transient. Linear probability growth until certainty (1) after a complete dwell time is also considered.

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
1. **7DT (Combined)**  
| Algorithm | Dwell Time (S) | Budget (S) | Sum Probability | Total Cost | Num Tiles in Path | Computation Time (S) | Path |
|-----------|----------------|-----------------|------------|-------------------|----------------------|------|---|
| Greedy Strategy   | 1 | 30 | 0.884264 | 29.8567 | 16 | 3.0e-06 |
| Ant Colony        | 1 | 30 | 0.923233 | 29.1531 | 19 | 6.7289 |
| Genetic Algorithm | 1 | 30 | 0.931198 | 29.9756   | 21 | 0.2734 |
| Greedy Strategy   |10 | 300 | 0.97527 | 292.229 | 28 | 1.5e-05 |
| Ant Colony        |10 | 300 | 0.978318 | 299.929 | 29 | 10.829 |
| Genetic Algorithm |10 | 300 | 0.978318 | 299.368 | 29 | 0.3762 |
| Greedy Strategy   |60 | 1800 | 0.981164 | 1762.79 | 30 | 3.5e-05 |
| Ant Colony        |60 | 1800 | 0.981164 | 1768.08 | 30 | 11.6985 |
| Genetic Algorithm |60 | 1800 | 0.981164 | 1778.42 | 30 | 0.4031 |
| Greedy Strategy   |600| 18000 | 0.981164 | 17422.8 | 30 | 1.5e-05 |
| Ant Colony        |600| 18000 | 0.981164 | 17429.8 | 30 | 11.6959 |
| Genetic Algorithm |600| 18000 | 0.981164 | 17440.8 | 30 | 0.4231 |

2. **7DT (Separate)**  
| Algorithm | Dwell Time (S) | Budget (S) | Sum Probability | Total Cost | Num Tiles in Path | Computation Time (S) | Path |
|-----------|----------------|-----------------|------------|-------------------|----------------------|------|---|
| Greedy Strategy   | 1 | 20 | 0.975682 | 19.0878 | 10 | 1.7e-06 |
| Ant Colony        | 1 | 20 | 0.984177 | 19.6886 | 11 | 0.3917 |
| Genetic Algorithm | 1 | 20 | 0.990751 | 19.7647   | 14 | 0.178 |
| Greedy Strategy   |10 | 200 | 0.995453 | 198.602 | 19 | 1.6e-06 |
| Ant Colony        |10 | 200 | 0.995453 | 197.284 | 19 | 0.7505 |
| Genetic Algorithm |10 | 200 | 0.995453 | 198.059 | 19 | 0.2305 |
| Greedy Strategy   |60 | 1200 | 0.995981 | 1159.55 | 20 | 1.8e-06 |
| Ant Colony        |60 | 1200 | 0.995981 | 1157.27 | 20 | 11.6985 |
| Genetic Algorithm |60 | 1200 | 0.995981 | 1174.21 | 20 | 0.8262 |
| Greedy Strategy   |600| 12000 | 0.995981 | 11419.5 | 20 | 1.9e-06 |
| Ant Colony        |600| 12000 | 0.995981 | 11417.5 | 20 | 0.7859 |
| Genetic Algorithm |600| 12000 | 0.995981 | 11434.4 | 20 | 0.2393 |

3. **Deep Slow Telescope**  
| Algorithm | Budget (S) | Sum Probability | Total Cost | Num Tiles in Path | Computation Time (S) | Path |
|-----------|----------------|-----------------|------------|-------------------|----------------------|------|---|
| Greedy Strategy   | 1000 | 0.966239 | 995.53 | 8 | 1.9e-06 |
| Ant Colony        | 1000 | 0.987839 | 995.53 | 11 | 0.3120 |
| Genetic Algorithm | 1000 | 0.951076 | 999.051| 8 | 0.1580 |
| Greedy Strategy   |2000 | 0.996134 | 1977.31 | 20 | 1.6e-06 |
| Ant Colony        |2000 | 0.995453 | 1963.19 | 19 | 0.7505 |
| Genetic Algorithm |2000 | 0.999826 | 1974.91 | 20 | 0.2305 |
