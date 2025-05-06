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
    TBD

## Experiment Results:

