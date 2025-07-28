# Artifact for RTSS 2025
This repository contains the source code, datasets, and analysis tools supporting the paper *"Probabilistic Response-Time-Aware Search for Transient Astrophysical Phenomena"*.


## Overview
This artifact addresses the problem of scheduling astronomical events follow-up observations under real-time constraints. The task is to select and order sky tiles for telescope observations to maximize the probability of detecting an afterglow, subject to a hard deadline and slew-time constraints.

### Included Components
- C++ implementations of all algorithms: GCP, ILP, Genetic, Greedy
- Precomputed results in CSV format for:
  - Small and large instances under a single deadline
  - Small instances incorporating worst-case execution time (WCET)
  - Large instances under multi-deadline settings
- Jupyter notebooks to reproduce plots and figures from the paper
- Long-running Python scripts to recompute results from scratch


## Directory Structure

```
.
├── include/                         # C++ header files
├── src/                             # C++ source implementations
├── CMakeLists.txt                   # CMake build configuration

├── Data/
│   ├── small/          # Tiling CSV files for small-scale sky maps 
│   ├── large/          # Tiling CSV files for large-scale sky maps 
│   ├── large_wcet/     # Subset of tiling CSV files for large-scale maps used in WCET-aware experiments
│   └── wcet/           # Worst-case execution time measurements

├── results/
│   ├── precomputed_results/        # Precomputed outputs used in the paper
│   │   ├── small/                  # FoM results (non-ILP) for small instances
│   │   ├── small_gurobi/          # ILP-based FoM results for small instances
│   │   ├── large/                 # FoM results for large instances
│   │   ├── small_with_wcet/       # WCET-aware FoM for small instances
│   │   ├── multi_deadline/        # Multi-deadline scenario results
│   │   └── analysis_precomputed_results.ipynb  # Notebook for Figures 7–10 in paper
│
│   ├── recomputed_results/        # Outputs from re-running experiments
│   │   ├── small/                 
│   │   ├── small_with_wcrt/       
│   │   ├── large/                 
│   │   ├── multi_deadline/        
│   │   └── analysis_recomputed_results.ipynb  # Visualization of recomputed results
│
│   ├── run_small_instances.py         # Batch run script for small instances
│   ├── run_large_instances.py         # Batch run script for large instances
│   ├── run_instances_wcet.py          # Batch run script (WCET-aware)
│   └── run_multi_deadline.py          # Batch run script for multi-deadline setup

├── requirements.txt                # Python dependencies
└── README.md


```

## Reproducing Paper Figures
### 1. Set Up the Python Environment

We recommend setting up a virtual environment for Python, This will install all required packages, including: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`

```bash
python3 -m venv rtss
source rtss/bin/activate
pip install -r requirements.txt
```
---

### 2. Reproducing Result in a Jupyter notebook

To reproduce all figures for experiment results:

```bash
cd results/recompute_results
jupyter notebook analysis_precomputed_results.ipynb
```
The notebook includes:
- Average percentage deviation from ILP baseline across deadlines for small instances(Fig.7)
- Percentage deviation in FoM from the GCP baseline and Computation time vs Deadline presented in log-scaled plot for large instances (Fig. 8)
- FoM percentage deviation from GCP after accounting for WCET. for five small instances (Fig 9)
- Expected FoM progression over the multi deadlines (Fig. 10)

All required `.csv` results are precomputed and stored in `results/precomputed_results`.

---


## Running Full Experiments

If you'd like to re-run the full evaluation (5+ hours), use the provided batch scripts:

---

### 1. Install Dependencies

#### 1.1 Gurobi Optimizer (Required)

- **Download**: [https://www.gurobi.com/downloads/](https://www.gurobi.com/downloads/)
- **Installation Guide**: [How to Install Gurobi](https://support.gurobi.com/hc/en-us/articles/4534161999889)
- **License Setup**: Gurobi requires a license to run.
  - [How to Retrieve and Set Up a Gurobi License](https://support.gurobi.com/hc/en-us/articles/12872879801105)

After installation, ensure the following environment variable is set in your shell:

```bash
export GUROBI_HOME=/path/to/gurobi
export PATH="${GUROBI_HOME}/bin:$PATH"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib:$LD_LIBRARY_PATH"
```
#### 1.2 LEMON Graph Library (Required)

- **Download**: [LEMON 1.3.1 Source](http://lemon.cs.elte.hu/pub/sources/lemon-doc-1.3.1.tar.gz)
- **Installation Guide**: [LEMON Installation (Linux)](http://lemon.cs.elte.hu/trac/lemon/wiki/InstallLinux)

After building LEMON, set the following environment variables in `CMakeLists.txt`:
```cmake
set(LEMON_SOURCE_DIR "/your/path/to/lemon-1.3.1")
set(LEMON_BUILD_DIR "/your/path/to/lemon-1.3.1/build")
```

Alternatively, you may set the following environment in your shell and comment the two line in `CMakeLists.txt`:
```bash
export LEMON_SOURCE_DIR=/path/to/lemon-1.3.1
export LEMON_BUILD_DIR=/path/to/lemon-1.3.1/build
```
---

### 2. Build the C++ Executables

Once all dependencies are installed, you can build the C++ project with:

```bash
mkdir build && cd build
cmake ..
make -j
```


### 3. Recompute Results from Scratch

Each Python script corresponds to a different experiment setting. 

#### 3.1 Small Instances (~ 3 hours)
```bash
python3 run_small_instances.py
```
Results will be saved to `results/recompute_results/small`.

#### 3.2 Large Instances (~ 30 minutes)

```bash
python3 run_large_instances.py
```
Results will be saved to `results/recompute_results/large`.

#### 3.3 Small Instances with WCET (~ 20 minutes)

```bash
python3 run_small_instances_wcet.py
```
Results will be saved to `results/recompute_results/instances_with_wcet`.

#### 3.4 Multi-Deadline Large Instances (~ 30 minutes)

```bash
python3 run_multi_deadline.py
```
Results will be saved to `results/recompute_results/multi_deadlines`.

---

#### 3.5. Visualizing the Results

```bash
cd results/recompute_results
jupyter notebook analysis_recomputed_results.ipynb
```
All required `.csv` results are stored in `results/recomputed_results`.

