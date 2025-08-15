# Artifact for RTSS 2025
This repository contains the source code, datasets, and analysis tools supporting the paper *"Probabilistic Response-Time-Aware Search for Transient Astrophysical Phenomena."*

## System Requirements

- OS: Linux (some instructions are Ubuntu-specific)
- CPU: Minimum dual-core, 8+ cores preferred
- RAM: 
- Required Storage: **Total: 4GB**
  - Repo: 150MB
  - Conda environments: 3.8GB
- Optional Storage:
  - Docker: 500MB
  - Docker Container: 3.8


## Overview
This artifact addresses the problem of scheduling follow-up observations of astrophysical transients under real-time constraints. The task is to select and order sky tiles for telescope observations to maximize the expected figure of merit (FoM) associated with successful detection of an afterglow within one or more hard deadlines and given dwell-time and slew-time constraints.

### Included Components
- C++ implementations of all algorithms: GCP, ILP, Genetic, Greedy
- Precomputed results in CSV format for:
  - Small and large problem instances under a single deadline
  - Small problem instances incorporating maximum-observed execution time (MOET)
  - Large problem instances under multi-deadline settings
- Jupyter notebooks to reproduce plots and figures from the paper
- Long-running Python scripts to recompute results from scratch


## Directory Structure

```
.
├── include/                         # C++ header files
├── src/                             # C++ source implementations
├── CMakeLists.txt                   # CMake build configuration

├── skytiling/          
│   └── precompute_tile_maps.sh     # Precomputes tile boundaries by projecting telescope FoV
│   └── flatten_healpix.py          # Convert multi-order hierarchical HEALPix to flat HEALPix format
│   └── produce_tilings.py          # Creates all tilings

├── data/               # Tiling CSV files for small-scale sky maps
│   ├── small/          # Tiling CSV files for small-scale sky maps 
│   ├── large/          # Tiling CSV files for large-scale sky maps 
│   ├── large_moet/     # Subset of large-scale maps used in MOET-aware experiments
│   └── moet/           # Maximum-observed execution time measurements

├── results/
│
│   ├── figures/                   # Figures 7–10 in paper
|
│   ├── precomputed_results/       # Precomputed outputs used in the paper
│   │   ├── small/                 # FoM results for small instances
│   │   ├── large/                 # FoM results for large instances
│   │   ├── small_with_moet/       # MOET-aware FoM for small instances
│   │   └── multi_deadline/        # Multi-deadline scenario results
│
│   ├── recomputed_results/        # Outputs from re-running experiments
│   │   ├── small/      
│   │   ├── large/                            
│   │   ├── small_with_moet/       
│   │   └── multi_deadline/        
│
│   ├── visualize_results.ipynb        # Notebook for Figures 7–10 in paper
│   ├── run_small_instances.py         # Batch run script for small instances
│   ├── run_large_instances.py         # Batch run script for large instances
│   ├── run_small_instances_moet.py    # Batch run script (MOET-aware)
│   └── run_multi_deadline.py          # Batch run script for multi-deadline setup

├── rtss25-sky-tiling.yml               # Python dependencies for tiling code
├── rtss25-telescope-search.yml         # Python dependencies
└── README.md


```
## 1 Environment Setup
You can run the artifact via **Docker (recommended)** or a **Local Setup**. A Gurobi license is needed only to run ILP-based experiments.
### 1.1 (Preliminary, Optional) Obtaining a Gurobi License

Our algorithms include an ILP implementation of the orienteering problem, which is solved using the commercially-available Gurobi Optimizer.

Gurobi requires a license (`gurobi.lic`) to run.
- Refer to: [How to Retrieve and Set Up a Gurobi License](https://support.gurobi.com/hc/en-us/articles/12872879801105)
- **Note:** If you do not have a Gurobi license, you can still run experiments 3.2, 3.3, and 3.4 below, which do not invoke the ILP solver.

If you are an academic user, Gurobi provides **free academic licenses**:
  - [Free Academic License](https://www.gurobi.com/academia/academic-program-and-licenses/)
  - **Note**: If you're using an academic license and intend to run the experiments using our provided Docker container, be sure to request an Academic WLS License (floating license). Named-User Academic Licenses are not compatible with Docker containers.


To obtain an academic license:

1. Navigate to the Gurobi portal. https://portal.gurobi.com/
2. Login or register to create a free Gurobi account using your academic email address.
3. Navigate to Gurobi's academic license request page. https://portal.gurobi.com/iam/licenses/request/?type=academic
4. Under "WLS Academic," click, "Generate Now!"
5. Download the generated `gurobi.lic` file.
6. Move or copy the file to the path of your choice. All commands listed hereafter assume it is in `~/gurobi.lic`.


### 1.2 (Option A, Preferred) Using the Provided Docker Container

#### 1.2.1 Install Docker

You may install Docker according to [these instructions](https://docs.docker.com/engine/install/). Here, we include the instructions for Ubuntu distributions:

1. Set up Docker's `apt` repository:

```bash
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] \
  https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
```

2. Install the latest Docker packages.

```bash
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

#### 1.2.2 Pull the Docker Image
```bash
sudo docker pull ghcr.io/daisy0419/rtss25-op-solver:1.0
```
All dependencies are pre-installed, and project binaries are precompiled in the image. You can jump to Reproducing Paper Figures or Running Full Experiments.

### 1.3 (Option B) Local Installation

#### 1.3.1 Clone Repository

Clone this repository to the path of your choice. All commands listed hereafter assume it is placed directly into your home directory.

```bash
cd ~
git clone https://github.com/Daisy0419/Telescope-Searching-Problem/releases/tag/rtss2025_artifact
```

#### 1.3.2 Python Environment Setup

We recommend setting up a [conda](https://docs.conda.io/en/latest/) environment for Python.

If you do not have conda installed locally:

```bash
cd ~/Telescope-Search-Problem
```
Download and install Miniconda (change ~/conda to your preferred location)
```bash
export CONDA_DIR=~/conda
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p "${CONDA_DIR}"
rm miniconda.sh
```

Make conda available in your shell
```bash
. "${CONDA_DIR}/etc/profile.d/conda.sh"
conda config --system --set channel_priority flexible
```
Once conda is available, create the python environments using provided yaml files:
```bash
conda env create -f rtss25-sky-tiling.yml
conda env create -f rtss25-telescope-search.yml
```

Activate **rtss25-sky-tiling** before running tiling scripts or **rtss25-telescope-search** before running all other scripts:
```bash
conda activate rtss25-sky-tiling
# or
conda activate rtss25-telescope-search
```

---
#### 1.3.3 C++ Environment Setup

**(1) Gurobi Optimizer (Required)**

The Gurobi Optimizer is used to solve our ILP approach to the orienteering problem.

1. Download and extract Gurobi to the directory of your choice. All commands listed hereafter assume it is placed directly in your home directory.

```bash
cd ~
wget https://packages.gurobi.com/12.0/gurobi12.0.3_linux64.tar.gz
tar xvfz gurobi12.0.3_linux64.tar.gz
```

2. Set the necessary environment variables in your shell (change `~/gurobi1203` to your preferred location).

```bash
export GUROBI_HOME=~/gurobi1203/linux64
export PATH="${GUROBI_HOME}/bin:$PATH"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib:$LD_LIBRARY_PATH"
```

**Note**: If you don't have Gurobi license and you do not plan to run experiments involving Gurobi, you still need to install Gurobi in order to compile the project code (due to build-time linking requirements).


#### (2) LEMON Graph Library (Required)

The Lemon Graph Library is used to compute the minimum-weight perfect matching used in our Greedy Christofides Pathfinding algorithm. 

1. Download, extract, and build the LEMON Graph Library to the directory of your choice. All commands listed hereafter assume it is placed directly in your home directory.

```bash
wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
tar xvfz lemon-1.3.1.tar.gz
cd lemon-1.3.1
mkdir build && cd build
cmake ..
make -j
sudo make install
```

2. Set the necessary environment variables in your shell (change `~/lemon-1.3.1` to your preferred location).

```bash
export LEMON_SOURCE_DIR=~/lemon-1.3.1
export LEMON_BUILD_DIR=~/lemon-1.3.1/build
```

**(3) Build the C++ Executables**

Once all dependencies are installed, you can build the C++ project with:

```bash
cd ~/Telescope-Search-Problem
mkdir build && cd build
cmake ..
make -j
```

This produces two binaries in build/. The only difference between them is the main program:
- **ts** — built from src/main.cpp. Entry point for batch experiments used in the paper, invoked by the Python scripts in results/.
- **op** — built from src/main_custom.cpp. Entry point for single-case runs, one algorithm on one (skymap, budget, slewrate) instance.

--- 

## 2 Reproducing Paper Figures
you can visualize the result via Jupyter notebook either via **container** or **locally**.
### 2.1 Run Jupyter notebook via Docker Container
```bash
docker run --rm -it -p 8888:8888 \
  -v "$PWD:/workspace" \
  ghcr.io/daisy0419/rtss25-op-solver:1.0 \
  bash -lc 'conda run -n rtss25-telescope-search \
    jupyter lab --ip=0.0.0.0 --port=8888 --no-browser \
      --IdentityProvider.token="" \
      --ServerApp.root_dir=/workspace \
      --allow-root'
```

and then open http://localhost:8888 in a browser and natigate to **results/visualize_results.ipynb** in sidebar.

### 2.2 Run Jupyter notebook Locally

```bash
conda activate rtss25-telescope-search
cd ~/Telescope-Searching-Problem/results
jupyter notebook visualize_results.ipynb
```

### 2.3 Reproducing Result in a Jupyter notebook
The notebook includes:
- Average percentage deviation from ILP baseline across deadlines for small instances(Fig.7)
- Percentage deviation in FoM from the GCP baseline and Computation time vs Deadline presented in log-scaled plot for large instances (Fig. 8)
- FoM percentage deviation from GCP after accounting for MOET. for five small instances (Fig 9)
- Expected FoM progression over the multi deadlines (Fig. 10)

All required `.csv` results are precomputed and stored in `results/precomputed_results`.

The notebook saves figures to `results/figures`.

---

## 3 Running Full Experiments

To re-run the full set of experiments (~5+ hours total runtime), you may either set up the environment locally or use the provided Docker container. Experiments can then be executed using the provided batch scripts.

### 3.1 Set Up the Run Environment
#### 3.1.1 Run Experiment in a Docker Container
If you have a Gurobi license on your local machine, mount the license into the container.

**License Note**: If you're using an academic license, be sure to request an Academic WLS License (floating license). Named-User Academic Licenses are not compatible with Docker containers.

```bash
sudo docker run --rm -it -v "~/gurobi.lic:/licenses/gurobi.lic:ro" -e GRB_LICENSE_FILE=/licenses/gurobi.lic ghcr.io/daisy0419/rtss25-op-solver:1.0
```

If you do not have a Gurobi license, you can still run experiments that do not rely on ILP-based solvers:
```bash
sudo docker run --rm -it ghcr.io/daisy0419/rtss25-op-solver:1.0
cd results
```

#### 3.1.2 Run Experiment Locally
```bash
cd ~/Telescope-Searching-Problem
```

### 3.2 Regenerate Sky Tilings (Data Preparation, Optional)

The sky tiles used for the small and large problem instances evaluated are already stored in the `data/small` and `data/large` directories.
Optionally, you may regenerate them.
All commands in this step are run from the `sky_tiling` directory.

```bash
cd sky_tiling
```

#### 3.2.1 Precompute FoV Projections (~ 10 minutes)

Each telescope's FoV is projected onto the unit sphere to precompute tile boundaries.

```bash
conda activate rtss25-sky-tiling
bash precompute_tile_maps.sh
```

This step is optional. 
Files already exist in `sky_tiling/tile_center_files` and `sky_tiling/tile_pixel_maps`.

Recomputed files appear in those same directories, with file names containing `recomputed`.

#### 3.2.2 Flatten HealPix (~ 10 seconds)

Multi-order hierarchical HEALPix maps in `data/ligo_healpix`
need to be flattened into a non-hierarchical format. 

```bash
conda activate rtss25-telescope-search
python flatten_healpix.py
```

This step is optional.
Files already exist in `data/ligo_healpix_flattened`;
they will be overwritten by this script.

#### 3.2.3 Generate Tiles (~ 15 seconds)

Intersection of flattened HEALPix maps and precomputed tilings are intersected
to produce the tiles, with probabilities, that serve as the inputs to the search problem.

```bash
conda activate rtss25-sky-tiling
python produce_tilings.py
```

Files in `data/small` and `data/large` will be overwritten.


### 3.3 Rerun all Experiments
Each Python script corresponds to a different experiment setting. 
If you do not have a Gurobi license, you can still run experiments 3.3.2, 3.3.3, and 3.3.4, which do not rely on ILP solvers.

#### 3.3.1 Small Instances (~ 3 hours)
```bash
python3 run_small_instances.py
```
Results will be saved to `results/recompute_results/small`.

#### 3.3.2 Large Instances (~ 30 minutes)

```bash
python3 run_large_instances.py
```
Results will be saved to `results/recompute_results/large`.

#### 3.3.3 Small Instances with MOET (~ 20 minutes)

```bash
python3 run_small_instances_moet.py
```
Results will be saved to `results/recompute_results/instances_with_moet`.

This script evaluates algorithms with budgets that account for MOET (Maximum Observed Execution Time). It runs in **three parts**:

(1) Collect raw runs to estimate MOET

For each skymap and each budget in {10,20,…,100}, the script runs Greedy, Genetic, and GCP 5 times each using the same original budget (no MOET applied yet). Results are stored in instances_with_moet/get_moet/.

(2) Compute MOET per (Skymap, Method, Budget) 

Across the 5 runs in (1), it then computes:

```bash
MOET(Skymap, Method, Budget) = max(TimeSec)
```

and writes save results to results/recomputed_results/instances_with_moet/moet/<skymap>.csv

(3) Re-run with MOET-adjusted budgets 

For each (Skymap, Method, Budget), the budget is reduced by the corresponding MOET:
```bash
adjusted_budget = max(0, original_budget - MOET(Skymap, Method, original_budget))
```
It then runs the instances again with these adjusted budgets and save result to

results/recomputed_results/instances_with_moet/result_with_moet/out_<skymap>.csv


#### 3.3.4 Multi-Deadline Large Instances (~ 30 minutes)

```bash
python3 run_multi_deadline.py
```
Results will be saved to `results/recompute_results/multi_deadlines`.

---

### 3.4. Visualizing the Results
Again, you can visualize the result via Jupyter notebook either **inside the container** or **locally**.
- **Inside the container**
```bash
conda activate rtss25-telescope-search
jupyter notebook --ip=0.0.0.0 --no-browser
```
Open http://localhost:8888 in a browser and navigate to **results/visualize_results.ipynb** in sidebar. 

- **Locally**
```bash
cd results
conda activate rtss25-telescope-search
jupyter notebook visualize_results.ipynb
```
All required `.csv` results are stored in `results/recomputed_results`.

---

## 4 Extensibility of Experiments

### 4.1 Running a Single Algorithm with Designated Tiling, Speed, and Time Budget

We provide a command-line interface (CLI) to run **one algorithm at a time** on a given problem instance.
As before, a Gurobi license is needed only to run the ILP-based algorithm.

#### Usage

```bash
cd build
./op <file> <budget> <alg> [slew_rate=50]
```

#### Arguments

\<file>: path to a CSV file specifying the set of tiles, with probabilities, for the problem instance (see Input format below).

\<budget>: total time budget (slew + dwell).

\<alg>: one of greedy | genetic | gcp | ilp.

\[slew_rate] (optional): slew time per angular distance unit (default 50 degree per second).


#### Examples
Greedy on a small instance, budget = 50 (degree/second), default slew_rate
```bash
./op ../data/small/filtered_GW191105_143521_7dt_separate.csv 50 greedy
```

Genetic on a large instance with a larger budget and custom slew_rate
```bash
./op ../data/large/filtered_GW191103_012549_7dt.csv 500 genetic 30
```

GCP on a large instance
```bash
./op ../data/large/filtered_GW191109_010717_7dt.csv 200 gcp
```

ILP (Gurobi) on a small instance
```bash
./op ../data/small/filtered_GW191105_143521_7dt_separate.csv 50 ilp 40
```
##### Input Format
The input is a CSV tiling file with the following first five columns in order (additional columns after these and are ignored):

(1) Rank — rank of filtered tiles (integer)

(2) index — tile ID in the original tiling (integer; not used but kept for index alignment)

(3) RA — right ascension of the tile center in degrees [0, 360)

(4) Dec — declination of the tile center in degrees [-90, 90]

(5) Probability — probability (or likelihood/prize, double) of this tile 

Header example:
```bash
Rank,index,RA,Dec,Probability[, ...]
```

##### Output

The result will be printed to stdout:
- Path (node sequence)
- Number of nodes
- Sum probability (FoM)
- Wall-clock runtime (seconds)

### 4.2 Add Your Algorithm
All algorithms share the following interface:
```cpp
std::vector<int> my_algo(const std::vector<std::vector<double>>& costs,
                         const std::vector<double>& prize,
                         double budget, int s, int t);
```


Once have your algoritm added, add a dispatch branch in main_custom.cpp:

```cpp
else if (alg_lc == "myalgo") {
    (void)run_and_report("MyAlgo", [&]{
        return my_algo(costs, probability, eff_budget, start_idx, end_idx);
    }, costs, probability, ranks, padding);
}
```
Rebuild:
```bash
cd build && make -j
./op Data/small/foo.csv 100 myalgo
```


