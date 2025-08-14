# Artifact for RTSS 2025
This repository contains the source code, datasets, and analysis tools supporting the paper *"Probabilistic Response-Time-Aware Search for Transient Astrophysical Phenomena."*

## System Requirements

- Linux-based OS
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
  - Small problem instances incorporating worst-case execution time (WCET)
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
│   ├── large_wcet/     # Subset of large-scale maps used in WCET-aware experiments
│   └── wcet/           # Worst-case execution time measurements

├── results/
│
│   ├── figures/                   # Figures 7–10 in paper
|
│   ├── precomputed_results/       # Precomputed outputs used in the paper
│   │   ├── small/                 # FoM results (non-ILP) for small instances
│   │   ├── small_gurobi/          # ILP-based FoM results for small instances
│   │   ├── large/                 # FoM results for large instances
│   │   ├── small_with_wcet/       # WCET-aware FoM for small instances
│   │   ├── multi_deadline/        # Multi-deadline scenario results
│   │   └── analysis_precomputed_results.ipynb  # Notebook for Figures 7–10 in paper
│
│   ├── recomputed_results/        # Outputs from re-running experiments
│   │   ├── small/                 
│   │   ├── small_with_wcet/       
│   │   ├── large/                 
│   │   ├── multi_deadline/        
│   │   └── analysis_recomputed_results.ipynb  # Visualization of recomputed results
│
│   ├── run_small_instances.py         # Batch run script for small instances
│   ├── run_large_instances.py         # Batch run script for large instances
│   ├── run_instances_wcet.py          # Batch run script (WCET-aware)
│   └── run_multi_deadline.py          # Batch run script for multi-deadline setup

├── rtss25-sky-tiling.yml               # Python dependencies for tiling code
├── rtss25-telescope-search.yml         # Python dependencies
└── README.md


```

## Reproducing Paper Figures
### 1. Python Environment Setup

We recommend setting up a [conda](https://docs.conda.io/en/latest/) environment for Python.

If you do not have conda installed locally:

```bash
cd path/to/Telescope-Search-Problem
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

### 2. Reproducing Result in a Jupyter notebook

To reproduce all figures for experiment results:

```bash
conda activate rtss25-telescope-search
cd ~/Telescope-Searching-Problem/results/precomputed_results
jupyter notebook analysis_precomputed_results.ipynb
```
The notebook includes:
- Average percentage deviation from ILP baseline across deadlines for small instances(Fig.7)
- Percentage deviation in FoM from the GCP baseline and Computation time vs Deadline presented in log-scaled plot for large instances (Fig. 8)
- FoM percentage deviation from GCP after accounting for WCET. for five small instances (Fig 9)
- Expected FoM progression over the multi deadlines (Fig. 10)

All required `.csv` results are precomputed and stored in `results/precomputed_results`.

The notebook saves figures to `results/figures`.

---


## Running Full Experiments

To re-run the full set of experiments (~5+ hours total runtime), you may either set up the environment locally or use the provided Docker container. Experiments can then be executed using the provided batch scripts.

---

### 1. C++ Environment Setup

### 1.0 (Preliminary, Optional) Obtaining a Gurobi License

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
6. Move or copy the file to the `Telescope-Searching-Problem` artifact repository.


### 1.1 (Option A, Preferred) Using the Provided Docker Container

#### (1) Install Docker

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

#### (2) Pull the Docker Image
```bash
sudo docker pull ghcr.io/daisy0419/rtss25-op-solver:1.0
```

#### (3) Run the Docker Container
If you have a Gurobi license on your local machine, mount the license into the container.

**License Note**: If you're using an academic license, be sure to request an Academic WLS License (floating license). Named-User Academic Licenses are not compatible with Docker containers.

```bash
cd ~/Telescope-Searching-Problem
sudo docker run --rm -it -v "./gurobi.lic:/licenses/gurobi.lic:ro" -e GRB_LICENSE_FILE=/licenses/gurobi.lic ghcr.io/daisy0419/rtss25-op-solver:1.0
```

If you do not have a Gurobi license, you can still run experiments that do not rely on ILP-based solvers:

```bash
sudo docker run --rm -it ghcr.io/daisy0419/rtss25-op-solver:1.0
```

#### (4) Build project code
The project executables have been precompiled inside the Docker image. You can proceed directly to running experiments in the next sections.

### 1.1 (Option B) Local Installation

#### (1) Gurobi Optimizer (Required)


1. Download and extract Gurobi to the directory of your choice.

```bash
wget https://packages.gurobi.com/12.0/gurobi12.0.3_linux64.tar.gz
tar xvfz gurobi12.0.3_linux64.tar.gz
```

2. Set the necessary environment variables in your shell (change `~/gurobi1203` to your preferred location)

```bash
export GUROBI_HOME=~/gurobi1203
export PATH="${GUROBI_HOME}/bin:$PATH"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib:$LD_LIBRARY_PATH"
```

**Note**: If you don't have Gurobi license and you do not plan to run experiments involving Gurobi, you still need to install Gurobi in order to compile the project code (due to build-time linking requirements).


#### (2) LEMON Graph Library (Required)

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

#### (3) Build the C++ Executables

Once all dependencies are installed, you can build the C++ project with:

```bash
mkdir build && cd build
cmake ..
make -j
```

### 2. Regenerate Sky Tilings

The sky tiles used for the small and large problem instances evaluated are already stored in the `data/small` and `data/large` directories.
Optionally, you may regenerate them.
All commands in this step are run from the `sky_tiling` directory.

```bash
cd ~/Telescope-Searching-Problem/sky_tiling
```

#### 2.1 Precompute FoV Projections (~ 10 minutes)

Each telescope's FoV is projected onto the unit sphere to precompute tile boundaries.

```bash
conda activate rtss25-sky-tiling
bash precompute_tile_maps.sh
```

This step is optional. 
Files already exist in `sky_tiling/tile_center_files` and `sky_tiling/tile_pixel_maps`.

Recomputed files appear in those same directories, with file names containing `recomputed`.

#### 2.2 Flatten HealPix (~ 10 seconds)

Multi-order hierarchical HEALPix maps in `data/ligo_healpix`
need to be flattened into a non-hierarchical format. 

```bash
conda activate rtss25-telescope-search
python flatten_healpix.py
```

This step is optional.
Files already exist in `data/ligo_healpix_flattened`;
they will be overwritten by this script.

#### 2.3 Generate Tiles (~ 15 seconds)

Intersection of flattened HEALPix maps and precomputed tilings are intersected
to produce the tiles, with probabilities, that serve as the inputs to the search problem.

```bash
conda activate rtss25-sky-tiling
python produce_tilings.py
```

Files in `data/small` and `data/large` will be overwritten.


### 3. Recompute Results from Scratch
Each Python script corresponds to a different experiment setting. 
If you do not have a Gurobi license, you can still run experiments 3.2, 3.3, and 3.4, which do not rely on ILP solvers.

```bash
cd ~/Telescope-Searching-Problem/results
```

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

### 4. Visualizing the Results

```bash
cd results/recompute_results
jupyter notebook analysis_recomputed_results.ipynb
```
All required `.csv` results are stored in `results/recomputed_results`.

## Extensition Test

### 1 Running a Single Algorithm (CLI) with Designates Data and Budget setting

We provide a CLI to run **one algorithm at a time** on a given instance. Again, if you don't have gurobi license, you won't be able to run ILP algorithm

##### Usage

```bash
cd build
./op <file> <budget> <alg> [slew_rate=50]
```

##### Arguments

\<file>: path to the CSV instance (see Dataset format below).

\<budget>: total time budget (slew + dwell).

\<alg>: one of greedy | genetic | gcp | ilp.

\[slew_rate] (optional): slew time per angular distance unit (default 50).


##### Examples
Greedy on a small instance, budget = 50, default slew_rate
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

##### Output

The result will be printed to stdout:
- Path (node sequence)
- Number of nodes
- Sum probability (FoM)
- Wall-clock runtime (seconds)

### 2 Add Your Algorithm
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
