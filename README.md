# Monte Carlo Polymer Simulation

This repository contains a Fortran-based Monte Carlo simulation of a polymer chain.
The project uses a `Makefile` that handles compiling, running, dependency management, and job submission to the university HPC cluster.

## Prerequisites
* **Fortran Compiler:** `gfortran`
* **Python 3:** For data visualization.
* **Make:** To run the automated pipeline.

### List of participants and tasks:

    JONATHAN BARRIENTOS ANDRADE (J-dot-Barrientos) - Initialization
    ADRIÀ BRÚ I CORTÉS (cooligula) - Main Leader / Post processing and statistics
    ADRIAN LLAMAS JARAMILLO (AdrianLLJ) - Energy

## Project structure
```
mc-polymer-simulation/
├── src/                	# Directory containing all source code.
│   ├── sequential/     	# Initial sequential code for the March 12th deliverable.
│   └── parallel/       	# The MPI version supporting replica, parallel energy, and mixed parallelism.
├── scripts/            	# Python scripts for post-processing and plotting.
├── results/             
│   ├── data/           	# Energy, torsion angle, etc.
│   ├── plots/          	# Plots of obtained data.
│   └── trajectory/     	# Monte Carlo trajectory.
├── build/					# Build directory (ignored).
├── Makefile            	# Main Makefile.
├── README.md				# Main README.
├── .gitignore				
├── .env					# Environment file (ignored): define CERQT2 Makefile variables.
├── input.dat				# Simulation input data.
├── submit.sub				# Cluster job description.
├── requirements.txt		# Python pip dependencies.
└── tests/              	# Differents tests.
```

## Local Usage

For testing, debugging, and short simulations, you can run the entire pipeline on your local machine.

### 1. Run the Simulation
To compile the Fortran code (if modified) and run the simulation locally, use:
```bash
make run
```

### 2. Generate Plots

To generate visualizations of the simulation data:
```bash
make plot
```

This command will:
1. Run the simulation (if the data is missing).
2. Create an isolated Python virtual environment (.venv).
3. Install the required dependencies (numpy, matplotlib).
4. Generate the plots and save them to results/plots/.

### 3. Cleanup
```bash
make clean      # Deletes data and build folders
make clean-all  # Deletes data, build folders, and the Python .venv
```

## Cluster Usage (High-Performance Computing)
For long production runs (e.g., $10^8$ steps), the simulation should be offloaded to the university cluster (cerqt2.qt.ub.edu). The Makefile fully automates the syncing, submission, and retrieval of these jobs.

### Step 1: Passwordless SSH Setup (One-Time)
To allow the Makefile to communicate with the cluster without prompting for a password every time, you must set up an SSH key pair.

#### 1. Generate the key locally:
```bash
ssh-keygen -t ed25519
```
#### 2. Copy the key to the cluster:
```
ssh-copy-id your_username@cerqt2.qt.ub.edu
```

### Step 2: Configure the .env File
For security, your cluster credentials are never hardcoded or uploaded to GitHub. You must create a file named .env in the root of the project.

Create a .env file with the following variables:
```bash
CLUSTER_USER=your_ub_username
CLUSTER_HOST=cerqt2.qt.ub.edu
CLUSTER_DIR=~/EIA_project
CLUSTER_EMAIL=your_email@ub.edu
```

### Step 3: Cluster Workflow Commands
Once your SSH keys and .env file are set up, you can control the cluster entirely from your local terminal.

#### Submit a Job:
```bash
make run-cluster
```

#### Check Job Status:
```bash
make check-cluster
```

#### Fetch and Plot Results:
```bash
make plot-cluster
```
Once you receive the completion email, run this command. It will download the generated data from the cluster to your local results/cluster_archive/ folder, locate the newest run, and generate the plots locally using your .venv.

#### Clean Cluster Workspace:
```bash
make clean-cluster
```