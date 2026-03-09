# Author: Adrià Brú

# Global config and mode selection
# Default mode is sequential. Override via CLI: "make MODE=parallel"
# Default number of cores is 4. Override via CLI: "make NCORES=8"
MODE ?= sequential
NCORES ?= 4

# Build directory
BUILD_DIR := build
OUT_DIR := results
DATA_DIR := $(OUT_DIR)/data
PLOT_DIR := $(OUT_DIR)/plots
INPUT_FILE := input.dat
PYTHON_REQS := requirements.txt

# Random data file that proves the simulation run has completed
PRIMARY_DATA := $(DATA_DIR)/energy.dat

# Module inclusion (dynamic)
ifeq ($(MODE),sequential)
    # Include sequential rules
    include src/sequential/sequential.mk
    
    # Mappings
    TARGET_EXEC := $(SEQ_EXEC)
    RUN_CMD     := ./$(TARGET_EXEC) $(DATA_DIR) < $(INPUT_FILE)
    CLEAN_CMD   := $(MAKE) -f src/sequential/sequential.mk clean-seq

else ifeq ($(MODE),parallel)
    # Include parallel rules (for the future)
    include src/parallel/parallel.mk
    
    # Mappings
    TARGET_EXEC := $(PAR_EXEC)
    RUN_CMD     := cd $(PAR_DIR) && mpirun -np $(NCORES) ./$(TARGET_EXEC) $(DATA_DIR) < $(INPUT_FILE)
    CLEAN_CMD   := $(MAKE) -f src/parallel/parallel.mk clean-par 

else
    $(error "Invalid MODE. Use MODE=sequential or MODE=parallel")
endif

# Phony targets (always assumed to be outdated)
.PHONY: all build run plot clean clean_data

all: build

# Compile TARGET_EXEC
build: $(TARGET_EXEC)
	@echo "[$(MODE)] Build complete: $(TARGET_EXEC)"

# We add this to avoid rerunning the simulation every time we want to plot
# This way, the simulation will only run if the files aren't there
$(PRIMARY_DATA): $(TARGET_EXEC) $(INPUT_FILE)
	@echo "[$(MODE)] Running simulation..."
	mkdir -p $(DATA_DIR)
	$(RUN_CMD)
	@echo "[$(MODE)] Simulation finished. Data is in $(DATA_DIR)/"

# "make run" will only run if the data doesn't exist yet
run: $(PRIMARY_DATA)

# Plotting
plot: run .venv
	mkdir -p $(PLOT_DIR)
	@echo "Generating figures..."
	./.venv/bin/python scripts/main.py $(DATA_DIR) $(PLOT_DIR)
	@echo "Plots saved to results/plots/"

plot-cluster: sync-down .venv
	@echo "Locating the most recent cluster job..."
	@LATEST_DIR=$$(ls -td results/cluster_archive/*/ 2>/dev/null | head -1); \
	if [ -z "$$LATEST_DIR" ]; then \
		echo "Error: No downloaded cluster data found."; exit 1; \
	fi; \
	echo "Plotting data from $$LATEST_DIR..."; \
	mkdir -p "$${LATEST_DIR}plots"; \
	./.venv/bin/python scripts/main.py "$${LATEST_DIR}data" "$${LATEST_DIR}plots"; \
	echo "Cluster plots successfully saved to $${LATEST_DIR}plots/"

# Create venv if not present
.venv:
	@echo "Creating local virtual environment..."
	python3 -m venv ./.venv
	./.venv/bin/pip install -r $(PYTHON_REQS)

# Clean
clean-all: clean clean-venv

clean: clean-data clean-build
	@echo "Cleaning project..."

clean-data:
	@echo "Cleaning results..."
	rm -rf $(OUT_DIR)

clean-build:
	@echo "Cleaning build..."
	rm -rf $(BUILD_DIR)

clean-venv:
	@echo "Cleaning Python venv"
	rm -rf .venv


# Cluster automation
# .env file MUST include the following variables: CLUSTER_USER, CLUSTER_HOST, CLUSTER_DIR
# Example:
# CLUSTER_USER := user_name
# CLUSTER_HOST := cerqt2.qt.ub.edu
# CLUSTER_DIR  := ~/MC_Project
-include .env

CLUSTER_USER ?= missing_user
CLUSTER_HOST ?= missing_host
CLUSTER_DIR  ?= missing_dir
CLUSTER_EMAIL ?= missing_email@ub.edu

.PHONY: sync-up cluster-run sync-down cluster-check clean-cluster

# Push local code to cluster
sync-up:
	@if [ "$(CLUSTER_USER)" = "missing_user" ]; then echo "Error: .env file missing or CLUSTER_USER not set" && exit 1; fi
	@echo "Uploading code to cluster..."
	rsync -avz --exclude '.venv' --exclude 'build' --exclude 'results' --exclude '.git' --exclude '.env' ./ $(CLUSTER_USER)@$(CLUSTER_HOST):$(CLUSTER_DIR)

# Push code and submit the job (including compilation ofc)
run-cluster: sync-up
	@echo "Compiling and submitting job on cluster..."
	ssh $(CLUSTER_USER)@$(CLUSTER_HOST) "source /etc/profile && cd $(CLUSTER_DIR) && qsub -M $(CLUSTER_EMAIL) submit.sub"
	@echo "Job submitted! Use 'make check-cluster' to check status."

# Pull results back from cluster
sync-down:
	@echo "Fetching results from cluster..."
	rsync -avz $(CLUSTER_USER)@$(CLUSTER_HOST):$(CLUSTER_DIR)/results_archive/ ./results/cluster_archive/

check-cluster:
	@echo "Checking job status..."
	ssh $(CLUSTER_USER)@$(CLUSTER_HOST) "source /etc/profile && qstat"

clean-cluster:
	@echo "Cleaning cluster project directory"
	ssh $(CLUSTER_USER)@$(CLUSTER_HOST) "rm -rf $(CLUSTER_DIR)"