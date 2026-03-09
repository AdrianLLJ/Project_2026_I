# Author: Adrià Brú

# Global config and mode selection
# Default mode is sequential. Override via CLI: "make MODE=parallel"
# Default number of cores is 4. Override via CLI: "make CORES=8"
MODE ?= sequential
CORES ?= 4

# Build directory
BUILD_DIR := build
OUT_DIR := results
DATA_DIR := $(OUT_DIR)/data
PLOT_DIR := $(OUT_DIR)/plots
INPUT_FILE := input.dat

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
    RUN_CMD     := cd $(PAR_DIR) && mpirun -np $(CORES) ./$(TARGET_EXEC) $(DATA_DIR) < $(INPUT_FILE)
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
	mkdir -p $(DATA_DIR) $(PLOT_DIR)
	$(RUN_CMD)
	@echo "[$(MODE)] Simulation finished. Data is in $(DATA_DIR)/"

# "make run" will only run if the data doesn't exist yet
run: $(PRIMARY_DATA)

# Plotting
plot: run .venv
	@echo "Generating figures..."
	./.venv/bin/python scripts/main.py $(DATA_DIR) $(PLOT_DIR)
	@echo "Plots saved to results/plots/"

# Create venv if not present
.venv:
	python3 -m venv ./.venv
	./.venv/bin/pip install numpy matplotlib

# Clean
clean: clean_data clean_build
	@echo "Cleaning project..."

clean_data:
	@echo "Cleaning results..."
	rm -rf $(OUT_DIR)

clean_build:
	@echo "Cleaning build..."
	rm -rf $(BUILD_DIR)