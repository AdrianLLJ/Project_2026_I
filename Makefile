# Author: Adrià Brú

# Global config and mode selection
# Default mode is sequential. Override via CLI: "make MODE=parallel"
# Default number of cores is 4. Override via CLI: "make CORES=8"
MODE ?= sequential
CORES ?= 4

# Build directory
BUILD_DIR := build
OUT_DIR := results/data
INPUT_FILE := input.dat

# Module inclusion (dynamic)
ifeq ($(MODE),sequential)
    # Include sequential rules
    include src/sequential/sequential.mk
    
    # Mappings
    TARGET_EXEC := $(SEQ_EXEC)
    RUN_CMD     := ./$(TARGET_EXEC) $(OUT_DIR) < $(INPUT_FILE)
    CLEAN_CMD   := $(MAKE) -f src/sequential/sequential.mk clean-seq

else ifeq ($(MODE),parallel)
    # Include parallel rules (for the future)
    include src/parallel/parallel.mk
    
    # Mappings
    TARGET_EXEC := $(PAR_EXEC)
    RUN_CMD     := cd $(PAR_DIR) && mpirun -np $(CORES) ./$(TARGET_EXEC) $(OUT_DIR) < $(INPUT_FILE)
    CLEAN_CMD   := $(MAKE) -f src/parallel/parallel.mk clean-par 

else
    $(error "Invalid MODE. Use MODE=sequential or MODE=parallel")
endif

# =Phony targets
.PHONY: all build run plot clean

all: build

# Compile TARGET_EXEC
build: $(TARGET_EXEC)
	@echo "[$(MODE)] Build complete: $(TARGET_EXEC)"

# The run target executes whatever RUN_CMD was set to
run: build
	@echo "[$(MODE)] Running simulation..."
	$(RUN_CMD)
	@echo "[$(MODE)] Simulation finished. Data is in results/data/"

# Plotting
plot: run
	@echo "Generating figures..."
	python3 scripts/main.py results/data/main.csv results/plots/
	@echo "Plots saved to results/plots/"

# Clean
clean:
	@echo "[$(MODE)] Cleaning project..."
	rm -rf $(BUILD_DIR)