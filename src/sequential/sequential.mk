# Authors: Jonathan, Adrián, Adrià

# Specific compiler
FC := gfortran

# Directory for this module
SEQ_DIR := src/sequential

# Build dir
BUILD_DIR := build

# Output dirs for sequential
SEQ_OBJ_DIR := $(BUILD_DIR)/obj/sequential
MOD_OUT_DIR := $(BUILD_DIR)/modules

# PATH TO MODULES (Relative to root)
MOD_DIR_INIT := $(SEQ_DIR)/initialization
MOD_DIR_ENER := $(SEQ_DIR)/energy
MOD_DIR_MC := $(SEQ_DIR)/MC_update

# Compiler flags specific to sequential
SEQ_FLAGS := -I$(MOD_DIR_INIT) \
             -I$(MOD_DIR_ENER) \
             -I$(MOD_DIR_MC)   \
             -J$(MOD_OUT_DIR)

# Executable name
SEQ_EXEC := $(BUILD_DIR)/mainglobal_seq.exe

# Source files
SEQ_SRC := $(MOD_DIR_INIT)/constants.f90 \
           $(MOD_DIR_INIT)/system.f90 \
           $(MOD_DIR_INIT)/io_module.f90 \
           $(MOD_DIR_INIT)/init_config.f90 \
           $(MOD_DIR_INIT)/centerPolymer.f90 \
           $(MOD_DIR_ENER)/nonbonded.f90 \
           $(MOD_DIR_ENER)/bonded.f90 \
           $(MOD_DIR_ENER)/energy.f90 \
           $(MOD_DIR_MC)/MC_proposal.f90 \
           $(MOD_DIR_MC)/MC_loop.f90 \
           $(SEQ_DIR)/mainglobal.f90

# Object files (in build directory)
SEQ_OBJ := $(patsubst $(SEQ_DIR)/%.f90, $(SEQ_OBJ_DIR)/%.o, $(SEQ_SRC))

# Rules specific to the sequential module
$(SEQ_EXEC): $(SEQ_OBJ)
	@mkdir -p $(dir $@)
	$(FC) -o $@ $(SEQ_OBJ)

# Pattern rule for compiling Fortran objects in this module
$(SEQ_OBJ_DIR)/%.o: $(SEQ_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_OUT_DIR)
	$(FC) $(SEQ_FLAGS) -c $< -o $@

# Local clean target
clean-seq:
	rm -rf $(SEQ_OBJ_DIR) $(SEQ_EXEC) $(MOD_OUT_DIR)/*.mod