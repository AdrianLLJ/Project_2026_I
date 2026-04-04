# =========================
# Makefile MPI – Parallel Build
# =========================

# Compiler
FC := mpifort

# Project root
PAR_DIR := src/parallel

# Build directories
BUILD_DIR := build
PAR_OBJ_DIR := $(BUILD_DIR)/obj/parallel
MOD_OUT_DIR := $(BUILD_DIR)/modules

# Source directories
MOD_DIR_INIT := $(PAR_DIR)/initialization
MOD_DIR_ENER := $(PAR_DIR)/energy
MOD_DIR_MC   := $(PAR_DIR)/MC_update

# Compiler flags
PAR_FLAGS := -O3 \
             -I$(MOD_DIR_INIT) \
             -I$(MOD_DIR_ENER) \
             -I$(MOD_DIR_MC) \
             -J$(MOD_OUT_DIR)

# Executable
PAR_EXEC := $(BUILD_DIR)/mainglobal_mpi.exe

# Source files
PAR_SRC := $(MOD_DIR_INIT)/constants.f90 \
           $(MOD_DIR_INIT)/system.f90 \
           $(MOD_DIR_INIT)/io_module.f90 \
           $(MOD_DIR_INIT)/init_config.f90 \
           $(MOD_DIR_INIT)/centerPolymer.f90 \
           $(MOD_DIR_ENER)/nonbonded.f90 \
           $(MOD_DIR_ENER)/bonded.f90 \
           $(MOD_DIR_ENER)/energy.f90 \
           $(MOD_DIR_MC)/MC_proposal.f90 \
           $(MOD_DIR_MC)/MC_loop.f90 \
           $(PAR_DIR)/mainglobal.f90

# Object files (mirroring folder structure)
PAR_OBJ := $(patsubst $(PAR_DIR)/%.f90, $(PAR_OBJ_DIR)/%.o, $(PAR_SRC))

# =========================
# Build executable
# =========================
$(PAR_EXEC): $(PAR_OBJ)
	@mkdir -p $(dir $@)
	$(FC) -o $@ $(PAR_OBJ)

# =========================
# Pattern rule for object files
# =========================
$(PAR_OBJ_DIR)/%.o: $(PAR_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_OUT_DIR)
	$(FC) $(PAR_FLAGS) -c $< -o $@

# =========================
# Clean build files
# =========================
clean:
	rm -rf $(BUILD_DIR)