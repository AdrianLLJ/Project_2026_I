# Compiler
FC = gfortran

# PATH TO MODULES
MOD_DIR_INIT = ./initialization
MOD_DIR_ENER = ./energy
MOD_DIR_MOD = ./modules

# Compiler flags
FLAGS = -I$(MOD_DIR_INIT) -I$(MOD_DIR_ENER) -I$(MOD_DIR_MOD) -J$(MOD_DIR_MOD)

# Executable name
EXEC = mainglobal.exe

# Source files
SRC = $(MOD_DIR_INIT)/constants.f90 \
      $(MOD_DIR_INIT)/system.f90 \
      $(MOD_DIR_INIT)/io_module.f90 \
      $(MOD_DIR_INIT)/init_config.f90 \
      $(MOD_DIR_INIT)/centerPolymer.f90 \
	  $(MOD_DIR_ENER)/nonbonded.f90 \
	   $(MOD_DIR_ENER)/bonded.f90 \
	  $(MOD_DIR_ENER)/energy.f90 \
	  mainglobal.f90

# Object files
OBJ = $(SRC:.f90=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) -o $@ $(OBJ)

%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC) *.mod $(MOD_DIR_INIT)/*.mod \
							   $(MOD_DIR_ENER)/*.mod \
							   $(MOD_DIR_MOD)/*.mod

.PHONY: all clean