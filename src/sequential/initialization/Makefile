# Makefile for the polymer initialization program

# Compiler
FC = gfortran

# Executable name
EXEC = polymer

# Source files
SRC = constants.f90 system.f90 io_module.f90 init_config.f90 main.f90

# Object files (replace .f90 with .o)
OBJ = $(SRC:.f90=.o)

# Default target
all: $(EXEC)

# Link the executable from object files
$(EXEC): $(OBJ)
	$(FC) -o $@ $(OBJ)

# Generic rule to compile .f90 to .o
%.o: %.f90
	$(FC) -c $<

# Clean compiled files
clean:
	rm -f $(OBJ) $(EXEC)

.PHONY: all clean
