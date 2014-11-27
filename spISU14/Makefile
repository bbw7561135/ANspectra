# the compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -Wall -g -fPIC
FCFLAGS+= -fdefault-real-8
FCFLAGS+= -ffpe-trap=invalid,zero,overflow
#FCFLAGS+= -O2
# tell the compiler where your .mod files are located
#FCFLAGS+= -I/usr/include
# tell the compiler where to put compiled modules
FCFLAGS+= -J./bin

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = spISU14

# "make" builds all
all: $(PROGRAMS)

splectra.o: +PhysMathConstants.o
spISU14: +PhysMathConstants.o splectra.o AN_spectra.o AN_SHE.o NR_Spline.o

BIN=bin/
SRC=src/
LIB=lib/

SOURCES=$(addprefix $(SRC),$^)
OBJECTS=$(addprefix $(BIN),$^)
LIBRARIES=$(addprefix $(LIB),$^)
GARBAGE=$(foreach dir,$(BIN) ./,$(wildcard $(dir)*.o $(dir)*.mod  $(dir)*.MOD $(dir)fort.*))

# General rule for building shared library
%.so:
	ld -shared  $(OBJECTS) -o $(LIB)$@ 
	# $(FC) -shared $(FCFLAGS) -fPIC -o $@ $^ $(LDFLAGS)

# General rule for building prog from prog.o; $^ (GNU extension) is used in order to list additional object files on which the executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is used in order to list only the first prerequisite (the source file) and not the additional prerequisites such as module or include files
%.o: $(SRC)%.for
	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@
%.o: $(SRC)%.f90
	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@

# Utility targets
clean:
	$(RM) $(GARBAGE)
	$(RM) $(PROGRAMS)
	$(RM) $(LIBRARIES)

# DO NOT DELETE THIS LINE! - make depend depends on it
