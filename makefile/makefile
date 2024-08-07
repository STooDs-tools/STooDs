#========================================================================
# Makefile to build STooDs with gfortran
#========================================================================
#

# Define the Fortran Compiler and its options
FC  = gfortran
FLAGS = -fcheck=all -Wall -O3 # DEBUG: -fcheck=all -Wall -g -O0

# Define directories
DMSL_DIR = ../../miniDMSL/src/
BaM_DIR = ../../BaM/src/Models/
BMSL_DIR = ../../BMSL/src/
STooDs_DIR = ../src/


# Define DMSL files
DMSL = \
	$(DMSL_DIR)kinds_dmsl_kit.f90              \
	$(DMSL_DIR)types_dmsl_kit.f90              \
	$(DMSL_DIR)utilities_dmsl_kit.f90          \
	$(DMSL_DIR)mt19937-64.f95                  \
	$(DMSL_DIR)uniran1_minidmsl_mod.f90        \
	$(DMSL_DIR)linalg_dmsl_kit.f90             \
	$(DMSL_DIR)numerix_dmsl_kit.f90

# Define BMSL files
BMSL = \
	$(BMSL_DIR)Geodesy_tools.f90                 \
	$(BMSL_DIR)EmpiricalStats_tools.f90          \
	$(BMSL_DIR)Distribution_tools.f90            \
	$(BMSL_DIR)MultivariateDistribution_tools.f90

# Define BaM files
BaM = \
	$(BaM_DIR)TextFile_model.f90

# Define BaM files
STooDs = \
	$(STooDs_DIR)STOODS_types.f90          \
	$(STooDs_DIR)STOODS_formulas.f90       \
	$(STooDs_DIR)STOODS_inference.f90      \
	$(STooDs_DIR)STOODS_mcmc.f90           \
	$(STooDs_DIR)Main.f90 

all:build_exe

build_exe: compile_bam
	$(FC)  $(FLAGS) $(DMSL) $(BMSL) $(BaM) $(STooDs) -o STooDs

compile_stoods: compile_bam
	$(FC) $(FLAGS) -c $(stoods)

compile_bam: compile_BMSL
	$(FC) $(FLAGS) -c $(BaM)

compile_BMSL: compile_dmsl
	$(FC) $(FLAGS) -c $(BMSL)

compile_dmsl:
	$(FC) $(FLAGS) -c $(DMSL)





