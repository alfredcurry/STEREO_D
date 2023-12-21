#the compiler
CC = g++

#compiler flags 
# needs to be version 11 for boost to work
CFLAGS = -std=c++11
#include directory of some other codes
IDIR = ./commoncodes
PDIR = -I ../packages
INCDIR = -I ./commoncodes/ -I ./codefiles $(PDIR)
#codefile directory
CDIR = ./codefiles
FDIR = ./frontcodefiles
FTDIR = ${FDIR}/tests
#make directory
MKDIR_P = mkdir -p
ODIR = ./objects

TDIR = ./Tests
TDIRs = ${TDIR}/logs ${TDIR}/RunFiles/Ideal ${TDIR}/RunFiles/Ideal_Tab

#group codefiles
nablafiles = $(addprefix $(CDIR)/, Boost_solver.hpp Fluxes.hpp adiabatnab.hpp calc_nabla.hpp ConductiveNab.hpp)
backfiles = $(addprefix $(CDIR)/, Jacob_Maker.hpp converge.hpp calc_nabla.hpp StructStruct.hpp PhysicalConsts.h EquationOfState.hpp MeltFraction.hpp Viscosity.hpp Boundary_Layer.hpp MassTerm.hpp Conductivity.hpp)
one_per_t_files = $(addprefix $(CDIR)/, Core_Latent.hpp Radioactivity.hpp)
#group objects
backobjects = $(addprefix $(ODIR)/, converge.o Jacob_Maker.o calc_nabla.o )
frontobjects = $(addprefix $(ODIR)/, TimeStructure.o GeneralisedStructure.o IrrEnvStructure2.o Struct_test_boost.o Ideal.o Tab_Ideal.o Poly_Time_test_boost.o Struct_test_boost_Irr.o Seager.o Stix.o Mix.o Abe.o StixStructure.o SeaStructure.o Anal_BC.o Tab_BC.o)

#executables
Idealexes = TimeStruct GenStruct IrrEnvStruct2 SeaStruct
Idealexes_tab = GenStruct_Tab IrrEnvStruct2_Tab
Rockyexes = StixStruct MixStruct BC_Struct AbeStruct
tests = Struct_test_boost Poly_Time_test_boost Struct_test_boost_Irr Struct_test_boost_Tab Poly_Time_test_boost_Tab Struct_test_boost_Irr_Tab
allexes = $(Idealexes) $(Idealexes_tab) $(Rockyexes)

all: directories test_directory $(Idealexes) $(Idealexes_tab) $(tests) $(Rockyexes)

test: test_directory $(tests) $(addprefix run-,$(tests))

directories: 
	$(MKDIR_P) ${ODIR}
test_directory:
	$(MKDIR_P) ${TDIRs}

SeaStruct: $(ODIR)/SeaStructure.o $(ODIR)/Seager.o $(ODIR)/Anal_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o SeaStruct $(ODIR)/SeaStructure.o $(ODIR)/Seager.o $(ODIR)/Anal_BC.o $(backobjects)

StixStruct: $(ODIR)/StixStructure.o $(ODIR)/Stix.o $(ODIR)/Anal_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o StixStruct $(ODIR)/StixStructure.o $(ODIR)/Stix.o $(ODIR)/Anal_BC.o $(backobjects)

MixStruct: $(ODIR)/StixStructure.o $(ODIR)/Mix.o $(ODIR)/Anal_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o MixStruct $(ODIR)/StixStructure.o $(ODIR)/Mix.o $(ODIR)/Anal_BC.o $(backobjects)

AbeStruct: $(ODIR)/StixStructure.o $(ODIR)/Abe.o $(ODIR)/Anal_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o AbeStruct $(ODIR)/StixStructure.o $(ODIR)/Abe.o $(ODIR)/Anal_BC.o $(backobjects)

BC_Struct: $(ODIR)/StixStructure.o $(ODIR)/Mix.o $(ODIR)/Tab_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o BC_Struct $(ODIR)/StixStructure.o $(ODIR)/Mix.o $(ODIR)/Tab_BC.o $(backobjects)

IrrEnvStruct2_Tab: $(ODIR)/IrrEnvStructure2.o $(ODIR)/Tab_Ideal.o $(ODIR)/Anal_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o IrrEnvStruct2_Tab $(ODIR)/IrrEnvStructure2.o $(ODIR)/Tab_Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

IrrEnvStruct2: $(ODIR)/IrrEnvStructure2.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)
	$(CC) $(CFLAGS) -o IrrEnvStruct2 $(ODIR)/IrrEnvStructure2.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

TimeStruct: $(ODIR)/TimeStructure.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Ideal.o
	$(CC) $(CFLAGS) -o TimeStruct $(ODIR)/TimeStructure.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

GenStruct_Tab: $(ODIR)/GeneralisedStructure.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Tab_Ideal.o
	$(CC) $(CFLAGS) -o GenStruct_Tab $(ODIR)/GeneralisedStructure.o $(ODIR)/Tab_Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

GenStruct: $(ODIR)/GeneralisedStructure.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Ideal.o
	$(CC) $(CFLAGS) -o GenStruct $(ODIR)/GeneralisedStructure.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

Struct_test_boost: $(ODIR)/Struct_test_boost.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Ideal.o
	$(CC) $(CFLAGS) -o Struct_test_boost $(ODIR)/Struct_test_boost.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

Struct_test_boost_Tab: $(ODIR)/Struct_test_boost.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Tab_Ideal.o
	$(CC) $(CFLAGS) -o Struct_test_boost_Tab $(ODIR)/Struct_test_boost.o $(ODIR)/Tab_Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

Poly_Time_test_boost: $(ODIR)/Poly_Time_test_boost.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Ideal.o
	$(CC) $(CFLAGS) -o Poly_Time_test_boost $(ODIR)/Poly_Time_test_boost.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

Poly_Time_test_boost_Tab: $(ODIR)/Poly_Time_test_boost.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Tab_Ideal.o
	$(CC) $(CFLAGS) -o Poly_Time_test_boost_Tab $(ODIR)/Poly_Time_test_boost.o $(ODIR)/Tab_Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

Struct_test_boost_Irr: $(ODIR)/Struct_test_boost_Irr.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Ideal.o
	$(CC) $(CFLAGS) -o Struct_test_boost_Irr $(ODIR)/Struct_test_boost_Irr.o $(ODIR)/Ideal.o $(ODIR)/Anal_BC.o $(backobjects)

Struct_test_boost_Irr_Tab: $(ODIR)/Struct_test_boost_Irr.o $(backobjects) $(ODIR)/Anal_BC.o $(ODIR)/Tab_Ideal.o
	$(CC) $(CFLAGS) -o Struct_test_boost_Irr_Tab $(ODIR)/Struct_test_boost_Irr.o $(ODIR)/Tab_Ideal.o $(ODIR)/Anal_BC.o $(backobjects)
	
$(ODIR)/Seager.o: $(CDIR)/SeagerEquationOfState.cpp $(backfiles) $(CDIR)/EquationOfStateSea.hpp
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Stix.o: $(CDIR)/Tab_Stix_EoS.cpp $(backfiles) $(IDIR)/2Dinterpolator.hpp $(IDIR)/File_error.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(ODIR)/Mix.o: $(CDIR)/Mix_EoS.cpp $(backfiles) $(IDIR)/2Dinterpolator.hpp $(IDIR)/File_error.h $(IDIR)/Vector_reader.h $(CDIR)/MeltFunctions.hpp $(CDIR)/MeltViscosity_Arrhenius.hpp 
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@  

$(ODIR)/Abe.o: $(CDIR)/Abe_simple_EoS.cpp $(backfiles) $(IDIR)/File_error.h $(IDIR)/Vector_reader.h $(CDIR)/MeltFunctions.hpp $(CDIR)/MeltViscosity_Arrhenius.hpp 
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@  

$(ODIR)/Tab_BC.o: $(CDIR)/Tabulated_BC.cpp $(CDIR)/BoundaryCondition.hpp $(CDIR)/PhysicalConsts.h $(IDIR)/2Dinterpolator.hpp $(IDIR)/File_error.h $(IDIR)/Vector_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Anal_BC.o: $(CDIR)/Analytic_BC.cpp $(CDIR)/fitted_BC.hpp $(CDIR)/BoundaryCondition.hpp $(CDIR)/PhysicalConsts.h $(CDIR)/SurfaceTemp.hpp $(IDIR)/Vector_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@   

$(ODIR)/Ideal.o: $(CDIR)/IdealEquationOfState.cpp $(backfiles) $(CDIR)/EquationOfStateIdeal.hpp 
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Tab_Ideal.o: $(CDIR)/Tab_Ideal_EoS.cpp $(backfiles) $(IDIR)/2Dinterpolator.hpp $(IDIR)/File_error.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/converge.o: $(CDIR)/converge.cpp $(backfiles) $(IDIR)/block_tridiag_solve.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Jacob_Maker.o: $(CDIR)/Jacob_Maker.cpp $(backfiles) $(CDIR)/SurfacePressure.hpp
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/calc_nabla.o: $(CDIR)/calc_nabla.cpp $(backfiles) $(nablafiles)
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/IrrEnvStructure2.o: $(FDIR)/IrrEnvStructure2.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/TimeStructure.o: $(FDIR)/TimeStructure.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/GeneralisedStructure.o: $(FDIR)/GeneralisedStructure.cc $(backfiles) $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(ODIR)/SeaStructure.o: $(FDIR)/StructureSea.cc $(backfiles) $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(ODIR)/StixStructure.o: $(FDIR)/StructureStix.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h $(CDIR)/MassLoss.hpp $(CDIR)/Kang_analytical.hpp $(IDIR)/2Dinterpolator.hpp 
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(ODIR)/Struct_test_boost.o: $(FTDIR)/Struct_test_boost.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Poly_Time_test_boost.o: $(FTDIR)/Poly_Time_test_boost.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

$(ODIR)/Struct_test_boost_Irr.o: $(FTDIR)/Struct_test_boost_Irr.cc $(backfiles) $(CDIR)/timestep.hpp $(one_per_t_files) $(CDIR)/Print_and_save.hpp $(IDIR)/Interpolator.hpp $(IDIR)/Matrix_reader.h
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@ 

run-%: %
	-./$^ > ./Tests/logs/$^.log

.PHONY : clean
clean:
	$(RM) $(backobjects) $(frontobjects) $(allexes) $(tests) $(addsuffix /* , ${TDIRs})

