#F77 = mpif90 -i4 -r4 -O2 -byteswapio
FC=mpif90

# Release
FFLAGS=-mcmodel=medium -O3 -implicitnone -warn truncated_source -warn errors

# Maximum Debug
#FFLAGS=-mcmodel=medium -O0 -implicitnone -init=snan,arrays -g -traceback -check all -warn all -warn errors -warn stderrors -std90

# Normal Debug
#FFLAGS=-mcmodel=medium -O0 -implicitnone -init=snan -g -traceback -check arg_temp_created -check bounds -check uninit -warn truncated_source -warn errors

# check fpe
#FFLAGS=-mcmodel=medium -O0 -implicitnone -fpe0 -g -traceback -check bounds -check uninit -warn truncated_source -warn errors

# check uninit
#FFLAGS=-mcmodel=medium -O0 -implicitnone -init=snan,arrays -g -traceback -check bounds -check uninit -warn truncated_source -warn errors

FILES = dimensions.f random_utils.f inputs.f global.f misc.f  boundary.f grid_interp.f gutsp_dd.f  gutsp_dd.f  gutsf.f part_init.f gutsp_buf.f chem_rates.f maind.f 
DEBUG = -check all -g -warn
INCLUDE = incurv.h para.h
OBJECTS = dimensions.o random_utils.o inputs.o global.o misc.o boundary.o grid_interp.o gutsp_dd.o   gutsf.o   part_init.o initial.o gutsp_buf.o chem_rates.o maind.o

hybrid:	$(OBJECTS) 
	$(FC) $(FFLAGS) -o hybrid $(OBJECTS) 

debug: $(OBJECTS) 
	$(FC) $(FFLAGS) -o hybrid_d $(OBJECTS) $(DEBUG)

clean:
	rm *.o *.mod hybrid 

maind.o:maind.f global.o dimensions.o inputs.o initial.o misc.o gutsp_dd.o gutsp_buf.o gutsf.o part_init.o grid_interp.o chem_rates.o $(INCLUDE)
	$(FC) $(FFLAGS) -c maind.f

gutsf.o:gutsf.f global.o boundary.o grid_interp.o  $(INCLUDE)
	$(FC) $(FFLAGS) -c gutsf.f

gutsp_dd.o:gutsp_dd.f global.o misc.o boundary.o grid_interp.o $(INCLUDE)
	$(FC) $(FFLAGS) -c gutsp_dd.f

misc.o:misc.f random_utils.o global.o boundary.o $(INCLUDE)
	$(FC) $(FFLAGS) -c misc.f

boundary.o:boundary.f global.o $(INCLUDE)
	$(FC) $(FFLAGS) -c boundary.f

part_init.o:part_init.f global.o dimensions.o misc.o gutsp_dd.o $(INCLUDE)
	$(FC) $(FFLAGS) -c part_init.f

initial.o:initial.f global.o inputs.o chem_rates.o boundary.o $(INCLUDE)
	$(FC) $(FFLAGS) -c initial.f

gutsp_buf.o:gutsp_buf.f dimensions.o global.o misc.o part_init.o $(INCLUDE)
	$(FC) $(FFLAGS) -c gutsp_buf.f

chem_rates.o:chem_rates.f global.o inputs.o gutsp_dd.o grid_interp.o $(INCLUDE)
	$(FC) $(FFLAGS) -c chem_rates.f

dimensions.o:dimensions.f $(INCLUDE)
	$(FC) $(FFLAGS) -c dimensions.f

random_utils.o:random_utils.f
	$(FC) $(FFLAGS) -c random_utils.f

inputs.o:inputs.f dimensions.o $(INCLUDE)
	$(FC) $(FFLAGS) -c inputs.f

global.o:global.f inputs.o dimensions.o $(INCLUDE)
	$(FC) $(FFLAGS) -c global.f

grid_interp.o:grid_interp.f global.o dimensions.o boundary.o $(INCLUDE)
	$(FC) $(FFLAGS) -c grid_interp.f
