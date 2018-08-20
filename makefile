#F77 = mpif90 -i4 -r4 -O2 -byteswapio
FC=mpif90
#FFLAGS=-mcmodel=medium -O2
FFLAGS=-mcmodel=medium -i4 -real-size 32 -O0 -g -traceback -check all -check uninit -ftrapuv
#FFLAGS=-mcmodel=medium -i4 -real-size 32 -O0 -g -traceback -check all -check uninit -ftrapuv -warn all

FILES = dimensions.f inputs.f global.f misc.f  boundary.f grid_interp.f gutsp_dd.f  gutsp_dd.f  gutsf.f part_init.f gutsp_buf.f chem_rates.f maind.f 
DEBUG = -check all -g -warn
INCLUDE = incurv.h para.h
OBJECTS = dimensions.o inputs.o global.o misc.o boundary.o grid_interp.o gutsp_dd.o   gutsf.o   part_init.o initial.o gutsp_buf.o chem_rates.o maind.o

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

misc.o:misc.f global.o $(INCLUDE)
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

inputs.o:inputs.f dimensions.o $(INCLUDE)
	$(FC) $(FFLAGS) -c inputs.f

global.o:global.f inputs.o dimensions.o $(INCLUDE)
	$(FC) $(FFLAGS) -c global.f

grid_interp.o:grid_interp.f global.o dimensions.o boundary.o $(INCLUDE)
	$(FC) $(FFLAGS) -c grid_interp.f
