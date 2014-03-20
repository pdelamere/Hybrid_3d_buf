#F77 = mpif90 -i4 -r4 -O2 -byteswapio
F77 = mpif90 -i4 -real-size 32 -O4

FILES = dimensions.f inputs.f global.f misc.f  boundary.f gutsp_dd.f  gutsp_dd.f  gutsf.f part_init.f gutsp_buf.f chem_rates.f maind.f
INCLUDE = incurv.h para.h
OBJECTS = dimensions.o inputs.o global.o misc.o boundary.o gutsp_dd.o   gutsf.o   part_init.o initial.o gutsp_buf.o chem_rates.o  maind.o

hybrid:	$(OBJECTS) 
	$(F77) -o hybrid $(OBJECTS) 

clean:
	rm *.o hybrid 

maind.o:maind.f $(INCLUDE);$(F77) -c maind.f
gutsf.o:gutsf.f $(INCLUDE);$(F77) -c gutsf.f
gutsp_dd.o:gutsp_dd.f $(INCLUDE);$(F77) -c gutsp_dd.f
misc.o:misc.f $(INCLUDE);$(F77) -c misc.f
boundary.o:boundary.f $(INCLUDE);$(F77) -c boundary.f
part_init.o:part_init.f $(INCLUDE);$(F77) -c part_init.f
initial.o:initial.f $(INCLUDE);$(F77) -c initial.f
gutsp_buf.o:gutsp_buf.f $(INCLUDE);$(F77) -c gutsp_buf.f
chem_rates.o:chem_rates.f $(INCLUDE);$(F77) -c chem_rates.f
dimensions.o:dimensions.f $(INCLUDE);$(F77) -c dimensions.f
inputs.o:inputs.f $(INCLUDE);$(F77) -c inputs.f
global.o:global.f $(INCLUDE);$(F77) -c global.f

