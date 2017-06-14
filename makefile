#F77 = mpif90 -i4 -r4 -O2 -byteswapio
FC=mpif90
FFLAGS=-mcmodel=medium -i4 -real-size 32 -O2
#FFLAGS=-mcmodel=medium -i4 -real-size 32 -g -traceback -check bounds

FILES = dimensions.f inputs.f global.f misc.f  boundary.f grid_interp.f gutsp_dd.f  gutsp_dd.f  gutsf.f part_init.f gutsp_buf.f chem_rates.f maind.f 
DEBUG = -check all -g -warn
INCLUDE = incurv.h para.h
OBJECTS = dimensions.o inputs.o global.o misc.o boundary.o grid_interp.o gutsp_dd.o   gutsf.o   part_init.o initial.o gutsp_buf.o chem_rates.o maind.o

CDF_PROG =	dummy

CDF_SRCS =	biomath_constants_mod.f90 biomath_interface_mod.f90 \
	biomath_mathlib_mod.f90 biomath_sort_mod.f90 biomath_strings_mod.f90 \
	cdf_aux_mod.f90 cdf_beta_mod.f90 cdf_binomial_mod.f90 \
	cdf_chisq_mod.f90 cdf_f_mod.f90 cdf_gamma_mod.f90 \
	cdf_nc_chisq_mod.f90 cdf_nc_f_mod.f90 cdf_nc_t_mod.f90 \
	cdf_neg_binomial_mod.f90 cdf_normal_mod.f90 cdf_poisson_mod.f90 \
	cdf_t_mod.f90 zero_finder.f90

CDF_OBJS =	biomath_constants_mod.o biomath_interface_mod.o biomath_mathlib_mod.o \
	biomath_sort_mod.o biomath_strings_mod.o cdf_aux_mod.o cdf_beta_mod.o \
	cdf_binomial_mod.o cdf_chisq_mod.o cdf_f_mod.o cdf_gamma_mod.o \
	cdf_nc_chisq_mod.o cdf_nc_f_mod.o cdf_nc_t_mod.o \
	cdf_neg_binomial_mod.o cdf_normal_mod.o cdf_poisson_mod.o cdf_t_mod.o \
	zero_finder.o

hybrid:	$(OBJECTS) CDF_all
	$(FC) $(FFLAGS) -o hybrid $(OBJECTS) $(CDF_OBJS)

debug: $(OBJECTS) 
	$(FC) $(FFLAGS) -o hybrid_d $(OBJECTS) $(DEBUG)

clean:
	rm *.o *.mod hybrid 

################# build cdflib90 #########################

CDF_all: $(CDF_PROG)

$(CDF_PROG): $(CDF_OBJS)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

biomath_interface_mod.o: biomath_constants_mod.o biomath_sort_mod.o \
	biomath_strings_mod.o
biomath_mathlib_mod.o: biomath_constants_mod.o
biomath_sort_mod.o: biomath_constants_mod.o
biomath_strings_mod.o: biomath_constants_mod.o
cdf_aux_mod.o: biomath_constants_mod.o biomath_interface_mod.o zero_finder.o
cdf_beta_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o cdf_aux_mod.o \
	zero_finder.o
cdf_binomial_mod.o: biomath_constants_mod.o cdf_aux_mod.o cdf_beta_mod.o \
	zero_finder.o
cdf_chisq_mod.o: biomath_constants_mod.o cdf_aux_mod.o cdf_gamma_mod.o \
	zero_finder.o
cdf_f_mod.o: biomath_constants_mod.o cdf_aux_mod.o cdf_beta_mod.o \
	zero_finder.o
cdf_gamma_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o cdf_aux_mod.o \
	zero_finder.o
cdf_nc_chisq_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o \
	cdf_aux_mod.o cdf_chisq_mod.o zero_finder.o
cdf_nc_f_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o cdf_aux_mod.o \
	cdf_beta_mod.o cdf_f_mod.o cdf_gamma_mod.o zero_finder.o
cdf_nc_t_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o cdf_aux_mod.o \
	cdf_beta_mod.o cdf_normal_mod.o cdf_t_mod.o zero_finder.o
cdf_neg_binomial_mod.o: biomath_constants_mod.o cdf_aux_mod.o cdf_beta_mod.o \
	zero_finder.o
cdf_normal_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o cdf_aux_mod.o
cdf_poisson_mod.o: biomath_constants_mod.o cdf_aux_mod.o cdf_gamma_mod.o
cdf_t_mod.o: biomath_constants_mod.o biomath_mathlib_mod.o cdf_aux_mod.o \
	cdf_beta_mod.o cdf_normal_mod.o zero_finder.o
zero_finder.o: biomath_constants_mod.o

##########################################################

maind.o:maind.f $(INCLUDE);$(FC) $(FFLAGS) -c maind.f
gutsf.o:gutsf.f $(INCLUDE);$(FC) $(FFLAGS) -c gutsf.f
gutsp_dd.o:gutsp_dd.f $(INCLUDE);$(FC) $(FFLAGS) -c gutsp_dd.f
misc.o:misc.f $(INCLUDE);$(FC) $(FFLAGS) -c misc.f
boundary.o:boundary.f $(INCLUDE);$(FC) $(FFLAGS) -c boundary.f
part_init.o:part_init.f $(INCLUDE) CDF_all;$(FC) $(FFLAGS) -c part_init.f
initial.o:initial.f $(INCLUDE);$(FC) $(FFLAGS) -c initial.f
gutsp_buf.o:gutsp_buf.f $(INCLUDE) CDF_all;$(FC) $(FFLAGS) -c gutsp_buf.f
chem_rates.o:chem_rates.f $(INCLUDE);$(FC) $(FFLAGS) -c chem_rates.f
dimensions.o:dimensions.f $(INCLUDE);$(FC) $(FFLAGS) -c dimensions.f
inputs.o:inputs.f $(INCLUDE) CDF_all;$(FC) $(FFLAGS) -c inputs.f
global.o:global.f $(INCLUDE);$(FC) $(FFLAGS) -c global.f
grid_interp.o:grid_interp.f $(INCLUDE);$(FC) $(FFLAGS) -c grid_interp.f
