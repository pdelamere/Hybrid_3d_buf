      MODULE INPUTS

      USE dimensions
      USE mpi
      save

      real b0_init
      integer ion_amu
      integer mp
      real nf_init
      real dt_frac
      integer nt
      integer nout
      integer part_nout
      integer output_wait
      real vsw
      real vth
      real Ni_tot_frac
      real dx_frac
      real nu_init_frac
      real mion,q
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      real lambda_i                  !ion inertial length

c grid parameters

      real dx,dy, delz, dx_buf, dt, dtsub_init

c time stepping parameters
      PARAMETER (ntsub = 10.0)        !number of subcycle time steps

c output directory
      character(50) out_dir
c      PARAMETER (out_dir='./tmp3/')

c logical variable for restart
      logical restart
      integer mbegin,mrestart
      PARAMETER (mbegin = 0)      !mstart

c neutral cloud expansion characteristics
      real vtop,vbottom

c misc constants
      real*8 mu0,epsilon0,pi,rtod,mO,km_to_m,kboltz,melec
      real*8 tempf0,m_pu
      real mpu
      real*8 mproton, eoverm, O_to_Ba
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)    !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)    !magnetic permeability of free space
      PARAMETER (epsilon0 = 8.85e-12) !dielectric constant

      PARAMETER (mproton = 1.67e-27)

      PARAMETER (km_to_m = 1e3)       !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)  !K

      real np_top,np_bottom
      real b0_top,b0_bottom,Lo,vth_top,vth_bottom,vth_max
      real m_top, m_bottom,m_heavy,np_bottom_proton


c electron ion collision frequency
      real nu_init,lww1,lww2

      PARAMETER (lww2 = 1.0)           !must be less than 1.0
      PARAMETER (lww1 = (1-lww2)/6.0)  !divide by six for nearest neighbor

c density scaling parameter, alpha, and ion particle array dims
       
      real*8 alpha, beta_pu  
      PARAMETER (beta_pu = 0.1)

      integer ri0 ! Number of grid cells pluto is offset from the center.

      real Qo, vrad, N_o, RIo, Rpluto, tau_photo, k_rec
      real S_radius
      PARAMETER (Qo = 3e27)       !neutral source rate
      PARAMETER (vrad = 0.05)     !escape velocity
c      PARAMETER (N_o = 5e34)     !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)      !Steady state neutral particle constant
      PARAMETER (RIo = 1200.0)    !Io radius
      PARAMETER (Rpluto = 1184.0) !Pluto radius
      PARAMETER (tau_photo = 1.5e9)
      PARAMETER (k_rec = 1e-5/1e15) !km^3 s^-1
      PARAMETER (S_radius = 30*Rpluto)

c domain decompostion parameters

      integer n_up, n_down, cart_dims
      integer io_proc

      parameter(cart_dims = 1) !dimensions for virtual topology
      parameter(n_up=1)
      parameter(n_down=2)

      logical periods(cart_dims), reorder
      data periods /.true./, reorder /.false./

      integer dims(cart_dims), tag, cart_coords(2)
      data tag /1/ !dimenions, /rows,columns/


c solar wind composition
      real b_sw_thermal_H, b_sw_thermal_He, b_sw_shell_H
      real f_sw_thermal_H, f_sw_thermal_He, f_sw_shell_H
      PARAMETER (b_sw_thermal_H = 1.0)
      PARAMETER (b_sw_thermal_He = 5.0)
      PARAMETER (b_sw_shell_H = 10.0)
      PARAMETER (f_sw_thermal_He = 0.08) 
      PARAMETER (f_sw_shell_H = 0.038) 
      PARAMETER (f_sw_thermal_H = 1.0 - f_sw_thermal_He - f_sw_shell_H) 
c      real f_mq_2,b_mq_2,f_shl,b_shl
c      PARAMETER (b_mq_2 = 2.0)
c      PARAMETER (f_mq_2 = 0.1*b_mq_2) 
c      PARAMETER (b_shl = 2.0)
c      PARAMETER (f_shl = 0.1*b_shl)

      real moment
      real surf_field
      real imf_theta
      real imf_phi

      ! PCE = PhotoChemical Equilibrium
      real ionization_frequency, recombination_rate, PCE_coef
      PARAMETER (ionization_frequency = 1.0e-9)!s^-1
      PARAMETER (recombination_rate = 1.0e-21) !km^3 s^-1
      PARAMETER (PCE_coef = 
     x       sqrt(ionization_frequency/recombination_rate)) !km^3 s^-1

      ! Inner radius of the optically thin atmosphere
      ! We globally cap ion density at the photochemical 
      ! equilibrium for this radius
      real r_thin
      PARAMETER (r_thin = 3.0*Rpluto)

      real max_ion_density

      real dummy_particle_tag
      real sw_thermal_H_tag
      real sw_thermal_He_tag
      real sw_shell_H_tag
      real pluto_photoionize_CH4_tag
      real pluto_stagnant_photoionize_CH4_tag
      real pluto_chex_CH4_tag
      PARAMETER (dummy_particle_tag = 0)
      PARAMETER (sw_thermal_H_tag = 1)
      PARAMETER (sw_thermal_He_tag = 2)
      PARAMETER (sw_shell_H_tag = 3)
      PARAMETER (pluto_photoionize_CH4_tag = 4)
      PARAMETER (pluto_stagnant_photoionize_CH4_tag = 5)
      PARAMETER (pluto_chex_CH4_tag = 6)
      CONTAINS

c----------------------------------------------------------------      
      subroutine readInputs()
c----------------------------------------------------------------      
      open(unit=100, file='inputs.dat', status='old')
      
      read(100,*) b0_init
      write(*,*) 'b0_init...........',b0_init
      read(100,*) ion_amu
      write(*,*) 'amu...............',ion_amu
      read(100,*) m_pu
      write(*,*) 'mpu...............',m_pu
      mpu = m_pu ! mpu is depreciated
      read(100,*) nf_init
      write(*,*) 'nf_init...........',nf_init
      read(100,*) dt_frac
      write(*,*) 'dt_frac...........',dt_frac
      read(100,*) nt
      write(*,*) 'nt................',nt
      read(100,*) nout
      write(*,*) 'nout..............',nout
      read(100,*) part_nout
      write(*,*) 'part_nout..............',part_nout
      read(100,*) output_wait
      write(*,*) 'output_wait..............',output_wait
      read(100,*) vsw
      write(*,*) 'vsw...............',vsw
      read(100,*) vth
      write(*,*) 'vth...............',vth
c      read(100,*) Ni_max
c      write(*,*) 'Ni_max....',Ni_max
      read(100,*) Ni_tot_frac
      write(*,*) 'Ni_tot_frac.......',Ni_tot_frac
      read(100,*) dx_frac
      write(*,*) 'dx_frac...........',dx_frac
      read(100,*) nu_init_frac
      write(*,*) 'nu_init_frac......',nu_init_frac
      read(100,*) out_dir
      write(*,*) 'output dir........',out_dir
      read(100,*) mrestart
      write(*,*) 'mrestart...........',mrestart
      read(100,*) ri0
      write(*,*) 'pluto offset.......',ri0
      read(100,*) surf_field
      write(*,*) 'surface field of pluto.......',surf_field
      read(100,*) imf_theta
      read(100,*) imf_phi
      write(*,*) 'IMF direction.......',imf_theta, imf_phi
     
      close(100)
      end subroutine readInputs
c----------------------------------------------------------------      

c----------------------------------------------------------------      
      subroutine check_inputs(my_rank)
c----------------------------------------------------------------      
      integer :: my_rank

      real*8 ak, btot, a1, a2, womega, phi, deltat

      !check input parameters

      if (my_rank .eq. 0) then 
      write(*,*) 'alpha...',alpha
      write(*,*) 'c/wpi...',lambda_i,dx,dy,delz
      write(*,*) 'dt......',dt,dtsub_init
      
      ak = 2./dx
      btot =  b0_init*q/mion
      a1 = ak**2*Btot/(alpha*(nf_init))
      a2 = (ak*Btot)**2/(alpha*(nf_init))
      womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
      phi = womega/ak
      deltat = dx/phi
      
      if (deltat/dtsub_init .le. 1) then 
         write(*,*) 'Field time step too long...'
         stop
      endif

      write(*,*) 'Courant check (>1)..',deltat/dtsub_init


      write(*,*) ' '
      write(*,*) 'Bottom parameters...'
      write(*,*) ' '
      va =  b0_init/sqrt(mu0*m_bottom*np_bottom/1e9)/1e3
      write(*,*) 'Alfven velocity.......',va
      write(*,*) 'Thermal velocity......',vth_top
      write(*,*) 'Mach number...........',vbottom/(va + vth_bottom)

      write(*,*) 'Thermal gyroradius..',m_bottom*vth_bottom/(q*b0_init),
     x            m_bottom*vth_bottom/(q*b0_init)/dx
      cwpi = 3e8/sqrt((np_bottom/1e9)*q*q/(epsilon0*m_bottom))
      write(*,*) 'Ion inertial length...',cwpi/1e3,cwpi/1e3/dx


      write(*,*) ' '
      write(*,*) 'Top parameters...'
      write(*,*) ' '

      va =  b0_init/sqrt(mu0*m_top*np_top/1e9)/1e3
      write(*,*) 'Alfven velocity.......',va
      write(*,*) 'Thermal velocity......',vth_top
      write(*,*) 'Mach number...........',vtop/(va + vth_top)

      write(*,*) 'Thermal gyroradius....',m_top*vth_top/(q*b0_init),
     x            m_top*vth_top/(q*b0_init)/dx
      cwpi = 3e8/sqrt((np_top/1e9)*q*q/(epsilon0*m_top))
      write(*,*) 'Ion inertial length...',cwpi/1e3,cwpi/1e3/dx


      endif


      end subroutine check_inputs
c----------------------------------------------------------------      

      
      END MODULE
