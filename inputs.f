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

c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_tot_0


c misc constants
      real*8 mu0,epsilon,pi,rtod,mO,mBa,km_to_m,kboltz,melec
      real*8 tempf0,m_pu
      real*8 mproton, eoverm, O_to_Ba
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)    !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)    !magnetic permeability of free space
      PARAMETER (epsilon0 = 8.85e-12) !dielectric constant

c      PARAMETER (m_pu = 64.0)
      PARAMETER (mproton = 1.67e-27)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
c      PARAMETER (mBa = m_pu*mO)    !mass of Ba (kg)
c      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move

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


      real Qo, vrad, N_o, RIo, Rpluto, tau_photo, k_rec, ri0
      real S_radius
      PARAMETER (Qo = 3e27)       !neutral source rate
      PARAMETER (vrad = 0.05)     !escape velocity
c      PARAMETER (N_o = 5e34)     !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)      !Steady state neutral particle constant
      PARAMETER (RIo = 1200.0)    !Io radius
      PARAMETER (Rpluto = 1184.0) !Pluto radius
      PARAMETER (tau_photo = 1.5e9)
      PARAMETER (k_rec = 1e-5/1e15) !km^3 s^-1
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 200*Rpluto)

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
      real f_mq_2,b_mq_2,f_shl,b_shl
      PARAMETER (b_mq_2 = 2.0)
      PARAMETER (f_mq_2 = 0.1*b_mq_2) 
      PARAMETER (b_shl = 2.0)
      PARAMETER (f_shl = 0.1*b_shl)

      real moment
      real surf_field
      real imf_theta
      real imf_phi
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
      subroutine initparameters()
c----------------------------------------------------------------      

      mion = ion_amu*1.67e-27

      lambda_i = (3e8/
     x            sqrt((nf_init/1e9)*q*q/(8.85e-12*mion)))/1e3

      dx = lambda_i*dx_frac 
      dy = lambda_i*dx_frac   !units in km
      delz = lambda_i*dx_frac          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz


      print*, 'dx....',dx

      dx_buf = 11*dx


      dt = dt_frac*mion/(q*b0_init)     !main time step
      dtsub_init = dt/ntsub !subcycle time step 

      print*, 'dt...',dt,dtsub_init

      vtop = vsw
      vbottom = -vsw

      Ni_tot_0 = Ni_max*Ni_tot_frac

      write(*,*) 'Ni_tot_0...',Ni_tot_0, Ni_max,Ni_tot_frac

      mO = mion    !mass of H (kg)

      mBa = m_pu*mO    !mass of Ba (kg)

      eoverm = q/mO

      m_heavy = 1.0
      np_top = nf_init
      np_bottom = nf_init/m_heavy
      f_proton_top = 0.5       !fraction relative to top
      b0_top = 1.0*b0_init
      b0_bottom = b0_init
      vth_top = vth
      vth_bottom = vth
      vth_max = 3*vth
      m_top = mion
      m_bottom = mion
      Lo = 4.0*dx             !gradient scale length of boundary

      nu_init = nu_init_frac*q*b0_init/mion

      alpha = (mu0/1e3)*q*(q/mion) !mH...determines particle scaling

      moment = surf_field * Rpluto**3

      imf_theta = (pi/180)*imf_theta
      imf_phi = (pi/180)*imf_phi

      end subroutine initparameters
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

c      write(*,*) 'Particles per cell....',Ni_tot_sys/(nx*nz)

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
