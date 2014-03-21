c para.h
c contains simulation parameter list for a barium release

c simulation domain dimensions
      PARAMETER (nx =121, ny =69, nz = 21)

c grid parameters
      PARAMETER (dx = 1800.0, dy = 1800.0)   !units in km
      PARAMETER (delz = 1800.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz

      PARAMETER (dx_buf = 15*dx)

c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
c      PARAMETER (dt = 0.05)     !main time step
      PARAMETER (nt = 500)        !number of time steps
c      PARAMETER (dtsub = 0.3) !subcycle time step 
      PARAMETER (dtsub_init = 0.2) !subcycle time step 
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (dt = dtsub_init*ntsub)     !main time step
c      PARAMETER (dtsub_init = 0.005) !subcycle time step 
c      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output 

c output directory
      character(50) out_dir
      PARAMETER (out_dir='/Volumes/MacD97-2/hybrid/3d_buf/run_test/')

c logical variable for restart
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 1600 )      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart


c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')

c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)

c neutral atmosphere
      real pwl,Ncol

      PARAMETER (pwl = 14)
      PARAMETER (Ncol = 6e16*1e10)  !km^-2

c neutral cloud expansion characteristics
      real vtop,vbottom,vth,vsw
      PARAMETER(vo = 20.0)
c      PARAMETER(vth = 10.0)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
      PARAMETER(vtop = -vsw)
      PARAMETER(vbottom = -vsw)

c total kinetic energy of released neutrals 
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2


c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 10000000)
      PARAMETER (Ni_max_buf = Ni_max/4)

c earths magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real*8 mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real*8 mproton,tempf0,m_pu,eoverm
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = mproton)    !mass of H (kg)
      PARAMETER (eoverm = q/mO)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (m_pu = 28.0)
      PARAMETER (mBa = m_pu*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K

c density scaling parameter, alpha, and ion particle array dims
       
      real*8 alpha, beta_pu

      PARAMETER (alpha = pi*4e-10*q*q/mO) !mH
      PARAMETER (beta_pu = 0.4)


      real Qo, vrad, N_o, RIo, tau_photo, k_rec 
      integer S_radius
	   PARAMETER (Qo = 2e27)  !neutral source rate
      PARAMETER (vrad = 0.1)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (RIo = 1200.0)  !Io radius
      PARAMETER (tau_photo = 1.2e9)
      PARAMETER (k_rec = 1e-6/1e15) !km^3 s^-1
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 30) !units of dx


c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef,np_top,np_bottom
      real b0_top,b0_bottom,Lo,vth_top,vth_bottom,vth_max
      real m_top, m_bottom
      real Mdot
      PARAMETER (Mdot = 0.0)
      PARAMETER (nf_init = 0.01e15)   !was 0.01
      PARAMETER (np_top = 1.0*nf_init)
      PARAMETER (np_bottom = nf_init)
      PARAMETER (b0_init = 0.2e-9)    !was 0.2
      PARAMETER (b0_top = 1.0*b0_init)
      PARAMETER (b0_bottom = 1.0*b0_init)
      PARAMETER (vth_top = 19.3)
      PARAMETER (vth_bottom = 19.3)
      PARAMETER (vth_max = 3*19.3)
      PARAMETER (m_top = mproton)
      PARAMETER (m_bottom = mproton)
      PARAMETER (nn_coef = 1e5)
      PARAMETER (Lo = 1.5*dx) !gradient scale length of boundary


c domain decompostion parameters

      integer n_up,n_down,n_left,n_right, cart_dims,comm_sz
      integer io_proc

      parameter(comm_sz = 10)
      parameter(io_proc = nint(comm_sz/2.))
      parameter(cart_dims = 2) !dimensions for virtual topology
      parameter(n_up=1)
      parameter(n_down=2)
      parameter(n_left=3)
      parameter(n_right=4)

      logical periods(cart_dims), reorder
      data periods /.true.,.false./, reorder /.false./

      integer dims(cart_dims), tag, cart_coords(2)
      data dims /comm_sz,1/, tag /1/ !dimenions, /rows,columns/



c electron ion collision frequency
      real nu_init, eta_init,lww1,lww2
      PARAMETER (nu_init = 0.1*q*b0_init/mproton)
      PARAMETER (eta_init = 0.0)
      PARAMETER (lww2 = 1.0)    !must be less than 1.0
      PARAMETER (lww1 = (1-lww2)/6.0)  !divide by six for nearest neighbor


c solar wind composition
      real f_mq_2,b_mq_2,f_shl,b_shl
      PARAMETER (f_mq_2 = 0.1)
      PARAMETER (b_mq_2 = 2.0) 
      PARAMETER (f_shl = 0.1)
      PARAMETER (b_shl = 2.0)








