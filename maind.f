      PROGRAM MAIND
     
c----------------------------------------------------------------------
c maind.f
c Parallel version with no ion fluid, Nov 24, 2004
c----------------------------------------------------------------------

      
      USE global
      USE dimensions
      USE inputs
      USE mpi
      USE initial
      USE misc
      USE gutsp_dd
      USE gutsp_buf
      USE gutsf
      USE part_init
      USE grid_interp
      USE chem_rates
      USE iso_fortran_env, only: output_unit,error_unit
      USE ifport, only: rename


c----------------------------------------------------------------------
c Listing of all declared variables
c
c Note that the format for specifying variables at different time
c levels is based upon the number of 1/2 time steps that the varible
c is behind the current time level.  For example, 
c uf2 is 1 time step behind uf, ufp2 is 1 time step ahead of uf,
c and b12 is 1 full time step behind b1 (not 12 1/2 steps behind b 
c since b does not exist...right). b1p2 is an exception...this is a
c temporary holder for b1 at m+1 in the predictor/corrector update
c of the magnetic field.
c----------------------------------------------------------------------
      save


      integer dumb(1)
      real b0(nx,ny,nz,3),            !ambient magnetic field
     x     b1(nx,ny,nz,3),    !1st order magnetic field
     x     b12(nx,ny,nz,3),   !b1 at previous time step
     x     b1p2(nx,ny,nz,3),  !temporary b1 at time level m+1
     x     bt(nx,ny,nz,3),    !total magnetic field..mc covarient
     x     btc(nx,ny,nz,3),   !btmf at cell center for particle move
     x     np(nx,ny,nz),      !ion number density at time 
                              !level n, n+1/2 including dummies
     x     np_tot(nx,ny,nz),  !ion number density at time 
                              !level n, n+1/2
     x     np_H(nx,ny,nz),    !proton number density at time 
                              !level n, n+1/2
     x     np_He(nx,ny,nz),   !He++ number density at time 
                              !level n, n+1/2
     x     np_H_shell(nx,ny,nz),   !shell number density at time 
                                 !level n, n+1/2
     x     np_He_shell(nx,ny,nz),   !shell number density at time 
                                 !level n, n+1/2
     x     np_sw(nx,ny,nz),   !solar wind number density at time 
                              !level n, n+1/2
     x     np_CH4(nx,ny,nz),  !CH4+ number density at time 
                              !level n, n+1/2
     x     np_dummy(nx,ny,nz),  !CH4+ number density at time 
                              !level n, n+1/2
     x     vp(Ni_max,3),      !particle velocity at t level n+1/2
     x     vp1(Ni_max,3),     !particle velocity at t level n
     x     vplus(Ni_max,3),   !v+ used in velocity update
     x     vminus(Ni_max,3),  !v- used in velocity update
     x     up(nx,ny,nz,3),    !particle flow at time level n, n+1/2
     x     up_tot(nx,ny,nz,3),
     x     up_H(nx,ny,nz,3),
     x     up_He(nx,ny,nz,3),
     x     up_H_shell(nx,ny,nz,3),
     x     up_He_shell(nx,ny,nz,3),
     x     up_sw(nx,ny,nz,3),
     x     up_CH4(nx,ny,nz,3),
     x     up_dummy(nx,ny,nz,3),  

     x     xp(Ni_max,3),      !coordinates of ion particles
     x     aj(nx,ny,nz,3),    !curlB/(alpha*n) 
     x     nu(nx,ny,nz),      !collision frequency
     x     nu_background(nx,ny,nz),
     x     Ep(Ni_max,3),      !Ion particle electric field
     x     E(nx,ny,nz,3)     !electric field from electron mom eqn


      real temp_p(nx,ny,nz)
      real temp_tot(nx,ny,nz)
      real temp_h(nx,ny,nz)
      real temp_he(nx,ny,nz)
      real temp_shell(nx,ny,nz)
      real temp_sw(nx,ny,nz)
      real temp_ch4(nx,ny,nz)

      real Evp,       !total particle kinetic energy
     x     EB1,       !total magnetic field energy
     x     EB1x,      !total b1x energy
     x     EB1y,      !total b1y energy
     x     EB1z,      !total b1z energy
     x     EE,        !total electric field energy
     x     EeP        !total electron pressure energy

      real pup(3),      !total particle momentum
     x     peb(3)       !total momentum carried by E and B fields

      real mr

      real chex_rate
      real bill_rate
      real satnp
      real mindt
      integer t1,t2,cnt_rt
      real time
      integer ierr

      integer ndiag
      integer ndiag_part
      integer n

      real ndot(nx,ny,nz)
      
      integer seedsize
      integer, dimension(:), allocatable :: cursteps

      real recvbuf
      integer count
      integer restart_counter

      integer i,j,k,l,m,mstart
      
c      character filenum
      character flnm
      character(len=:), allocatable::filenum
      character(len=10) :: arg
      character(len=10) :: acc
      character(len=7) :: stat

      ! extra clock stuff
      integer(kind=8) :: sys_count, sys_rate, sys_max
      real :: cpu_start, cpu_end

      logical ex
      integer para_dat_version
      para_dat_version = 4

c----------------------------------------------------------------------

      call readInputs()
      write(*,*) 'initializing variables...'
      call initparameters()

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

      ! Do some initial setup of the working directory.
      if(my_rank .eq. 0) then
          call execute_command_line('mkdir -p '//trim(out_dir)//'grid')
          call execute_command_line(
     x           'mkdir -p '//trim(out_dir)//'particle')
          call execute_command_line(
     x           'cp --backup=numbered inputs.dat '//trim(out_dir))
      endif

      filenum = int_to_str(my_rank+1)

      write(*,*) 'filenum...',filenum,my_rank+1

      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      allocate( cursteps(procnum) )
      io_proc = nint(procnum/2.)
      dims(1) = procnum

c create virtual topology (set dimensions in para.h)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_CART_CREATE(MPI_COMM_WORLD, cart_dims, dims, periods, 
     x     reorder,cartcomm, ierr)

      call MPI_COMM_RANK(cartcomm, cart_rank, ierr)
      call MPI_CART_COORDS(cartcomm, cart_rank, cart_dims, cart_coords, 
     x                     ierr)
      call MPI_CART_SHIFT(cartcomm,0,1,up_proc,down_proc,ierr)

      call system_clock(t1,cnt_rt)


      call get_command_argument(number=1,value=arg,status=ierr)
      restart = (trim(arg) == "restart")

      print *,'Ni_tot_0, Ni_tot',Ni_tot_0,Ni_tot
      
      if (my_rank .eq. 0) then
         call check_inputs(my_rank)
         write(*,*) 'Total Particles per cell....',Ni_tot/num_cells
         write(*,*) ' '
      endif
         

      ndiag = nout-1
      ndiag_part = part_nout-1
      nuei = 0.0

c initialize seed for each processor

      seed = t1 +my_rank*100

      ! Make sure all ranks are initialized and sychronized 
      ! before attempting to make a system call
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
      
      call get_command_argument(number=1,value=arg,status=ierr)
      if(ierr .eq. 0 .and. trim(arg) .eq. "debug") then
            call debug_random_initialize()
      else 
            call random_initialize()
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

      input_E = 0.0
      input_chex = 0.0
      input_bill = 0.0

      call grd_no_strech()
      call get_nu(nu_background)
      nu = nu_background
      call grd6_setup(b0,b1,b1p2,bt)
      call get_beta()


      if (.not.(restart)) then
         call cpu_time(cpu_start)
         call sw_part_setup(np,vp,vp1,xp,up)
         call cpu_time(cpu_end)
         write(*,*) "Particle setup time = ",cpu_end-cpu_start,
     x              "seconds"

         call f_update_tlev(b1,b12,b1p2,bt,b0)
      endif

      call initialize_buffer()
c----------------------------------------------------------------------




c----------------------------------------------------------------------
c check for restart flag
c----------------------------------------------------------------------
      inquire (file=trim(out_dir)//'para.dat',exist=ex)
      
      ! sanity check to ensure no files get overwritten
      if(restart .and. (.not. ex)) then
          write(*,*) 'restart is true, but '//trim(out_dir)//'/para.dat
     x      does not exist'
          write(*,*) 'stopping'
          call MPI_FINALIZE(ierr)
          stop
      else if((.not. restart) .and. ex) then
          write(*,*) 'not a restart, but '//trim(out_dir)//'/para.dat
     x        exists'
          write(*,*) 'stopping'
          call MPI_FINALIZE(ierr)
          stop
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      ! Read restart files
      mstart = 0
      write(*,*) 'restart status....',restart
      if (restart) then 
         write(*,*) 'opening restart.vars......'


          open(1000+my_rank,file=trim(out_dir)//'restart.vars'//filenum,
     x          status='unknown',
     x          form='unformatted')
          write(*,*) 'reading restart.vars......',filenum
         
          read(1000+my_rank) b1,b12,
     x         E,input_E,mstart,
     x         Evp,EB1,EE, Ni_tot, Ni_tot_0

          bt = b0 + b1

          close(1000+my_rank)
          open(1000+my_rank,file=trim(out_dir)//'restart.part'//filenum,
     x         status='unknown',form='unformatted')
          write(*,*) 'reading restart.part......',filenum
          read(1000+my_rank) vp(:Ni_tot,:),vp1(:Ni_tot,:),
     x         xp(:Ni_tot,:),
     x         mrat(:Ni_tot), tags(:Ni_tot), beta_p(:Ni_tot)

          close(1000+my_rank)

      endif
      restart_counter = mstart + mrestart



c----------------------------------------------------------------------
c write para.h file

      if (my_rank .eq. 0) then


         open(109, file=trim(out_dir)//'para.dat',
     x        status='unknown',form='unformatted')
         
         write(109) para_dat_version
         write(109) nx,ny,nz,dx,dy,delz
         write(109) nt,dtsub_init,ntsub,dt,nout
         write(109) out_dir
         write(109) vtop,vbottom
         write(109) Ni_max
         write(109) mproton,m_pu,m_heavy
         write(109) np_top,np_bottom
         write(109) b0_top,b0_bottom
         write(109) vth_top,vth_bottom
         write(109) alpha,beta
         write(109) RIo

         write(109) b0_init
         write(109) ion_amu
         write(109) m_pu
         write(109) nf_init
         write(109) dt_frac
         write(109) vsw
         write(109) vth
         write(109) Ni_tot_frac
         write(109) dx_frac
         write(109) nu_init_frac
         write(109) mrestart
         write(109) ri0

         write(109) part_nout

         close(109)

      endif
 
c----------------------------------------------------------------------



c----------------------------------------------------------------------
c Initialize diagnostic output files
c----------------------------------------------------------------------

      if(restart) then
          acc = 'append'
          stat = 'unknown'
      else
          acc = 'sequential'
          stat = 'new'
      endif


      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      open(110,file=trim(out_dir)//'grid/'//
     x     'c.np_'//filenum//'.dat', access= acc,
     x     status='unknown',form='unformatted')
      open(111,file=trim(out_dir)//'grid/'//
     x     'c.np_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')

      open(115,file=trim(out_dir)//'grid/'//
     x     'c.np_H_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(116,file=trim(out_dir)//'grid/'//
     x     'c.np_He_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(117,file=trim(out_dir)//'grid/'//
     x     'c.np_H_shell_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(1175,file=trim(out_dir)//'grid/'//
     x     'c.np_He_shell_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(118,file=trim(out_dir)//'grid/'//
     x     'c.np_sw_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(119,file=trim(out_dir)//'grid/'//
     x     'c.np_CH4_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(120,file=trim(out_dir)//'grid/'//
     x     'c.np_tot_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(121,file=trim(out_dir)//'grid/'//
     x     'c.np_dummy_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')

      open(9115,file=trim(out_dir)//'grid/'//
     x     'c.up_H_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(9116,file=trim(out_dir)//'grid/'//
     x     'c.up_He_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(9117,file=trim(out_dir)//'grid/'//
     x     'c.up_H_shell_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(91175,file=trim(out_dir)//'grid/'//
     x     'c.up_He_shell_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(9118,file=trim(out_dir)//'grid/'//
     x     'c.up_sw_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(9119,file=trim(out_dir)//'grid/'//
     x     'c.up_CH4_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(9120,file=trim(out_dir)//'grid/'//
     x     'c.up_tot_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(9121,file=trim(out_dir)//'grid/'//
     x     'c.up_dummy_3d_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')

      open(130,file=trim(out_dir)//'grid/'//
     x     'c.b1_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(131,file=trim(out_dir)//'grid/'//
     x     'c.b1_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(132,file=trim(out_dir)//'grid/'//
     x     'c.b0_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(133,file=trim(out_dir)//'grid/'//
     x     'c.bt_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(134,file=trim(out_dir)//'grid/'//
     x     'c.aj_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(140,file=trim(out_dir)//'grid/'//
     x     'c.aj_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(150,file=trim(out_dir)//'grid/'//
     x     'c.E_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(180,file=trim(out_dir)//'grid/'//
     x     'c.up_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(181,file=trim(out_dir)//'grid/'//
     x     'c.up_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(300,file=trim(out_dir)//'grid/'//
     x     'c.temp_p_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(301,file=trim(out_dir)//'grid/'//
     x     'c.temp_tot_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(302,file=trim(out_dir)//'grid/'//
     x     'c.temp_h_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(303,file=trim(out_dir)//'grid/'//
     x     'c.temp_he_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(304,file=trim(out_dir)//'grid/'//
     x     'c.temp_shell_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(305,file=trim(out_dir)//'grid/'//
     x     'c.temp_sw_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(306,file=trim(out_dir)//'grid/'//
     x     'c.temp_ch4_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(307,file=trim(out_dir)//'particle/'//
     x     'c.xp_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(310,file=trim(out_dir)//'particle/'//
     x     'c.vp_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(311,file=trim(out_dir)//'particle/'//
     x     'c.tags_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(315,file=trim(out_dir)//'particle/'//
     x     'c.beta_p_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(320,file=trim(out_dir)//'particle/'//
     x     'c.mrat_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(9000,file=trim(out_dir)//'particle/'//
     x     'c.outflowing_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(9900,file=trim(out_dir)//
     x     'c.time_stepping_error_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

c----------------------------------------------------------------------


c======================================================================
c  MAIN LOOP!
c======================================================================

      do 1 m = mstart+1, nt

         dumb(1) = m
         call MPI_GATHER(dumb, 1, MPI_INT, cursteps, 1, MPI_INT, 0,
     x                MPI_COMM_WORLD, ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         if (my_rank .eq. 0) then
         do i=1,procnum
         if (cursteps(i) .ne. m) then
             call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr)
         endif
         enddo
         write(*,*) 'All the same steps'
         call flush(output_unit)
         endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


         ndiag = ndiag + 1
         ndiag_part = ndiag_part + 1

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         if (my_rank .eq. 0) then
            write(*,*) 'time...', m, dt,mstart
            call flush(output_unit)
         endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         write(*,*) 'Ni_tot...',Ni_tot,Ni_max,my_rank,m

         call get_interp_weights(xp)
         call update_np(xp, vp, vp1, np)             !np at n+1/2
         call update_np_boundary(np)
         call update_up(vp,np,up)       !up at n+1/2

         call ionization(np,xp,vp,vp1)


         call curlB(b1,np,aj)
         call edge_to_center(bt,btc,b0_us)

         call extrapol_up(up,vp,vp1,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call improve_up(vp1,vplus,vminus,up,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call get_vp_final(Ep,vp,vp1,vplus)
         
         call move_ion_half(xp,vp,vp1,Ep)

         call part_setup_buf(xp_buf,vp_buf)
         call get_Ep_buf(Ep_buf,b0,xp_buf,up)
         call get_vplus_vminus_buf(Ep_buf,vp_buf,vplus_buf,
     x        vminus_buf,b0)
         call get_vp_buf_final(Ep_buf,vp_buf,vplus_buf)
         call move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
         call exchange_ion_in(xp,vp,vp1,xp_buf,vp_buf)
         call exchange_ion_out(xp,vp,vp1,xp_buf,vp_buf,
     x        E,Bt,9000)
!         call reflect_boundary(xp,vp,vp1)

         call exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)



         call get_interp_weights(xp)
         call update_np(xp, vp, vp1, np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2
         call update_np_boundary(np)
         ! These variables do not impact the simulation and only need
         ! to be updated on the output step
         if (ndiag .eq. nout) then
               call get_temperature(xp,vp,np,temp_p)

               ! total
               call separate_temperature(xp,vp,np,temp_tot,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag 
     x           .or. tags(:Ni_tot)==sw_thermal_He_tag
     x           .or. tags(:Ni_tot)==sw_shell_H_tag
     x           .or. tags(:Ni_tot)==pluto_photoionize_CH4_tag
     x           .or. tags(:Ni_tot)==pluto_stagnant_photoionize_CH4_tag
     x           .or. tags(:Ni_tot)==pluto_chex_CH4_tag))
               ! Thermal H
               call separate_temperature(xp,vp,np,temp_H,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag ))
               ! He
               call separate_temperature(xp,vp,np,temp_He,
     x                   (tags(:Ni_tot)==sw_thermal_He_tag ))
               ! Shell H
               call separate_temperature(xp,vp,np,temp_shell,
     x                   (tags(:Ni_tot)==sw_shell_H_tag ))
               ! Solar wind
               call separate_temperature(xp,vp,np,temp_sw,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag 
     x               .or. tags(:Ni_tot)==sw_thermal_He_tag
     x               .or. tags(:Ni_tot)==sw_shell_H_tag))
               ! CH4
               call separate_temperature(xp,vp,np,temp_CH4,
     x                   (tags(:Ni_tot)==pluto_photoionize_CH4_tag
     x          .or. tags(:Ni_tot)==pluto_stagnant_photoionize_CH4_tag
     x          .or. tags(:Ni_tot)==pluto_chex_CH4_tag))

               ! total
               call separate_np(np_tot,
     x                   (tags(:Ni_tot) .ne. dummy_particle_tag ))
               call separate_np(np_H,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag ))
               call separate_np(np_He,
     x                   (tags(:Ni_tot)==sw_thermal_He_tag))
               call separate_np(np_H_shell,
     x                   (tags(:Ni_tot)==sw_shell_H_tag))
               call separate_np(np_He_shell,
     x                   (tags(:Ni_tot)==sw_shell_He_tag))
               call separate_np(np_sw,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag 
     x               .or. tags(:Ni_tot)==sw_thermal_He_tag
     x               .or. tags(:Ni_tot)==sw_shell_H_tag
     x               .or. tags(:Ni_tot)==sw_shell_He_tag))
               call separate_np(np_CH4,
     x                   (tags(:Ni_tot)==pluto_photoionize_CH4_tag
     x          .or. tags(:Ni_tot)==pluto_stagnant_photoionize_CH4_tag
     x          .or. tags(:Ni_tot)==pluto_chex_CH4_tag))
               call separate_np(np_dummy,
     x                   (tags(:Ni_tot)==dummy_particle_tag ))
               call separate_up(vp,up_tot,
     x                   (tags(:Ni_tot) .ne. dummy_particle_tag ))
               call separate_up(vp,up_H,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag ))
               call separate_up(vp,up_He,
     x                   (tags(:Ni_tot)==sw_thermal_He_tag))
               call separate_up(vp,up_H_shell,
     x                   (tags(:Ni_tot)==sw_shell_H_tag))
               call separate_up(vp,up_He_shell,
     x                   (tags(:Ni_tot)==sw_shell_He_tag))
               call separate_up(vp,up_sw,
     x                   (tags(:Ni_tot)==sw_thermal_H_tag 
     x               .or. tags(:Ni_tot)==sw_thermal_He_tag
     x               .or. tags(:Ni_tot)==sw_shell_H_tag
     x               .or. tags(:Ni_tot)==sw_shell_He_tag))
               call separate_up(vp,up_CH4,
     x                   (tags(:Ni_tot)==pluto_photoionize_CH4_tag
     x          .or. tags(:Ni_tot)==pluto_stagnant_photoionize_CH4_tag
     x          .or. tags(:Ni_tot)==pluto_chex_CH4_tag))
               call separate_up(vp,up_dummy,
     x                   (tags(:Ni_tot)==dummy_particle_tag ))

         endif
         
c**********************************************************************
c SUBCYCLING LOOP!
c**********************************************************************

         ! Set subcycle timestep and number of steps back to their
         ! fixed initial values (dtsub_init and ntsub don't change)
         dtsub = dtsub_init
         ntf = ntsub
         
         call check_time_step(b0,b1,bt,np,m,9900)

         count = 1

         call MPI_ALLREDUCE(ntf,recvbuf,count,
     x        MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)

         ntf = recvbuf

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do 2 n = 1, int(ntf)

         !convert main cell covarient bt to main cell contravarient
         call curlB(b1,np,aj)

         ! update nu here
         !call update_nu(nu, nu_background, aj, bt)

         call predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu) 
         call correct_B(b0,b1,b1p2,E,aj,up,np,nu)


         call f_update_tlev(b1,b12,b1p2,bt,b0)


 2     continue
c**********************************************************************


         call move_ion_half(xp,vp,vp1,Ep)

         call part_setup_buf(xp_buf,vp_buf)
         call move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
         
         call exchange_ion_in(xp,vp,vp1,xp_buf,vp_buf)
         call exchange_ion_out(xp,vp,vp1,xp_buf,vp_buf,
     x        E,Bt,9000)
!         call reflect_boundary(xp,vp,vp1)

         call exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)


         call check_min_den(np,xp,vp,vp1,up,bt)




c----------------------------------------------------------------------
c diagnostic output
c----------------------------------------------------------------------

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


         !energy diagnostics
         call get_bndry_Eflux(b1,E)
         call Energy_diag(vp,b0,b1,E,Evp,EB1,EE,
     x                    nu,up,np)

         !if ((ndiag .ge. nout) .and. (m .ge. output_wait)) then
         if ((ndiag .ge. nout) .and. (m .ge. output_wait)) then

c save 3d arrays------------------------
               ! Output grid data for the whole domain
               write(115) m
               write(115) np_H
               write(116) m
               write(116) np_He
               write(117) m
               write(117) np_H_shell
               write(1175) m
               write(1175) np_He_shell
               write(118) m
               write(118) np_sw
               write(119) m
               write(119) np_CH4
               write(120) m
               write(120) np_tot
               write(121) m
               write(121) np_dummy

               write(9115) m
               write(9115) up_H
               write(9116) m
               write(9116) up_He
               write(9117) m
               write(9117) up_H_shell
               write(91175) m
               write(91175) up_He_shell
               write(9118) m
               write(9118) up_sw
               write(9119) m
               write(9119) up_CH4
               write(9120) m
               write(9120) up_tot
               write(9121) m
               write(9121) up_dummy

               write(133) m
               write(133) bt
               write(134) m
               write(134) aj
               write(150) m
               write(150) E
               write(300) m
               write(300) temp_p/1.6e-19
               write(301) m
               write(301) temp_tot/1.6e-19
               write(302) m
               write(302) temp_h/1.6e-19
               write(303) m
               write(303) temp_he/1.6e-19
               write(304) m
               write(304) temp_shell/1.6e-19
               write(305) m
               write(305) temp_sw/1.6e-19
               write(306) m
               write(306) temp_ch4/1.6e-19

               ndiag = 0

         endif
               ! Only output particle data only near pluto
         if ( m .eq. nt-1 ) then
               if ( my_rank .gt. procnum/2 - 30 .and.
     x              my_rank .lt. procnum/2 + 30) then
                   write(307) m
                   write(307) ((xp(i,j), i=1,Ni_tot),j=1,3)
                   write(310) m
                   write(310) ((vp(i,j), i=1,Ni_tot),j=1,3)
                   write(311) m
                   write(311) tags(:Ni_tot)
                   write(315) m
                   write(315) beta_p(:Ni_tot)
                   write(320) m
                   write(320) mrat(:Ni_tot)
                   ndiag_part = 0
               endif
               endif

c----------------------------------------------------------------------


c----------------------------------------------------------------------
c Write restart file
c----------------------------------------------------------------------

         if (m .eq. restart_counter) then
          ierr = rename(trim(out_dir)//'restart.vars'//filenum, 
     x           trim(out_dir)//'restart.vars'//filenum//'.old')
          ierr = rename(trim(out_dir)//'restart.part'//filenum, 
     x           trim(out_dir)//'restart.part'//filenum//'.old')


          write(*,*) 'writing restart file....',
     x                 'restart.part'//filenum//'.new',my_rank,cart_rank
          open(1000+my_rank,file=trim(out_dir)//'restart.vars'//filenum,
     x              status='unknown',
     x              form='unformatted')
         
          write(1000+my_rank)  b1,b12,
     x             E,input_E,m,
     x             Evp,EB1,EE, Ni_tot, Ni_tot_0

          close(1000+my_rank)
          open(1000+my_rank,file=trim(out_dir)//'restart.part'//filenum,
     x             status='unknown',form='unformatted')
          write(1000+my_rank) vp(:Ni_tot,:),vp1(:Ni_tot,:),
     x             xp(:Ni_tot,:),
     x             mrat(:Ni_tot), tags(:Ni_tot), beta_p(:Ni_tot)

          close(1000+my_rank)

          ! Barrier to make sure all processors are done writing.
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ! Delete the .old files.
          open(1000+my_rank,
     x      file=trim(out_dir)//'restart.vars'//filenum//'.old',
     x      status='unknown',form='unformatted')
          close(1000+my_rank, status='delete')
          open(1000+my_rank,
     x      file=trim(out_dir)//'restart.part'//filenum//'.old',
     x      status='unknown',form='unformatted')
          close(1000+my_rank, status='delete')
                  
          restart_counter = restart_counter + mrestart


         endif

c----------------------------------------------------------------------




 1     continue
c======================================================================


          close(110)
          close(115)
          close(116)
          close(117)
          close(118)
          close(120)
          close(130)
          close(140)
          close(150)
          close(160)
          close(170)
          close(172)
          close(175)
          close(180)
          close(190)
          close(192)
          close(210)
          close(211)
          close(220)
          close(221)
          close(300)
          close(301)
          close(302)
          close(303)
          close(304)
          close(305)
          close(310)
          close(315)
          close(320)
          close(350)


          deallocate(filenum)

       call system_clock(t2,cnt_rt)
       time = (real(t2) - real(t1))/real(cnt_rt)
       if (my_rank .eq. 0) then
          write(*,*) 
          write(*,*) 'Elapsed time....',time,' sec'
          write(*,*)
       endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


       call MPI_FINALIZE(ierr)

       stop 
       end
