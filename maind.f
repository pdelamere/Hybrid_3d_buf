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
      USE iso_fortran_env, only: error_unit


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


      real b0(nx,ny,nz,3),            !ambient magnetic field
     x     b1(nx,ny,nz,3),    !1st order magnetic field
     x     b12(nx,ny,nz,3),   !b1 at previous time step
     x     b1p2(nx,ny,nz,3),  !temporary b1 at time level m+1
     x     bt(nx,ny,nz,3),    !total magnetic field..mc covarient
     x     btc(nx,ny,nz,3),   !btmf at cell center for particle move
     x     np(nx,ny,nz),      !particle ion den at time level n, n+1/2
     x     np_1(nx,ny,nz),
     x     np_2(nx,ny,nz),
     x     vp(Ni_max,3),      !particle velocity at t level n+1/2
     x     vp1(Ni_max,3),     !particle velocity at t level n
     x     vplus(Ni_max,3),   !v+ used in velocity update
     x     vminus(Ni_max,3),  !v- used in velocity update
     x     up(nx,ny,nz,3),    !particle flow at time level n, n+1/2
     x     xp(Ni_max,3),      !coordinates of ion particles
     x     aj(nx,ny,nz,3),    !curlB/(alpha*n) 
     x     nu(nx,ny,nz),      !collision frequency
     x     Ep(Ni_max,3),      !Ion particle electric field
     x     E(nx,ny,nz,3)     !electric field from electron mom eqn


      real xp_buf(Ni_max_buf,3)
      real vp_buf(Ni_max_buf,3)
      real Ep_buf(Ni_max_buf,3)
      real vplus_buf(Ni_max_buf,3)
      real vminus_buf(Ni_max_buf,3)

      real xp_out_buf(Ni_max_buf,3)
      real vp_out_buf(Ni_max_buf,3)
      real E_out_buf(Ni_max_buf,3)
      real B_out_buf(Ni_max_buf,3)
      real mrat_out_buf(Ni_max_buf)

      real temp_p(nx,ny,nz)

      real Evp,       !total particle kinetic energy
     x     Euf,       !total fluid kinetic energy
     x     EB1,       !total magnetic field energy
     x     EB1x,      !total b1x energy
     x     EB1y,      !total b1y energy
     x     EB1z,      !total b1z energy
     x     EE,        !total electric field energy
     x     EeP        !total electron pressure energy

      real pup(3),      !total particle momentum
     x     puf(3),      !total fluid momentum
     x     peb(3),      !total momentum carried by E and B fields
     x     input_p(3)   !input momentum

      real mr

      real chex_rate
      real bill_rate
      real satnp
      real mindt
      integer*4 t1,t2,cnt_rt
      real time
      integer ierr

      real ndot(nx,ny,nz)
      
      integer seedsize
      integer, dimension(:), allocatable :: seeder

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

      logical ex
      integer para_dat_version = 3




c----------------------------------------------------------------------

      call readInputs()
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

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      write(*,*) 'filenum...',filenum,my_rank+1

      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
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

c----------------------------------------------------------------------
c Initialize all variables
c----------------------------------------------------------------------
      write(*,*) 'initializing variables...'


      call get_command_argument(number=1,value=arg,status=ierr)
      restart = (trim(arg) == "restart")

      if (.not. restart) then
         Ni_tot = Ni_tot_0
         Ni_tot_sw = Ni_tot
         Ni_tot_sys = Ni_tot
         print *,'Ni_tot_sys, Ni_tot..',Ni_tot_sys,Ni_tot,Ni_tot_sw
      endif
      
      if (my_rank .eq. 0) then
         call check_inputs(my_rank)
         write(*,*) 'Particles per cell....',Ni_tot_sys/(nx*ny*nz)
         write(*,*) ' '
      endif
         

      ndiag = 0
      ndiag_part = 0
      prev_Etot = 1.0
      nuei = 0.0

c initialize seed for each processor

      seed = t1 +my_rank*100

      ! Make sure all ranks are initialized and sychronized before attempting
      ! to make a system call
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
      
      call get_command_argument(number=1,value=arg,status=ierr)
      if(ierr .eq. 0 .and. trim(arg) .eq. "debug") then
            call debug_random_initialize()
      else 
            call random_initialize()
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

      if (.not.(restart)) then
         do 66 i=1,nx
            do 66 j=1,ny
               do 66 k=1,nz
                  input_E = 0.0
                  input_p = 0.0
                  input_chex = 0.0
                  input_bill = 0.0
 66               continue
               endif


      if (.not.(restart)) then
         mrat(1:Ni_tot) = 1.0
         mrat(Ni_tot+1:) = 1.0/m_pu !mass N_2+ = 28.0
         beta_p(1:Ni_tot) = 1.0
         beta_p(Ni_tot+1:) = beta_pu
      endif

      call grd8()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu)


      if (.not.(restart)) then
         call get_beta()
      endif


      if (.not.(restart)) then
         call sw_part_setup_maxwl(np,vp,vp1,xp,input_p,up)

         call part_setup_buf(xp_buf,vp_buf)
         
         call part_setup_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf,b0)
                  

         call f_update_tlev(b1,b12,b1p2,bt,b0)
      endif

c----------------------------------------------------------------------




c----------------------------------------------------------------------
c check for restart flag
c----------------------------------------------------------------------
      inquire (file=trim(out_dir)//'para.dat',exist=ex)
      
      ! sanity check to ensure no files get overwritten
      if(restart .and. (.not. ex)) then
          write(*,*) 'restart is true, but '//trim(out_dir)//' does not
     x      exist'
          write(*,*) 'stopping'
          call MPI_FINALIZE(ierr)
          stop
      else if((.not. restart) .and. ex) then
          write(*,*) 'not a restart, but '//trim(out_dir)//' exists'
          write(*,*) 'stopping'
          call MPI_FINALIZE(ierr)
          stop
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      mstart = 0
      write(*,*) 'restart status....',restart
      if (restart) then 
         write(*,*) 'opening restart.vars......'


          open(1000+my_rank,file=trim(out_dir)//'restart.vars'//filenum,
     x          status='unknown',
     x          form='unformatted')
          write(*,*) 'reading restart.vars......',filenum
         
          read(1000+my_rank)  b0,b1,b12,b1p2,bt,btc,np,
     x         up,aj,nu,E,input_E,input_p,mstart,input_EeP,
     x         prev_Etot,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,
     x         beta_p,beta_p_buf,wght,beta

          close(1000+my_rank)
          open(1000+my_rank,file=trim(out_dir)//'restart.part'//filenum,
     x         status='unknown',form='unformatted')
          write(*,*) 'reading restart.part......',filenum
          read(1000+my_rank) vp,vp1,vplus,vminus,
     x         xp,Ep,Ni_tot,
     x         Ni_tot_sys,ijkp,
     x         mrat,
     x         xp_buf,vp_buf,Ep_buf,vplus_buf,
     x         vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
     x         B_out_buf,mrat_out_buf,
     x         in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
     x         mrat_buf

          close(1000+my_rank)

         if(my_rank .eq. 0) then
           write(error_unit,*) 'python'
           call execute_command_line(
     x     'python3 fileShrinker.py '//
     x         trim(out_dir)//'grid/ '//int_to_str(mstart)
     x         ,exitstat=ierr)
           if(ierr .ne. 0) then
             write(*,*) 'failed to shrink files1'
             call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr)
             stop
           endif
           call execute_command_line(
     x     'python3 fileShrinker.py '//
     x         trim(out_dir)//'particle/ '//int_to_str(mstart)
     x         ,exitstat=ierr)
           if(ierr .ne. 0) then
             write(*,*) 'failed to shrink files'
             call MPI_ABORT(MPI_COMM_WORLD,ierr,ierr)
             stop
           endif
         endif
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
         write(109) mpu
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
     x     'c.np_3d_1_'//filenum//'.dat', access= acc,
     x     status=stat,form='unformatted')
      open(116,file=trim(out_dir)//'grid/'//
     x     'c.np_3d_2_'//filenum//'.dat', access= acc,
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
     x     'c.temp_p_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')
      open(301,file=trim(out_dir)//'grid/'//
     x     'c.temp_p_3d_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(305,file=trim(out_dir)//'particle/'//
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

      open(330,file=trim(out_dir)//'grid/'//
     x     'c.temp_p_3d_1_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

      open(331,file=trim(out_dir)//'grid/'//
     x     'c.temp_p_3d_2_'//filenum//'.dat',
     x     status=stat, access= acc,
     x     form='unformatted')

c----------------------------------------------------------------------


c======================================================================
c  MAIN LOOP!
c======================================================================

      if (Ni_tot .lt. 0.95*Ni_max) then
        do i=1,100
         call Ionize_pluto_mp(np,np_2,vp,vp1,xp,m,input_p,up)
        enddo
      endif
      do 1 m = mstart+1, nt
         if (my_rank .eq. 0) then
            write(*,*) 'time...', m, dt,mstart
         endif
c         if (m .eq. 300) then
c           call Neut_Center(cx,cy,cz)
c           do i = 1,nx
c              do j = 1,ny 
c                 do k = 1,nz
c                    x = qx(i)-cx
c                    y = qy(j)-cy
c                    z = gz(k)-cz ! global z
c
c                    r = sqrt(x**2 + y**2 + z**2)
c
c                    b0(i,j,k,1) = (3*moment*x*y/r**5)*eoverm
c     x                              + b0(i,j,k,1)
c                    b0(i,j,k,2) = (moment*(3*y**2-r**2)/r**5)*eoverm
c     x                              + b0(i,j,k,1)
c                    b0(i,j,k,3) = (3*moment*z*y/r**5)*eoverm
c     x                              + b0(i,j,k,3)
c                 enddo
c              enddo
c           enddo
c           call f_update_tlev(b1,b12,b1p2,bt,b0)
c         endif

         !Calculate neutral density


         !Ionize cloud and calculate ion density
         write(*,*) 'Ni_tot...',Ni_tot,Ni_max,my_rank

         mr = 1.0/m_pu
         if (Ni_tot .lt. 0.90*Ni_max) then
            call Ionize_pluto_mp(np,np_2,vp,vp1,xp,m,input_p,up)
         endif

         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2
         call update_np_boundary(np)

         !energy diagnostics
         
         call get_bndry_Eflux(b1,E)
         call Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,
     x                    EeP,nu,up,np)

         call curlB(b1,np,aj)
         
         call edge_to_center(bt,btc)

         call extrapol_up(up,vp,vp1,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call improve_up(vp1,vplus,vminus,up,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call get_vp_final(Ep,vp,vp1,vplus)
         
         call move_ion_half(xp,vp,vp1,input_p,Ep)

         call get_Ep_buf(Ep_buf,b0,xp_buf,up)
         call get_vplus_vminus_buf(Ep_buf,vp_buf,vplus_buf,
     x        vminus_buf,b0)
         call get_vp_buf_final(Ep_buf,vp_buf,vplus_buf)
         call move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
         call exchange_ion_in(xp,vp,vp1,input_p,xp_buf,vp_buf)
         call exchange_ion_out(xp,vp,vp1,input_p,xp_buf,vp_buf,
     x        E,Bt,xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf)

         call exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)

         call part_setup_buf(xp_buf,vp_buf)


         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2
         ndiag = ndiag + 1
         if (ndiag .eq. nout) then         
            call get_temperature(xp,vp,np,temp_p)
            mr = 1.0
            call separate_np(np_1,mr)
            mr = 1.0/m_pu
            call separate_np(np_2,mr)
         endif
         call update_np_boundary(np)

         
c**********************************************************************
c SUBCYCLING LOOP!
c**********************************************************************

         dtsub = dtsub_init
         ntf = ntsub
call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
         call check_time_step(bt,np)

         count = 1

         call MPI_ALLREDUCE(ntf,recvbuf,count,
     x        MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)

         ntf = recvbuf

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do 2 n = 1, int(ntf)

         !convert main cell covarient bt to main cell contravarient
         call curlB(b1,np,aj)     


         call predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu) 
         call correct_B(b0,b1,b1p2,E,aj,up,np,nu)


         call f_update_tlev(b1,b12,b1p2,bt,b0)


 2     continue
c**********************************************************************


         call move_ion_half(xp,vp,vp1,input_p,Ep)

         call move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
         
         call exchange_ion_in(xp,vp,vp1,input_p,xp_buf,vp_buf)
         call exchange_ion_out(xp,vp,vp1,input_p,xp_buf,vp_buf,
     x        E,Bt,xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf)

         call exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)

         call part_setup_buf(xp_buf,vp_buf)

         call check_min_den(np,xp,vp,vp1,up,bt)

         if (Ni_tot .lt. 0.9*Ni_max) then
            call res_chex(xp,vp,vp1)
         endif



c----------------------------------------------------------------------
c diagnostic output
c----------------------------------------------------------------------

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)



         ndiag_part = ndiag_part + 1
         if (ndiag .eq. nout) then


            nproc_2rio = nint(100*rio/(delz*nz))

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c save 3d arrays------------------------
               write(111) m
               write(111) np
c               write(115) m
c               write(115) np_1
c               write(116) m
c               write(116) np_2
               write(131) m
               write(131) b1
               write(132) m
               write(132) b0
               write(133) m
               write(133) bt
c               write(140) m
c               write(140) aj
c               write(150) m
c               write(150) E
               write(181) m
               write(181) up
               write(301) m
               write(301) temp_p/1.6e-19
               if ( ndiag_part .eq. part_nout ) then
                   write(305) m
                   write(305) xp
                   write(310) m
                   write(310) vp
                   write(311) m
                   write(311) tags
                   write(315) m
                   write(315) beta_p
                   write(320) m
                   write(320) mrat
                   ndiag_part = 0
                endif
c               write(330) m 
c               write(330) temp_p_1/1.6e-19
c               write(331) m 
c               write(331) temp_p_2/1.6e-19


c save 2d arrays----------------------
c               write(110) m
c               write(110) np(:,ny/2,:),np(:,:,2)
c               write(115) m
c               write(115) np_1(:,ny/2,:),np_1(:,:,2)
c               write(116) m
c               write(116) np_2(:,ny/2,:),np_2(:,:,2)
c               write(130) m
c               write(130) bt(:,ny/2,:,:),bt(:,:,2,:)
c               write(140) m
c               write(140) aj
c               write(150) m
c               write(150) E
c               write(180) m
c               write(180) up(:,ny/2,:,:),up(:,:,2,:)
c               write(300) m
c               write(300) temp_p(:,ny/2,:)/1.6e-19,
c     x                    temp_p(:,:,2)/1.6e-19
               ndiag = 0

         endif

c----------------------------------------------------------------------


c----------------------------------------------------------------------
c Write restart file
c----------------------------------------------------------------------

         if (m .eq. restart_counter) then


          write(*,*) 'writing restart file....',
     x                 'restart.part'//filenum//'.new',my_rank,cart_rank
          open(1000+my_rank,file=trim(out_dir)//'restart.vars'//filenum,
     x              status='unknown',
     x              form='unformatted')
         
          write(1000+my_rank)  b0,b1,b12,b1p2,bt,btc,np,
     x             up,aj,nu,E,input_E,input_p,m,input_EeP,
     x             prev_Etot,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,
     x             beta_p,beta_p_buf,wght,beta

          close(1000+my_rank)
          open(1000+my_rank,file=trim(out_dir)//'restart.part'//filenum,
     x             status='unknown',form='unformatted')
          write(1000+my_rank) vp,vp1,vplus,vminus,
     x             xp,Ep,Ni_tot,
     x             Ni_tot_sys,ijkp,
     x             mrat,
     x             xp_buf,vp_buf,Ep_buf,vplus_buf,
     x             vminus_buf,xp_out_buf,vp_out_buf,E_out_buf,
     x             B_out_buf,mrat_out_buf,
     x             in_bounds,Ni_tot_buf,in_bounds_buf,Ni_tot_out_buf,
     x             mrat_buf

          close(1000+my_rank)
                  
          restart_counter = restart_counter + mrestart


         endif

c----------------------------------------------------------------------


         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


 1     continue
c======================================================================


          close(110)
          close(115)
          close(116)
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
          close(305)
          close(310)
          close(315)
          close(320)
          close(330)
          close(331)
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
