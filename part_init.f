      MODULE part_init


      USE global
      USE dimensions
      USE misc
      USE gutsp_dd
      USE mpi
      USE iso_fortran_env, only: output_unit,error_unit

      contains


c----------------------------------------------------------------------
      SUBROUTINE Energy_diag(vp,b0,b1,E,Evp,EB1,
     x                       EE,nu,up,np)
c----------------------------------------------------------------------

      real vp(Ni_max,3)
      real b0(nx,ny,nz,3)
      real b1(nx,ny,nz,3)
      real E(nx,ny,nz,3)
      real Evp                  !kinetic energy of particles
      real EB1                  !Magnetic field energy 
      real EE                   !Electric field energy 
      real nu(nx,ny,nz)
      real up(nx,ny,nz,3)
      real np(nx,ny,nz)


      real actual_E
      real supposed_E
      real vol                  !volume of cell

      real recvbuf              !buffer for mpi calls

      integer c                 !counter
      real mO_q                 !proton mass to charge ratio

      integer i,j,k,l,m
      integer ierr

      real S_EB1
      real S_EE
      real S_Evp
      real S_input_E
      real S_bndry_Eflux

      c = 1
      mO_q = mion/q

      ! Initialize energies to zero before accumulating
      EB1 = 0.0
      EE = 0.0
      Evp = 0.0

      ! First compute field energies
      do 10 i=1,nx-1
         do 10 j=1,ny-1
            do 10 k=1,nz-1
               vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
               do 10 m=1,3
                  EB1 = EB1 + 
     x              (vol/(2.0*mu0))*(mO_q*b1(i,j,k,m))**2
                  EE = EE + (epsilon0*vol/2.0)*
     x                      (mO_q*E(i,j,k,m)*km_to_m)**2
 10               continue

      ! Then compute particle (kinetic) energy
      do 15 l=1,Ni_tot
         do 15 m=1,3
            Evp = Evp + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
 15   continue

      ! Add up global magnetic field energy
      call MPI_ALLREDUCE(EB1,recvbuf,c,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_EB1 = recvbuf


      ! Add up global electric field energy
      call MPI_ALLREDUCE(EE,recvbuf,c,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_EE = recvbuf


      ! Add up global particle energy
      call MPI_ALLREDUCE(Evp,recvbuf,c,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_Evp = recvbuf


      ! Add up global input energy
      call MPI_ALLREDUCE(input_E,recvbuf,c,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_input_E = recvbuf


      ! Add up boundary energy flux
      call MPI_ALLREDUCE(bndry_Eflux,recvbuf,c,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_bndry_Eflux = recvbuf

      ! Total energy that exists in the simulation on this step
      ! (kinetic + magnetic + electric)
      actual_E = S_Evp + S_EB1 + S_EE

      ! Total energy that should exist in the simulation
      ! based on what was put in and taken out.
      ! (kinetic + electromagnetic)
      supposed_E = S_input_E + S_bndry_Eflux

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(output_unit)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (my_rank .eq. 0) then

      NE_part = S_Evp/S_input_E
      NE_total = actual_E/supposed_E
      write(*,*) 'Particle energy.........', S_Evp
      write(*,*) 'Input particle energy...', S_input_E
      write(*,*) 'Magnetic energy.........', S_EB1
      write(*,*) 'Electric energy.........', S_EE
      write(*,*) 'Total boundary E flux...', S_bndry_Eflux
      write(*,*) 'Normalized particle energy...',NE_part
      write(*,*) 'Normalized total energy......',NE_total

      if (isnan(NE_part) .or. isnan(NE_total)) then
          write(error_unit,*) 'Energy is NaN'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          stop
      endif

      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(output_unit)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      return
      end SUBROUTINE Energy_diag
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE maxwl_init(vthm,vx,vy,vz)
c----------------------------------------------------------------------

      real vthm
      real v,vx,vy,vz
      real rnd,f

      vx = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
      vy = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
      vz = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())

      return
      end SUBROUTINE maxwl_init
c----------------------------------------------------------------------

      SUBROUTINE locate_new(l,x,y,z)
      integer l
      real x,y,z
      integer i,j,k
         i=0
 11      continue
         i = i + 1
         if (x .gt. qx(i)) go to 11 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         j=0
 13      continue
         j = j + 1
         if (y .gt. qy(j)) go to 13 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j

         
         k=2
 12      continue
         k = k + 1
         if (z .gt. qz(k)) go to 12 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k
      end SUBROUTINE locate_new

      SUBROUTINE insert_maxwl_buf(start_N, stop_N, drift, vthrm, vp_buf,
     x        xp_buf, mr, b, t)
      integer start_N, stop_N
      real vp_buf(Ni_max_buf,3)
      real xp_buf(Ni_max_buf,3)
      integer l, m
      real mr, b, t !mrat, beta, and tag
      real drift, vthrm

      do l = start_N,stop_N
         xp_buf(l,1) = qx(nx)+(1.0-pad_ranf())*dx_buf
         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         call maxwl_init(vthrm,vp_buf(l,1),vp_buf(l,2),vp_buf(l,3))
         vp_buf(l,1) = vp_buf(l,1) + drift

         mrat_buf(l) = mr
         beta_p_buf(l) = b
         tags_buf(l) = t
      enddo
      end SUBROUTINE insert_maxwl_buf

      SUBROUTINE insert_maxwl(start_N, stop_N, drift, vthrm, vp, vp1,
     x        xp, mr, b, t)
      integer start_N, stop_N
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      integer l, m
      real mr, b, t !mrat, beta, and tag
      real drift, vthrm
      integer ierr


      do l = start_N,stop_N
         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         call maxwl_init(vthrm,vp(l,1),vp(l,2),vp(l,3))
         vp(l,1) = vp(l,1) + drift

         call locate_new(l, xp(l,1), xp(l,2), xp(l,3))

         mrat(l) = mr
         beta_p(l) = b
         tags(l) = t

         do 20 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
 20      continue
      enddo
      end SUBROUTINE insert_maxwl

      SUBROUTINE insert_shell_buf(start_N, stop_N, drift, vinj, vp_buf,
     x        xp_buf, mr, b, t)
      integer start_N, stop_N
      real vp_buf(Ni_max_buf,3)
      real xp_buf(Ni_max_buf,3)
      integer l, m
      real mr, b, t !mrat, beta, and tag
      real drift, vinj

      do l = start_N,stop_N
         xp_buf(l,1) = qx(nx)+(1.0-pad_ranf())*dx_buf
         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         call sphere(vinj,vp_buf(l,1),vp_buf(l,2),vp_buf(l,3))
         vp_buf(l,1) = vp_buf(l,1) + drift

         mrat_buf(l) = mr
         beta_p_buf(l) = b
         tags_buf(l) = t
      enddo
      end SUBROUTINE insert_shell_buf

      SUBROUTINE insert_shell(start_N, stop_N, drift, vinj, C, vp, vp1,
     x        xp, mr, b, t)
      integer start_N, stop_N
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      integer l, m
      real mr, b, t !mrat, beta, and tag
      real drift, vinj
      real C

      do l = start_N,stop_N
         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         ! Thin shell!
         call sphere(vinj, vp(l,1),vp(l,2),vp(l,3))
         vp(l,1) = vp(l,1) + drift
         call locate_new(l, xp(l,1), xp(l,2), xp(l,3))

         mrat(l) = mr
         beta_p(l) = b
         tags(l) = t

         do 20 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
 20      continue
      enddo
      end SUBROUTINE insert_shell

c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup(np,vp,vp1,xp,up)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real np(nx,ny,nz)
      real up(nx,ny,nz,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)

      real phi, theta, tmp
      real vx,vy,vz
      integer i,j,k,l,m

      vth = vth_top

c Initialize Thermal H+
      call insert_maxwl(1, Ni_thermal_H, -vsw, vth, vp, vp1,
     x        xp, 1.0, b_sw_thermal_H, sw_thermal_H_tag)

cc initialize Thermal He++
c      call insert_maxwl(Ni_thermal_H+1, Ni_thermal_H+Ni_thermal_He,
c     x  -vsw, vth, vp, vp1, xp, 0.5, b_sw_thermal_He, 
c     x  sw_thermal_He_tag)
c
cc initialize Shell H+
c      call insert_shell(Ni_thermal_H+Ni_thermal_He+1,
c     x                  Ni_thermal_H+Ni_thermal_He+Ni_shell_H,
c     x                  -vsw, vsw, 0.2544, vp, vp1, xp, 1.0, 
c     x                  b_sw_shell_H, sw_shell_H_tag)
c
cc initialize Shell He+
c      call insert_shell(Ni_thermal_H+Ni_thermal_He+Ni_shell_H+1,
c     x                  Ni_tot,
c     x                  -vsw, vsw, 0.02038, vp, vp1, xp, 0.25, 
c     x                  b_sw_shell_He, sw_shell_He_tag)


      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(xp, vp, vp1, np)
      write(*,*) 'update_up...'
      call update_up(vp,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup
c----------------------------------------------------------------------

      end MODULE part_init









