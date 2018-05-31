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
            Evp = Evp + 0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
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

      write(*,*) 'Normalized particle energy...',S_Evp/S_input_E
      write(*,*) 'Normalized total energy......',actual_E/supposed_E

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
c      include 'incurv.h'

      real vthm
      integer flg
      real v,vx,vy,vz
      real rnd,f

      vx = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
      vy = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
      vz = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())

      return
      end SUBROUTINE maxwl_init
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl(np,vp,vp1,xp,up)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v,tmp
      real rand
      real vx,vy,vz
      real dvx,dvz,v1

      integer flg
      real nprat
      real bwght
      real mpart

      v1 = 1.0


      nprat = np_bottom/np_top

      Ni_tot = Ni_tot + Ni_tot*nint((1/nprat)-1)/2

c initialize protons

      bwght = 1.0
      vth = vth_top
      mpart = mproton

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         

         i=0
 11      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 11 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         j=0
 13      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 13 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j

         
         k=2
 12      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 12 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k


         call maxwl_init(vth,vx,vy,vz)

         ii = ijkp(l,1)
         kk = ijkp(l,3)

         vp(l,1) = -vsw + vx 
         vp(l,2) = vy 
         vp(l,3) = vz

         mrat(l) = mproton/mpart
         beta_p(l) = bwght

         do 20 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*bwght)
 20      continue
                 
 10      continue

c initialize He++

      Ni_tot_1 = Ni_tot + 1
      
      Ni_tot = Ni_tot + f_mq_2*Ni_tot

      do 30 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         

         i=0
 31      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         j=0
 33      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 33 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j


         k=2
 32      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 32 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k

         vth = vth_top

         call maxwl_init(vth,vx,vy,vz)

         vp(l,1) = -vsw + vx 
         vp(l,2) = vy 
         vp(l,3) = vz

         mrat(l) = 1.0/2.0
         beta_p(l) = b_mq_2

         do 40 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
 40      continue
 30   continue


c add shell distribution

c         Ni_tot_1 = Ni_tot + 1
c
c         Ni_tot = Ni_tot + f_shl*Ni_tot
c
c         do 69 l = Ni_tot_1,Ni_tot
c
c
c            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
c            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
c            xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
c            
c            mrat(l) = 1.0
c            beta_p(l) = b_shl
c
c
c            i=0
c 71         continue
c            i = i + 1
c            if (xp(l,1) .gt. qx(i)) go to 71 !find i on non-uniform 
c            i = i-1
c            ijkp(l,1)= i
c            
c            
c            j=0
c 73         continue
c            j = j + 1
c            if (xp(l,2) .gt. qy(j)) go to 73 !find i on non-uniform 
c            j = j-1
c            ijkp(l,2)= j
c
c            
c            k=2
c 72         continue
c            k = k + 1
c            if (xp(l,3) .gt. qz(k)) go to 72 !find k on non-uniform 
c            k = k-1
c            ijkp(l,3)= k
c
c            
c            vz = pad_ranf()*2 - 1
c            tmp = sqrt(1-vz**2)
c            phi = 2*PI*pad_ranf()
c            vx = tmp*cos(phi)
c            vy = tmp*sin(phi)
c
c            vp(l,1) = -vsw+vsw*vx
c            vp(l,2) = vsw*vy
c            vp(l,3) = vsw*vz
c
c            do  m=1,3
c               vp1(l,m) = vp(l,m)
c               input_E = input_E + 
c     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
c     x              (beta*b_shl)
c            enddo
c            
c
c 69      enddo
      tags(1:Ni_tot) = 1



      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(xp, vp, vp1, np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end SUBROUTINE sw_part_setup_maxwl
c----------------------------------------------------------------------

      end MODULE part_init









