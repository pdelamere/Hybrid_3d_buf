      MODULE part_init


      USE global
      USE dimensions
      USE misc
      USE gutsp_dd
      USE mpi

      contains


c----------------------------------------------------------------------
      SUBROUTINE Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,
     x                       EE,EeP,nu,up,np)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vp(Ni_max,3),
c     x     uf(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
c     x     etemp(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz)

      real mO_q


      real Evp                  !kinetic energy of particles
      real Euf                  !kinetic energy of fluid flow
      real EB1,EB1x,EB1y,EB1z   !Magnetic field energy 
      real EE                   !Electric field energy 
      real EeP                  !Electron pressure energy
      real total_E              !total energy
      real aveEvp               !average particle energy
      real norm_E               !normalized energy
      real vol                  !volume of cell
      real denf                 !fluid density

      real recvbuf
      integer count
      count = 1

      mO_q = mion/q

      Euf = 0.0
      EB1 = 0.0
      EB1x = 0.0
      EB1y = 0.0
      EB1z = 0.0
      EE = 0.0
      EeP = 0.0


      do 10 i=1,nx-1
c         j = 2
         do 10 j=1,ny-1
            do 10 k=1,nz-1
               vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
               EB1x = EB1x + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,1))**2 
               EB1y = EB1y + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,2))**2 
               EB1z = EB1z + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,3))**2 
c               EeP = EeP + kboltz*etemp(i,j,k)
               do 10 m=1,3
                  denf = np(i,j,k)/(km_to_m**3)
                  Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
c                  EB1 = EB1 + 
c     x              (vol/(2.0*mu0))*(mO_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                  EB1 = EB1 + 
     x              (vol/(2.0*mu0))*(mO_q*b1(i,j,k,m))**2
                  EE = EE + (epsilon0*vol/2.0)*
     x                      (mO_q*E(i,j,k,m)*km_to_m)**2
 10               continue

c      input_EeP = input_EeP + EeP

c      write(*,*) 'Energy diag...',Ni_tot,m_arr(2000000)
 
      Evp = 0.0
      do 15 l=1,Ni_tot
         do 15 m=1,3
            Evp = Evp + 0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x           (beta*beta_p(l))
 15   continue

c      write(*,*) 'Energy diag 2...',Ni_tot,m_arr(2000000)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(Evp,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_Evp = recvbuf

c      write(*,*) 'recvbuf...',recvbuf,Evp

      call MPI_ALLREDUCE(input_E,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_input_E = recvbuf

      call MPI_ALLREDUCE(bndry_Eflux,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_bndry_Eflux = recvbuf

c      total_E = S_Evp+EE+EB1
      total_E = S_Evp+EB1
      aveEvp = S_Evp/S_input_E
      
      if (my_rank .eq. 0) then

c      write(*,*) 'Input energy (J).............',S_input_E
cc      write(*,*) 'Input EeP energy (J).........',input_EeP
c      write(*,*) 'Total vp energy (J)..........',S_Evp
c      write(*,*) 'Total up energy (J)..........',Euf
c      write(*,*) 'Total B energy (J)...........',EB1/S_input_E
c      write(*,*) 'Total E energy (J)...........',EE/S_input_E
cc      write(*,*) 'Total EeP energy (J).........',EeP
c      write(*,*) 'Total energy (J).............',total_E
cc      write(*,*) 'Total energy w/ eP (J).......',total_E+EeP
c      write(*,*) 'Energy thru boundaries.......',bndry_Eflux/S_input_E
      write(*,*) 'Normalized particle energy...',aveEvp
      write(*,*) 'Normalized energy............',total_E/S_input_E,
     x   my_rank
      write(*,*) 'Normalized energy (bndry)....',
c     x                S_bndry_Eflux/total_E
     x                (total_E)/(S_input_E+S_bndry_Eflux)
c      write(*,*) 'Normalized energy (no b1z)...',(S_Evp+Euf+EE+EB1x+
c     x                                            EB1y)/S_input_E
cc      write(*,*) 'Normalized energy (w/ eP)....',
cc     x                             (total_E+EeP)/(input_E + input_EeP)
c      write(*,*) ' '

      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      norm_E = total_E/S_input_E

c      if (prev_Etot .eq. 0.0) then prev_Etot = norm_E
c      do 20 i=1,nx 
c         do 20 j=1,ny
c            do 20 k=1,nz
c               nu(i,j,k) = nu(i,j,k) + 
c     x                 nu(i,j,k)*2.0*((norm_E - prev_Etot)/norm_E)
c 20            continue
      prev_Etot = norm_E

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

c      vx = vthm*sqrt(-2*alog(1.0-pad_ranf()))*cos(2*PI*pad_ranf())
c      vy = vthm*sqrt(-2*alog(1.0-pad_ranf()))*cos(2*PI*pad_ranf())
c      vz = vthm*sqrt(-2*alog(1.0-pad_ranf()))*cos(2*PI*pad_ranf())


      vx = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
      vy = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
      vz = vthm*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())


c      flg = 0
c      do 40 while (flg .eq. 0)
c         v = (2*vth_max*pad_ranf())-vth_max
c         f = exp(-(v)**2 / vthm**2)
c         rnd = pad_ranf()
c         if (f .ge. rnd) then 
c            flg = 1
c            vx = v
c         endif
c 40   continue
c      flg = 0
c      do 42 while (flg .eq. 0)
c         v = (2*vth_max*pad_ranf())-vth_max
c         f = exp(-(v)**2 / vthm**2)
c         rnd = pad_ranf()
c         if (f .ge. rnd) then 
c            flg = 1
c            vy = v
c         endif
c 42   continue
c      flg = 0
c      do 44 while (flg .eq. 0)
c         v = (2*vth_max*pad_ranf())-vth_max
c         f = exp(-(v)**2 / vthm**2)
c         rnd = pad_ranf() 
c         if (f .ge. rnd) then 
c            flg = 1
c            vz = v
c         endif
c 44   continue

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
c      integer np_t_flg(Ni_max)
c      integer np_b_flg(Ni_max)

      integer flg
      real nprat
      real bwght
      real mpart

      v1 = 1.0

c      np_t_flg(:) = 0
c      np_b_flg(:) = 0


      nprat = np_bottom/np_top

c      Ni_tot_O = Ni_tot*(2./3.)
c      Ni_tot_S = Ni_tot*(1./3.)

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


c         ijkp(l,2) = floor(xp(l,2)/dy) 

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


c         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c         ijkp(l,2) = nint(xp(l,2)/dy)
         
c         k=1
c         do 12 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c            ijkp(l,3) = k       !grid
c            k=k+1
c 12      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif

         call maxwl_init(vth,vx,vy,vz)

         ii = ijkp(l,1)
         kk = ijkp(l,3)

         vp(l,1) = -vsw + vx 
         vp(l,2) = vy 
         vp(l,3) = vz 

c         m_arr(l) = mpart
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


c         ijkp(l,2) = floor(xp(l,2)/dy) 

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

c         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c         ijkp(l,2) = nint(xp(l,2)/dy)
         
c         k=1
c         do 32 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c            ijkp(l,3) = k       !grid
c            k=k+1
c 32      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         call maxwl_init(vth,vx,vy,vz)

c         ii = ijkp(l,1)
c         kk = ijkp(l,3)

         vp(l,1) = -vsw + vx 
         vp(l,2) = vy 
         vp(l,3) = vz 

c         m_arr(l) = 2*mproton
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

         Ni_tot_1 = Ni_tot + 1

         Ni_tot = Ni_tot + f_shl*Ni_tot

         do 69 l = Ni_tot_1,Ni_tot


            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
            xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
            
c            m_arr(l) = mproton
            mrat(l) = 1.0
            beta_p(l) = b_shl


            i=0
 71         continue
            i = i + 1
            if (xp(l,1) .gt. qx(i)) go to 71 !find i on non-uniform 
            i = i-1
            ijkp(l,1)= i
            
            
c            ijkp(l,2) = floor(xp(l,2)/dy) 


            j=0
 73         continue
            j = j + 1
            if (xp(l,2) .gt. qy(j)) go to 73 !find i on non-uniform 
            j = j-1
            ijkp(l,2)= j

            
            k=2
 72         continue
            k = k + 1
            if (xp(l,3) .gt. qz(k)) go to 72 !find k on non-uniform 
            k = k-1
            ijkp(l,3)= k

c            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
c            ijkp(l,2) = nint(xp(l,2)/dy)
            
c            k=1
c            do 70 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c               ijkp(l,3) = k    !grid
c               k=k+1
c 70         continue
c            k=ijkp(l,3)
c            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c               ijkp(l,3) = k+1
c            endif
            
c            ii = ijkp(l,1)
c            kk = ijkp(l,3)

            vz = pad_ranf()*2 - 1
            tmp = sqrt(1-vz**2)
            phi = 2*PI*pad_ranf()
            vx = tmp*cos(phi)
            vy = tmp*sin(phi)
c            vp(l,1) = 1.0*vsw*cos(theta) + vx !+dvx
c            vp(l,2) = vy 
c            vp(l,3) = 1.0*vsw*sin(theta) + vz        !+dvz 

            vp(l,1) = -vsw+vsw*vx !+dvx
            vp(l,2) = vsw*vy !+dvz 
            vp(l,3) = vsw*vz

c            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
c            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              (beta*b_shl)
            enddo
            

 69      enddo
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

c----------------------------------------------------------------------
      SUBROUTINE load_Maxwellian(np,vp,vp1,xp,up,vth,Ni_tot_1)
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real up(nx,ny,nz,3)
c      real phi,theta,rnd,f,v
c      real rand
c      real vx,vy,vz
c      real dvx,dvz,v1
c      integer np_t_flg(Ni_max)
c      integer np_b_flg(Ni_max)
      integer Ni_tot_1
      real vth

c      integer flg
c      real nprat

      v1 = 1.0

c      np_t_flg(:) = 0
c      np_b_flg(:) = 0

      do 10 l = 1,Ni_tot_1

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
c         m_arr(l) = mion
         mrat(l) = 1.0

c         ijkp(l,1) = floor(xp(l,1)/dx) 

         i=0
 31      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


c         ijkp(l,2) = floor(xp(l,2)/dy) 


         j=0
 33      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 33 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j
        
         k=0
 30      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 30 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k
         
c         vth = vth_bottom

         vx = vsw+vth*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
         vy = vth*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())
         vz = vth*sqrt(-alog(1.0-pad_ranf()))*cos(PI*pad_ranf())

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
     x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
         vp(l,2) = vy 
         vp(l,3) = vz 

c         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
c         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /beta
 45      continue
 10      continue

      call get_interp_weights(xp)
      call update_np(xp, vp, vp1, np)
      call update_up(vp,np,up)

      return
      end SUBROUTINE load_Maxwellian
c----------------------------------------------------------------------





      end MODULE part_init









