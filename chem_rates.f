
      MODULE chem_rates

      USE global
      USE inputs
      USE gutsp_dd
      USE grid_interp
      
      contains

c----------------------------------------------------------------------
      real function neut_corona(a,r)
          real a,r
          real pi
          pi = 3.14159
          neut_corona = 
     x         atmosphere(r)*((r/Rpluto)**(2*(tanh(5*(a-pi/2))+1)/2))
      end function neut_corona

      real function atmosphere(r)
          real r
          atmosphere = 1e15*(Rpluto/r)**25.0 + 5e9*(Rpluto/r)**8.0
      end function atmosphere

c---------------------------------------------------------------------
      real FUNCTION neutral_density(i,j,k)
c----------------------------------------------------------------------
c      include 'incurv.h'

      integer i,j,k
      real x,y,z
      real cx,cy,cz
      real r,a
      real nn0
      real cap_r

      cap_r = 1.6*Rpluto

      call Neut_Center(cx,cy,cz)
      x = qx(i)-cx
      y = qy(j)-cy
      z = gz(k)-cz ! global z
      rho2 = y**2 + z**2
      r = sqrt(x**2+rho2)
      a = atan2(sqrt(rho2), x)
      neutral_density=atmosphere(max(r,cap_r)) !neut_corona(a, max(r,cap_r))

      neutral_density = neutral_density*1e15
      return
      end FUNCTION neutral_density
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE res_chex(xp,vp,vp1)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      real cx,cy,cz,r,vrel
      real nn0,nn !,neutral_density
      real chex_tau,chex_prob

      real sigma_chex
      PARAMETER (sigma_chex = 1e-25)  !10 A^2 in units of km^2

      nn0 = Ncol*(pwl+1)/RIo

      call Neut_Center(cx,cy,cz)
      
      do l = 1,Ni_tot 
         vrel = sqrt(vp(l,1)**2 + vp(l,2)**2 + vp(l,3)**2)
c         if (r .ge. RIo) then 
         nn = neutral_density(ijkp(l,1),ijkp(l,2),ijkp(l,3))
c            nn = nn0*(RIo/r)**(pwl)
c         else
c            nn = nn0
c         endif
         chex_tau = 1./(nn*sigma_chex*vrel)
         chex_prob = dt/chex_tau

         if (pad_ranf() .lt. chex_prob) then
            do m=1,3
               input_E = input_E - 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              beta*beta_p(l) 
c               input_p(m) = input_p(m) - m_arr(l)*vp(l,m)/beta
            enddo                     
            
            vp(l,:) = 0.0
            vp1(l,:) = 0.0
c            m_arr(l) = m_pu*mproton
            mrat(l) = 1./m_pu
c            beta_p(l) = beta_pu
            beta_p(l) = 1.0
c            write(*,*) 'chex...',l,chex_prob
         endif


      enddo

      return
      end SUBROUTINE res_chex
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE get_ndot(ndot)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real ndot(nx,ny,nz)
      real ndot2(nx,ny,nz)

      real r, theta, phi, cx,cy,cz
      real x1,y1,z1,dvol,ndot_intgl


      real recvbuf
      integer count
      count = 1


      call Neut_Center(cx,cy,cz)

      ndot_intgl = 0.0
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = gz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
c               theta = atan2(sqrt(x1**2 + y1**2),z1)
c               phi = atan2(y1,x1)
               
c               dvol = dx*dy*dz_grid(k)
c shell distribution
c               ndot(i,j,k) = dvol*exp(-(r - 1.4*RIo)**2/(0.3*RIo)**2)*
c     x                       sin(theta)
c power law distribution
c               if (r .gt. RIo) then
c                  ndot(i,j,k) = dvol*(Rio/r)**3.5*
c     x                 sin(theta)
c               endif
c               if (r .le. RIo) then 
c                  ndot(i,j,k) = 0.0
c               endif

               npofr = vol*beta*neutral_density(i,j,k)*dt/tau_photo
               ndot(i,j,k) = exp(-(r - 1.4*RIo)**2/(0.2*RIo)**2)*
     x                       sin(theta)*(cos(phi)+1)/2 !+
c     x                    0.2*dvol*exp(-(r - 1.2*RIo)**2/(0.1*RIo)**2)*
c     x                       sin(theta)
               ndot_intgl = ndot_intgl + ndot(i,j,k)
            enddo
         enddo
      enddo

      write(*,*) 'ndot_intgl....',ndot_intgl,my_rank
      

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      write(*,*) 'ndot_intgl_global....',recvbuf,my_rank

c      write(*,*) 'ndot_intgl...',ndot_intgl
c      ndot = (Mdot/(m_pu*mproton*recvbuf))*ndot
      ndot = (Mdot/(m_pu*mproton*recvbuf))*ndot/(dx*dy*delz)
      if (ndot_intgl .lt. 0.001) ndot = 0.0

c add a power law contribution to ndot

c      ndot_intgl = 0.0
c      do i = 1,nx
c         do j = 1,ny
c            do k = 1,nz
c               x1 = qx(i) - cx
c               y1 = qy(j) - cy
c               z1 = gz(k) - cz !global coordinate
c               r = sqrt(x1**2 + y1**2 + z1**2)
c               theta = atan2(sqrt(x1**2 + y1**2),z1)
c               phi = atan2(y1,x1)
               
c               dvol = dx*dy*dz_grid(k)
c               if (r .gt. RIo) then
c                  ndot2(i,j,k) = dvol*(Rio/r)**12*
c     x                 sin(theta)*(cos(phi+PI)+1)/2
c               endif
c               if (r .le. RIo) then 
c                  ndot2(i,j,k) = 0.0
c               endif

c               ndot_intgl = ndot_intgl + ndot2(i,j,k)*dvol
c            enddo
c         enddo
c      enddo

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

c      write(*,*) 'ndot_intgl_global....',recvbuf,my_rank

c      write(*,*) 'ndot_intgl...',ndot_intgl
c      ndot2 = (Mdot/(mO*recvbuf))*ndot2
c      if (ndot_intgl .lt. 0.001) ndot2 = 0.0

c combine shell with power law

c      ndot = ndot + ndot2

c      write(*,*) 'Mdot...',sum(ndot)*dx*dy*delz*mO,my_rank

      write(*,*) 'Mdot...',sum(ndot)*(dx*dy*delz)*m_pu*mproton,my_rank,
     x         (ndot(nx/2,ny/2,nz/2))

c      write(*,*) 'Max ndot...',maxval(ndot)
c      stop

      return
      end SUBROUTINE get_ndot
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ndot_gauss(ndot)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real ndot(nx,ny,nz)
      real ndot2(nx,ny,nz)

      real r, theta, phi, cx,cy,cz
      real x1,y1,z1,dvol,ndot_intgl


      real recvbuf
      integer count
      count = 1


      call Neut_Center(cx,cy,cz)

      ndot_intgl = 0.0
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = qz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
               theta = atan2(sqrt(x1**2 + y1**2),z1)
               phi = atan2(y1,x1)
               dvol = dx*dy*dz_grid(k)
               ndot(i,j,k) = exp(-(r)**2/(RIo)**2)
               ndot_intgl = ndot_intgl + ndot(i,j,k)
            enddo
         enddo
      enddo

      write(*,*) 'ndot_intgl....',ndot_intgl,my_rank
      

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      write(*,*) 'ndot_intgl_global....',recvbuf,my_rank

      ndot = (Mdot/(m_pu*mproton*recvbuf))*ndot/(dx*dy*delz)
      if (ndot_intgl .lt. 0.01) ndot = 0.0

      write(*,*) 'Mdot...',sum(ndot)*(dx*dy*delz)*m_pu*mproton,my_rank,
     x         (ndot(nx/2,ny/2,nz/2))

      return
      end SUBROUTINE get_ndot_gauss
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ndot_Xianzhe(ndot,nn)
c     ndot is a density rate (cm-3 s-1)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real ndot(nx,ny,nz)
      real ndot2(nx,ny,nz)
      real nn(nx,ny,nz)

      real r, theta, phi, cx,cy,cz
      real x1,y1,z1,dvol,ndot_intgl

c      real RIo_off !now defined as a parameter in para.h
c      parameter (RIo_off = 0.0)

      real recvbuf
      integer count
      count = 1

      if(my_rank.eq.0) then
          write(*,*)'Scaled ionization rate XIANZHE (maind/) '
          write(*,*)'----------------------------------------------'
          endif

      call Neut_Center(cx,cy,cz)

      ndot_intgl = 0.0
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               x1 = qx(i) - cx
               y1 = qy(j) - cy
               z1 = qz(k) - cz !global coordinate
               r = sqrt(x1**2 + y1**2 + z1**2)
c               theta = atan2(sqrt(x1**2 + y1**2),z1)
c               phi = atan2(y1,x1)
               dvol = dx*dy*dz_grid(k)

               if ((r.ge.(Rio+RIO_off)) .and. r .lt. 3*Rio)   then  ! to compare to Linker values

c try Linker again
c                   ndot(i,j,k)= 5.e6*(r/Rio)**(-3.5)*1.e15/(25.*3600.)  ! no pwl to make it smooth

                     ndot(i,j,k) =  nn(i,j,k)/ (25. * 3600.) !  timescale for Ioniz = 25 hours
                     ndot(i,j,k) = ndot(i,j,k) *49./64. !  to scale to Xianzhe 7e27 s-1
                      endif

               if (r.lt.(Rio+RIO_off)) then
                     ndot(i,j,k)=0.
                     endif

               if( (r.ge.(Rio+RIO_off)).and.(r.le.(6.*Rio))  ) then  ! to compare to Linker values
                     ndot_intgl = ndot_intgl + ndot(i,j,k)*dvol !cm-3 integrated on a cell volume
                     endif
            enddo
         enddo
      enddo

      write(*,*) '    Proc ionioz rate XIANZHE(in 6Rio)= ',
     x           ndot_intgl, my_rank


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(ndot_intgl,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(my_rank.eq.0) then
      write(*,*) '   TOTAL XIANZHE (6 Rio) ndot_intgl_global ',
     x           recvbuf,my_rank
      endif

c     DOLS I remove the scaling to Mdot from para.h
c Linker give a rate and that is it (kg cm-3 s-1) but the code wants only km-3s-1
      return
      end SUBROUTINE get_ndot_Xianzhe
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_nuin(nuin,nn,uf,nf,ndot)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real nuin(nx,ny,nz)
      real nn(nx,ny,nz)
      real uf(nx,ny,nz,3)
      real nf(nx,ny,nz)
      real ndot(nx,ny,nz)
      real ufc(nx,ny,nz,3) !gather at cell center

c      real sigma_in
c      parameter (sigma_in = 3.0e-26)

c      call periodic(uf)
c      call periodic_scalar(nf)

      call face_to_center(uf,ufc)

c      call periodic(ufc)
      
      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
c               uf2(1) = 0.5*(uf(i,j,k,1) + uf(i-1,j,k,1)) 
c               uf2(2) = 0.5*(uf(i,j,k,2) + uf(i,j-1,k,2)) 
c               uf2(3) = 0.5*(uf(i,j,k,3) + uf(i,j,k-1,3)) 
c               nuin(i,j,k) = sqrt(uf2(1)**2 + uf2(2)**2 + 
c     x                       uf2(3)**2)*sigma_in*nn(i,j,k)
               nuin(i,j,k) = sqrt(ufc(i,j,k,1)**2 + ufc(i,j,k,2)**2 + 
     x                       ufc(i,j,k,3)**2)*sigma_in*nn(i,j,k)
c               nuin(i,j,k) = 10.0*ndot(i,j,k)/nf(i,j,k)
            enddo
         enddo
      enddo

      call periodic_scalar(nuin)

      return
      end SUBROUTINE get_nuin
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE Ionize_pluto_mp(np,np_2,vp,vp1,xp,m_tstep,input_p,up)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------
CVD$R VECTOR
c      include 'incurv.h'

      real np(nx,ny,nz),
     x     np_2(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)
c     x     gz(nz)

      real function ranf      
      real uptot1(3),uptot2(3)
c      integer*4 dNi         !# of new born ions created in dt
c      real vol              !volume of grid cell

      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
c      real r_xyz       !radial distance
c      real src(200)         !particle source distribution
c      real Nofr(200)        !number of neutrals as func of r
      real nnofr       !neutral density vs. r
      real npofr       !plasma production rate vs. r
      real Np_total         !total number of ions /s
      real vol         !volume of shell vs. r
      real vol_shell
      real vol_shell_min
      real intgl            !integral
      integer rijk
      real ddni
      integer cnt, l1
c      real neutral_density
      real npmax

      real rho2
      real x,y,z
      real small_beta_r
      real max_r
      small_beta_r = 1.6*Rpluto
      max_r = 200*Rpluto

c      integer*4 ion_cnt(nx,ny,nz)  !keeps running count of ions 
                                   !in cell for calculating the bulk
                                   !flow velocity
c      real vsum(nx,ny,nz,3) !total particle velocity in cell of the
                            !new born ions...for get bulk flow

      call Neut_Center(cx,cy,cz)

c get source density

      
      vol = dx**3
      cnt = 0
      l1 = Ni_tot+1
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1

               x = qx(i)-cx
               y = qy(j)-cy
               z = gz(k)-cz ! global z
               rho2 = y**2 + z**2
               r = sqrt(x**2+rho2)
             
               npmax = sqrt(neutral_density(i,j,k)/(tau_photo*k_rec))

c               if ((r .le. dx*S_radius) .and.
c     x              (np_2(i,j,k) .lt. npmax)) then

                  if (r .gt. max_r) cycle

                  if (r .le. small_beta_r) then
                     bpu = 0.01
                     npofr = vol*beta*bpu*
     x                    neutral_density(i,j,k)*dt/tau_photo
                  else 
                     bpu = 2.0
                     npofr = vol*beta*bpu*
     x                    neutral_density(i,j,k)*dt/tau_photo
                  endif

                  if ((npofr .ge. 1) .and. (npofr+l1 .lt. Ni_max)) then
                     do ll = 1,nint(npofr)
                        l = Ni_tot + 1
                        vp(l,1) = 0.0
                        vp(l,2) = 0.0
                        vp(l,3) = 0.0                        

                        xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                        xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                        xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)

                        ii=0
 26                     continue
                        ii = ii + 1
                        if (xp(l,1) .gt. qx(ii)) go to 26 !find i on non-uniform 
                        ii = ii-1
                        ijkp(l,1)= ii

                        jj=0
 18                     continue
                        jj = jj + 1
                        if (xp(l,2) .gt. qy(jj)) go to 18 !find j on non-uniform 
                        jj = jj-1
                        ijkp(l,2)= jj

                        kk=2
 15                     continue
                        kk = kk + 1
                        if (xp(l,3) .gt. qz(kk)) go to 15 !find k on non-uniform 
                        kk = kk-1
                        ijkp(l,3)= kk

                        mrat(l) = 1.0/m_pu
c                        m_arr(l) = mproton*m_pu
                        beta_p(l) = bpu
                        tags(l) = 1
                        Ni_tot = l
                        cnt = cnt + 1
                        do m=1,3
                           vp1(l,m) = vp(l,m)
                           input_E = input_E + 
     x                      0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x                          (beta*beta_p(l))
                    input_p(m) = input_p(m) + (mion/mrat(l))*vp(l,m)/
     x                          (beta*beta_p(l))
                        enddo                     

                     enddo
                  endif
                  if ((npofr .lt. 1).and.(npofr + l1 .lt. Ni_max)) then
                     if (npofr .gt. pad_ranf()) then
                        l = Ni_tot + 1
                        vp(l,1) = 0.0
                        vp(l,2) = 0.0
                        vp(l,3) = 0.0                        
                        
                        xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                        xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                        xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)


                        ii=0
 27                     continue
                        ii = ii + 1
                        if (xp(l,1) .gt. qx(ii)) go to 27 !find i on non-uniform 
                        ii = ii-1
                        ijkp(l,1)= ii

c                        ijkp(l,1) = floor(xp(l,1)/dx) !particle grid location index
c                        ijkp(l,2) = floor(xp(l,2)/dy)
                        

                        jj=0
 17                     continue
                        jj = jj + 1
                        if (xp(l,2) .gt. qy(jj)) go to 17 !find i on non-uniform 
                        jj = jj-1
                        ijkp(l,2)= jj                     

                        kk=2
 16                     continue
                        kk = kk + 1
                        if (xp(l,3) .gt. qz(kk)) go to 16 !find k on non-uniform 
                        kk = kk-1
                        ijkp(l,3)= kk

c                        kk=1
c                       do 16 while((xp(l,3).gt.qz(kk)).and.(kk .le. nz)) !find ck
c                           ijkp(l,3) = kk !grid
c                           kk=kk+1
c 16                     continue
c                        kk=ijkp(l,3)
cc                        if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
c                          ijkp(l,3) = kk+1
c                        endif
                        
                        mrat(l) = 1.0/m_pu
c                        m_arr(l) = mproton*m_pu
                        beta_p(l) = bpu
                        tags(l) = 1
                        Ni_tot = l
                        cnt = cnt + 1
                        do m=1,3
                           vp1(l,m) = vp(l,m)
                           input_E = input_E + 
     x                          0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2/
     x                          (beta*beta_p(l))
                           input_p(m)=input_p(m)+(mion/mrat(l))*vp(l,m)/
     x                          (beta*beta_p(l))
                        enddo                     
                        
                        
                     endif
                     
                  endif

                  
c               endif
            enddo
         enddo
      enddo

      
c      write(*,*) 'total new ions....',my_rank,cnt,Ni_tot         
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c      stop

      do 60 l = l1,Ni_tot
         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x        (ijkp(l,3) .gt. nz)) then
            write(*,*) 'removing ion...',my_rank,ijkp(l,:)
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue

c      stop

c      write(*,*) 'Ni_tot after wake....',Ni_tot

c      call get_interp_weights(xp)
c      call update_np(np)
c      call update_up(vp,np,up)


      
      return
      end SUBROUTINE Ionize_pluto_mp
c----------------------------------------------------------------------

      end MODULE chem_rates
