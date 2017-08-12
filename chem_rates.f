
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

      cap_r = 1.4*Rpluto

      call Neut_Center(cx,cy,cz)
      x = qx(i)-cx
      y = qy(j)-cy
      z = gz(k)-cz ! global z
      rho2 = y**2 + z**2
      r = sqrt(x**2+rho2)
      a = atan2(sqrt(rho2), x)

      if( r .lt. S_radius ) then
          neutral_density=atmosphere(max(r,cap_r)) !neut_corona(a, max(r,cap_r))
      else
          neutral_density=0
      endif

      neutral_density = neutral_density*1e15
      return
      end FUNCTION neutral_density
c----------------------------------------------------------------------

      SUBROUTINE Ionization(np,vp,vp1,xp,up)
      end SUBROUTINE Ionization

c----------------------------------------------------------------------
      SUBROUTINE res_chex(xp,vp,vp1)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      real cx,cy,cz,r,vrel
      real nn !,neutral_density
      real chex_tau,chex_prob

      real sigma_chex
      PARAMETER (sigma_chex = 1e-25)  !10 A^2 in units of km^2

      call Neut_Center(cx,cy,cz)
      
      do l = 1,Ni_tot 
         vrel = sqrt(vp(l,1)**2 + vp(l,2)**2 + vp(l,3)**2)
         nn = neutral_density(ijkp(l,1),ijkp(l,2),ijkp(l,3))

         chex_tau = 1./(nn*sigma_chex*vrel)
         chex_prob = dt/chex_tau

         if (pad_ranf() .lt. chex_prob) then
            do m=1,3
               input_E = input_E - 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              beta*beta_p(l) 
            enddo                     
            
            vp(l,:) = 0.0
            vp1(l,:) = 0.0
            mrat(l) = 1./m_pu
            beta_p(l) = 1.0
         endif


      enddo

      return
      end SUBROUTINE res_chex
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE Ionize_pluto_mp(np,vp,vp1,xp,m_tstep,input_p,up)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------

      real np(nx,ny,nz),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     xp(Ni_max,3),
     x     input_p(3),
     x     up(nx,ny,nz,3)
c     x     gz(nz)

      real function ranf      
      real uptot1(3),uptot2(3)
      real r                !dist of particle to neutral cloud center
      real t                !run time
      real v                !neutral cloud velocity, r/t
      real cx,cy,cz         !x,y,z coord of neutral cloud center
      real theta,phi        !spherical coords used to get velocities
      integer flg           !flag for while loop
      real zsum             !sum along z to find cell for new_borns
      real rnd              !random number
      real n_source,vol_source
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
      real npmax

      real rho2
      real x,y,z
      real small_beta_r
      real max_r
      small_beta_r = 1.4*Rpluto


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
                        
                        mrat(l) = 1.0/m_pu
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

                  
            enddo
         enddo
      enddo

      

      do 60 l = l1,Ni_tot
         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x        (ijkp(l,3) .gt. nz)) then
            write(*,*) 'removing ion...',my_rank,ijkp(l,:)
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue

      return
      end SUBROUTINE Ionize_pluto_mp
c----------------------------------------------------------------------

      end MODULE chem_rates
