
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
          real, intent(in) :: r
          real rr
          rr = max(r, 1.1*Rpluto)
          atmosphere = 1e15*(Rpluto/rr)**25.0 + 5e9*(Rpluto/rr)**8.0
          atmosphere = atmosphere*1e15
      end function atmosphere

c---------------------------------------------------------------------
      real FUNCTION neutral_density(i,j,k)
      integer i,j,k
      neutral_density = neutral_density_continuous(qx(i),qy(j),qz(k))
      return
      end FUNCTION neutral_density
c---------------------------------------------------------------------

      real FUNCTION neutral_density_continuous(x,y,z)
      real x,y,z
      real xx,yy,zz
      real cx,cy,cz

      call Neut_Center(cx,cy,cz)
      xx = x - cx
      yy = y - cy
      zz = z - cz
      r = sqrt(xx**2 + yy**2 + zz*2)
      if( r .lt. S_radius ) then
          neutral_density_continuous = atmosphere(r)
      else
          neutral_density_continuous = 0
      endif

      return
      end FUNCTION neutral_density_continuous
c----------------------------------------------------------------------

      SUBROUTINE Ionization(np,xp,vp,vp1)
          real np(nx,ny,nz)
          real xp(Ni_max,3)
          real vp(Ni_max,3)
          real vp1(Ni_max,3)
          call charge_exchange_ionization(xp,vp,vp1)
          call photoionization(np,xp,vp,vp1)
      end SUBROUTINE Ionization

c----------------------------------------------------------------------
      SUBROUTINE charge_exchange_ionization(xp,vp,vp1)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      real cx,cy,cz,r,vrel
      real nn !,neutral_density
      real chex_inv_tau,chex_prob

      real sigma_chex
      PARAMETER (sigma_chex = 1e-25)  !10 A^2 in units of km^2

      call Neut_Center(cx,cy,cz)
      
      do l = 1,Ni_tot 
         vrel = sqrt(vp(l,1)**2 + vp(l,2)**2 + vp(l,3)**2)
         nn = neutral_density_continuous(xp(l,1),xp(l,2),xp(l,3))

         ! one over the time constant
         chex_inv_tau = (nn*sigma_chex*vrel)

         chex_prob = dt*chex_inv_tau

         if (pad_ranf() .lt. chex_prob) then
            do m=1,3
               input_E = input_E - 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              beta*beta_p(l) 
            enddo                     
            
            vp(l,:) = 0.0
            vp1(l,:) = 0.0
            mrat(l) = 1./m_pu
            ! beta_p stays the same
            tags(l) = pluto_chex_CH4_tag
         endif


      enddo

      return
      end SUBROUTINE charge_exchange_ionization
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE photoionization(np,xp,vp,vp1)
c Ionizes the neutral cloud with a 28 s time constant and fill particle
c arrays, np, vp, up (ion particle density, velocity, 
c and bulk velocity).   
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

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
      real new_macro   !fractional part indicates probabalisitc particle
      real new_micro
      real cur_micro
      
      real Np_total         !total number of ions /s
      real vol         !volume of shell vs. r
      real vol_shell
      real vol_shell_min
      real intgl            !integral
      integer rijk
      integer l1

      real rho2
      real x,y,z

      real pu_beta_p
      real pu_tag



      call Neut_Center(cx,cy,cz)

c get source density

      
      vol = dx**3
      l1 = Ni_tot+1
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1

               if(np(i,j,k) .ge. max_ion_density) then
                   ! Note that this check is not needed. whenever it is
                   ! true new_macro would be negative and no new
                   ! particles would be created anyway.
                   continue
               endif

               x = qx(i)-cx
               y = qy(j)-cy
               z = gz(k)-cz ! global z
               rho2 = y**2 + z**2
               r = sqrt(x**2+rho2)
             
               if(r .le. 2) then
                   pu_beta_p = 0.1
                   pu_tag = pluto_stagnant_photoionize_CH4_tag
               else
                   pu_beta_p = 2.0
                   pu_tag = pluto_photoionize_CH4_tag
               endif

               cur_micro = np(i,j,k)*vol
               new_micro = vol*neutral_density(i,j,k)*dt/tau_photo
               if((cur_micro+new_micro)/vol .le. max_ion_density) then
                   new_macro = new_micro*beta*pu_beta_p
               else
                   new_macro = (max_ion_density*vol-cur_micro)*beta*pu_beta_p
               endif

               do ll = 1,min(nint(new_macro), Ni_max - Ni_tot)
                  l = Ni_tot + 1
                  vp(l,1) = 0.0
                  vp(l,2) = 0.0
                  vp(l,3) = 0.0                        

                  xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                  xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                  xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)

                  ii=0
 26               continue
                  ii = ii + 1
                  if (xp(l,1) .gt. qx(ii)) go to 26 !find i on non-uniform 
                  ii = ii-1
                  ijkp(l,1)= ii

                  jj=0
 18               continue
                  jj = jj + 1
                  if (xp(l,2) .gt. qy(jj)) go to 18 !find j on non-uniform 
                  jj = jj-1
                  ijkp(l,2)= jj

                  kk=2
 15               continue
                  kk = kk + 1
                  if (xp(l,3) .gt. qz(kk)) go to 15 !find k on non-uniform 
                  kk = kk-1
                  ijkp(l,3)= kk

                  mrat(l) = 1.0/m_pu
                  beta_p(l) = pu_beta_p
                  tags(l) = pu_tag
                  Ni_tot = l
                  do m=1,3
                     vp1(l,m) = vp(l,m)
                     ! Add the kinetic energy of the particle
                     ! (often zero)
                     input_E = input_E + 
     x               0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x                     (beta*beta_p(l))
                  enddo                     

               enddo
               new_macro = new_macro - nint(new_macro)
               if (Ni_tot+1 .le. Ni_max) then
                   if (new_macro .gt. pad_ranf()) then
                      l = Ni_tot + 1
                      vp(l,1) = 0.0
                      vp(l,2) = 0.0
                      vp(l,3) = 0.0                        
                      
                      xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                      xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                      xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)


                      ii=0
 27                   continue
                      ii = ii + 1
                      if (xp(l,1) .gt. qx(ii)) go to 27 !find i on non-uniform 
                      ii = ii-1
                      ijkp(l,1)= ii


                      jj=0
 17                   continue
                      jj = jj + 1
                      if (xp(l,2) .gt. qy(jj)) go to 17 !find i on non-uniform 
                      jj = jj-1
                      ijkp(l,2)= jj                     

                      kk=2
 16                   continue
                      kk = kk + 1
                      if (xp(l,3) .gt. qz(kk)) go to 16 !find k on non-uniform 
                      kk = kk-1
                      ijkp(l,3)= kk
                      
                      mrat(l) = 1.0/m_pu
                      beta_p(l) = pu_beta_p
                      tags(l) = pu_tag
                      Ni_tot = l
                      do m=1,3
                         vp1(l,m) = vp(l,m)
                         ! Add the kinetic energy of the particle
                         ! (often zero)
                         input_E = input_E + 
     x                          0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2/
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
            call remove_ion(xp,vp,vp1,l)
            
         endif
 60   continue

      return
      end SUBROUTINE photoionization
c----------------------------------------------------------------------

      end MODULE chem_rates
