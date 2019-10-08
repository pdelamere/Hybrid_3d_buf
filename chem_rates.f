
      MODULE chem_rates

      USE global
      USE inputs
      USE gutsp_dd
      USE grid_interp
      
      contains

c----------------------------------------------------------------------
      real function atmosphere(r)
          real, intent(in) :: r
          !! gaussian barium release
          real N0, N
          real vth_n
          real A
          real t
          real tstart
          tstart = 0.2
          N0 = 2.631e24
          vth_n = 0.5
          t = simulated_time + tstart
          N = N0*exp(-t/tau_photo)
          A = N/((sqrt(pi)*vth_n*t)**3)

          atmosphere = A*exp(-r**2/((vth_n*t)**2))
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
      real r

      call Neut_Center(cx,cy,cz)
      xx = x - cx
      yy = y - cy
      ! zz is computed differently since we need to convert z (which is
      ! local) to a global z.
      zz = (z + (procnum-(cart_rank+1))*qz(nz-1)) - cz

      r = sqrt(xx**2 + yy**2 + zz**2)
      neutral_density_continuous = atmosphere(r)

      return
      end FUNCTION neutral_density_continuous
c----------------------------------------------------------------------

      SUBROUTINE ionization(np,xp,vp,vp1)
          real np(nx,ny,nz)
          real xp(Ni_max,3)
          real vp(Ni_max,3)
          real vp1(Ni_max,3)
          !call charge_exchange_ionization(xp,vp,vp1)
          call photoionization(np,xp,vp,vp1)
      end SUBROUTINE ionization

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

      real x,y,z

      real pu_beta_p
      real pu_tag
      integer i,j,k,l,m
      integer ii,jj,kk,ll



      call Neut_Center(cx,cy,cz)

c get source density

      
      vol = dx**3
      l1 = Ni_tot+1
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1

               x = qx(i)-cx
               y = qy(j)-cy
               z = gz(k)-cz ! global z
               r = sqrt(x**2+y**2+z**2)
             
               pu_beta_p = b_sw_thermal_H 
               pu_tag = pluto_photoionize_CH4_tag

               cur_micro = np(i,j,k)*vol
               new_micro = vol*neutral_density(i,j,k)*dt/tau_photo
               new_macro = new_micro*beta*pu_beta_p

               do ll = 1,min(nint(new_macro), Ni_max - Ni_tot)
                  l = Ni_tot + 1
                  vp(l,1) = x/simulated_time
                  vp(l,2) = y/simulated_time
                  vp(l,3) = z/simulated_time

                  xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                  xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                  xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)

                  !find i on non-uniform 
                  ii=0
 26               continue
                  ii = ii + 1
                  if (xp(l,1) .gt. qx(ii)) go to 26
                  ii = ii-1
                  ijkp(l,1)= ii

                  !find j on non-uniform 
                  jj=0
 18               continue
                  jj = jj + 1
                  if (xp(l,2) .gt. qy(jj)) go to 18
                  jj = jj-1
                  ijkp(l,2)= jj

                  !find k on non-uniform 
                  kk=2
 15               continue
                  kk = kk + 1
                  if (xp(l,3) .gt. qz(kk)) go to 15
                  kk = kk-1
                  ijkp(l,3)= kk

                  mrat(l) = ion_amu/m_pu
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


                      !find i on non-uniform 
                      ii=0
 27                   continue
                      ii = ii + 1
                      if (xp(l,1) .gt. qx(ii)) go to 27
                      ii = ii-1
                      ijkp(l,1)= ii


                      !find i on non-uniform 
                      jj=0
 17                   continue
                      jj = jj + 1
                      if (xp(l,2) .gt. qy(jj)) go to 17
                      jj = jj-1
                      ijkp(l,2)= jj                     

                      !find k on non-uniform 
                      kk=2
 16                   continue
                      kk = kk + 1
                      if (xp(l,3) .gt. qz(kk)) go to 16
                      kk = kk-1
                      ijkp(l,3)= kk
                      
                      mrat(l) = ion_amu/m_pu
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
