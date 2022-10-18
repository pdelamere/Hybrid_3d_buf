
      MODULE chem_rates

      USE global
      USE inputs
      USE gutsp_dd
      USE grid_interp
      
      contains

c----------------------------------------------------------------------
      function equilibrium_density(t, r) result(nn)
          ! gaussian barium release with equilibrium model of
          ! photoionization
          real, intent(in) :: t
          real, intent(in) :: r
          real :: nn
          real N
          real A

          ! N is the number of neutrals available to be ionized
          N = N0*exp(-t/tau_photo)

          ! We assume that the neutrals are distributed as
          ! a gaussian.
          A = N/((sqrt(pi)*vth_n*t)**3)
          nn = A*exp(-r**2/((vth_n*t)**2))
      end function equilibrium_density

      function excitation_density(t, r) result(nn)
          ! gaussian barium release with excitation model of
          ! photoionization
          real, intent(in) :: t
          real, intent(in) :: r
          real :: nn
          real N
          real A
          integer i

          ! This loop computes the number of ions at time t in the
          ! excitation model
          N = 0
          do i=1,4
            N = N + N0*vc(i)*exp(lambda(i)*t)
          enddo
          ! The number of neutrals is N0 - Ni
          N = N0 - N

          ! We assume that the neutrals are distributed as
          ! a gaussian.
          A = N/((sqrt(pi)*vth_n*t)**3)
          nn = A*exp(-r**2/((vth_n*t)**2))
      end function excitation_density

      function equilibrium_ionization_rate(t, N) result(r)
          real :: r
          real :: t
          real, optional :: N
          real :: NN
          if(present(N)) then
              NN = N
          else 
              NN = N0
          endif

          r = NN/tau_photo * exp(-t/tau_photo)
      end function equilibrium_ionization_rate

      function excitation_ionization_rate(t, N) result(r)
          real :: r
          real :: t
          integer :: i
          real, optional :: N
          real :: NN
          if(present(N)) then
              NN = N
          else 
              NN = N0
          endif

          r = 0
          do i=1, 4
            r = r + N0*vc(i)*lambda(i)*exp(lambda(i)*t)
          enddo
      end function excitation_ionization_rate

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
      real r

      call NC_coords(xx,yy,zz, x,y,z)

      r = sqrt(xx**2 + yy**2 + zz**2)
      neutral_density_continuous = excitation_density(simulated_time, r)

      return
      end FUNCTION neutral_density_continuous
c----------------------------------------------------------------------

      SUBROUTINE ionization(np,xp,vp,vp1)
          real np(nx,ny,nz)
          real xp(Ni_max,3)
          real vp(Ni_max,3)
          real vp1(Ni_max,3)
          !call photoionization(np,xp,vp,vp1)
          call photoionization_switching(np,xp,vp,vp1)
          !call fake_charge_exchange(xp, vp, vp1)
          call charge_exchange_ionization(xp,vp,vp1)
      end SUBROUTINE ionization

      subroutine photoionization_switching(np,xp,vp,vp1)
      real np(nx,ny,nz)
      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real sigma
      real dsigma
      dsigma = 1.27 ! km/s Don's expansion rate
      sigma = simulated_time*dsigma

      if((mod(procnum, 2) .eq. 0) .or. (3*sigma .gt. qz(nz)/2)) then
          call photoionization(np,xp,vp,vp1)
      else
          if(my_rank .eq. procnum/2) then
              call photoionization2(np,xp,vp,vp1)
          endif
      endif
      end subroutine
c----------------------------------------------------------------------
      SUBROUTINE photoionization(np,xp,vp,vp1)
c Ionizes the neutral cloud and fills particle
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
      real N

      call Neut_Center(cx,cy,cz)
      t = simulated_time

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
             
               pu_beta_p = 50*b_sw_thermal_H 
               pu_tag = pluto_photoionize_CH4_tag

               N = vol*neutral_density(i,j,k)
               new_micro = excitation_ionization_rate(t, N)*dt
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
                      vp(l,1) = x/simulated_time
                      vp(l,2) = y/simulated_time
                      vp(l,3) = z/simulated_time
                      
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
      return
      end SUBROUTINE photoionization
c----------------------------------------------------------------------

      subroutine photoionization2(np,xp,vp,vp1)
      ! Photoionization2 does exact normal distribution sampling. It
      ! should only be called by the one proc that contains the center
      ! of the neutral cloud and only before the edge of the cloud
      ! reaches the boundary.
      real np(nx,ny,nz)
      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      real t
      real sigma
      real dsigma
      real new_micro   
      integer new_macro   
      real pu_beta_p
      integer ierr
      integer l,m
      real r

      real cx,cy,cz         !x,y,z coord of neutral cloud center
      call Neut_Center2(cx,cy,cz) ! cz is in local coordinates

      dsigma = 1.27 ! km/s Don's expansion rate
      t = simulated_time
      sigma = t*dsigma

      pu_beta_p = 50*b_sw_thermal_H 
      new_micro = excitation_ionization_rate(t)*dt
      new_macro = nint(new_micro*beta*pu_beta_p)

      if((Ni_tot + new_macro) .gt. Ni_max) then
          call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
          stop
      endif

      ! This part is a little weird,
      ! randn_fill2 just fills xp with normally distributed random
      ! numbers with mu=0 and sigma=sigma. That's the coorect thing to
      ! do if xp's coordinates were centered on the cloud, but they
      ! aren't. We shift a little later (in the loop). First, vp and vp1
      ! are set. vp has a simple form when the position is in
      ! coordinates centered on the cloud (that's why we didn't shift it
      ! yet). vp1 = vp since the assumption is that the particle was
      ! always moving with this velocity.
      call randn_fill2(xp(Ni_tot+1:Ni_tot+new_macro,:), 0.0, sigma)
      vp(Ni_tot+1:Ni_tot+new_macro,:)=xp(Ni_tot+1:Ni_tot+new_macro,:)/t
      vp1(Ni_tot+1:Ni_tot+new_macro,:)=vp(Ni_tot+1:Ni_tot+new_macro,:)

      do l=Ni_tot+1, Ni_tot+new_macro
        do 
          ! test if the particle is in the domain (assume x and y are in
          ! the domain)
            if(abs(xp(l,3)) .lt. qz(nz)/2) then
                exit
            else
                call randn_one(xp(l,3))
                vp(l,3) = xp(l,3)/t
                vp1(l,3) = vp(l,3)
            endif
        enddo
        xp(l,1) = xp(l,1) + cx
        xp(l,2) = xp(l,2) + cy
        xp(l,3) = xp(l,3) + cz

        ijkp(l,1) = search(xp(l,1), qx)
        ijkp(l,2) = search(xp(l,2), qy)
        ijkp(l,3) = searchk(xp(l,3), qz)

        mrat(l) = ion_amu/m_pu
        beta_p(l) = pu_beta_p
        tags(l) = pluto_photoionize_CH4_tag
        do m=1,3
           ! Add the kinetic energy of the particle
           input_E = input_E + 
     x               0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x                     (beta*beta_p(l))
        enddo                     
      enddo

      Ni_tot = Ni_tot + new_macro
      end subroutine photoionization2


      SUBROUTINE charge_exchange_ionization(xp,vp,vp1)
          real xp(Ni_max,3)
          real vp(Ni_max,3)
          real vp1(Ni_max,3)
          integer l

          do l = 1, Ni_tot
          if((tags(l) .eq. pluto_photoionize_CH4_tag) .or.
     x       (tags(l) .eq. pluto_chex_CH4_tag) .or.
     x       (tags(l) .eq. sw_thermal_H_tag)) then
             call Ba_chex(l, xp, vp, vp1)
          endif

          enddo

      end SUBROUTINE charge_exchange_ionization

      SUBROUTINE Ba_chex(l, xp,vp,vp1)
         ! Handles both Ba + Ba^+ -> Ba^+ + Ba
         ! and Ba + O^+ -> Ba^+ + O
         ! I.e. these are charge exchange reactions that consume neutral
         ! barium and produce ionized barium. The only difference
         ! between the two reactions is the cross section and the
         ! identity of the ion being converted.
         integer l
         real xp(Ni_max,3)
         real vp(Ni_max,3)
         real vp1(Ni_max,3)

         real vrel
         real nn
         real x,y,z
         real nvx,nvy,nvz
         real rvx,rvy,rvz 
         real sigma
         real chex_inv_tau
         real chex_prob

         integer m
         nn = neutral_density_continuous(xp(l,1), xp(l,2), xp(l,3))

         call NC_coords(x,y,z, xp(l,1), xp(l,2), xp(l,3))
         nvx = x/simulated_time
         nvy = y/simulated_time
         nvz = z/simulated_time

         rvx = vp(l,1) - nvx
         rvy = vp(l,2) - nvy
         rvz = vp(l,3) - nvz

         vrel = sqrt(rvx**2 + rvy**2 + rvz**2)

         if((tags(l) .eq. pluto_photoionize_CH4_tag) .or.
     x      (tags(l) .eq.  pluto_chex_CH4_tag)) then
            sigma = sigma_Ba_res_chex
         elseif(tags(l) .eq. sw_thermal_H_tag) then
            sigma = sigma_Ba_chex
         endif
         chex_inv_tau = nn*sigma*vrel
         chex_prob = dt*chex_inv_tau

         if (pad_ranf() .lt. chex_prob) then
             ! Remove kinetic energy of the old particle
             do m=1,3
                input_E = input_E - 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              beta*beta_p(l) 
             enddo
             vp(l,1) = nvx
             vp(l,2) = nvy
             vp(l,3) = nvz
             vp1(l,:) = vp(l,:)
             mrat(l) = ion_amu/m_pu
             tags(l) = pluto_chex_CH4_tag
             ! Add kinetic energy of the new particle
             do m=1,3
                input_E = input_E + 
     x              0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x              beta*beta_p(l) 
             enddo
         endif
      end SUBROUTINE Ba_chex 

      end MODULE chem_rates
