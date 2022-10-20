
      MODULE chem_rates

      USE global
      USE inputs
      USE gutsp_dd
      USE grid_interp
      
      contains

c----------------------------------------------------------------------
      function gaussian_density(N, t, r) result(nn)
          ! Assuming an expasion rate given by vth_n
          real, intent(in) :: N, t, r
          real :: nn
          real :: A

          A = N/((sqrt(pi)*vth_n*t)**3)
          nn = A*exp(-r**2/((vth_n*t)**2))
      end function gaussian_density

      function equilibrium_N(t) result(N)
          ! total number of neutrals at time t in the equilibrium model
          real, intent(in) :: t
          real :: N

          N = N0*exp(-t/tau_photo)
      end function equilibrium_N

      function excitation_N(t) result(N)
          ! total number of neutrals at time t in the excitation model
          real, intent(in) :: t
          real :: N
          integer i

          ! This loop computes the number of *ions* at time t in the
          ! excitation model
          N = 0
          do i=1,4
            N = N + N0*vc(i)*exp(lambda(i)*t)
          enddo
          ! The number of neutrals is N0 - Ni
          N = N0 - N
      end function excitation_N

      function equilibrium_ionization_rate(t) result(r)
          ! Ionization rate of neutral micro particles 
          ! in the equilibrium model. Number of new ions per second
          ! is N_n * rate
          real :: r
          ! t is unused (argument is kept for similarity with 
          ! excitation model
          real :: t 

          r = 1/tau_photo
      end function equilibrium_ionization_rate

      function excitation_ionization_rate(t) result(r)
          ! Ionization rate of neutral micro particles 
          ! in the excitation model. Number of new ions per second
          ! is N_n * rate
          real :: r
          real :: t
          integer :: i

          r = 0
          do i=1, 4
            r = r + vc(i)*lambda(i)*exp(lambda(i)*t)
          enddo
      end function excitation_ionization_rate

      !!!!!!!!!!! Make sure these are both using the same model e.g.
      ! both excitation or both equilibrium
      function ionization_rate(t) result(r)
          real t
          real r
          r = equilibrium_ionization_rate(t)
      end function
      function neutral_N(t) result(N)
          real t
          real N
          N = equilibrium_N(t)
      end function

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
      real N,t,r
      t = simulated_time

      call NC_coords(xx,yy,zz, x,y,z)
      r = sqrt(xx**2 + yy**2 + zz**2)

      N = neutral_N(t)
      neutral_density_continuous=gaussian_density(N,t,r)

      return
      end FUNCTION neutral_density_continuous
c----------------------------------------------------------------------

      SUBROUTINE ionization(np,xp,vp,vp1)
          real np(nx,ny,nz)
          real xp(Ni_max,3)
          real vp(Ni_max,3)
          real vp1(Ni_max,3)
          call photoionization(np,xp,vp,vp1)
          !call photoionization_switching(np,xp,vp,vp1)
          !call charge_exchange_ionization(xp,vp,vp1)
      end SUBROUTINE ionization

      subroutine photoionization_switching(np,xp,vp,vp1)
      real np(nx,ny,nz)
      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real sigma
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
      real new_micro
      integer new_macro 
      
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
      integer ierr

      call Neut_Center(cx,cy,cz)
      t = simulated_time

      vol = dx**3
      l1 = Ni_tot+1
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               pu_beta_p = 5*b_sw_thermal_H 
               pu_tag = pluto_photoionize_CH4_tag

               N = vol*neutral_density(i,j,k)
               new_micro = N * ionization_rate(t) * dt
               new_macro = nint(new_micro*pu_beta_p)

               if(Ni_tot + new_macro .gt. Ni_max) then
                   call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
                   stop
               endif

               do ll = 1,new_macro
                  l = Ni_tot + 1

                  xp(l,1) = qx(i) + (pad_ranf()-0.5)*dx_grid(i)
                  xp(l,2) = qy(j) + (pad_ranf()-0.5)*dy_grid(j)
                  xp(l,3) = qz(k) + (pad_ranf()-0.5)*dz_grid(k)

                  call NC_coords(x,y,z, xp(l,1), xp(l,2), xp(l,3))

                  vp(l,1) = x/simulated_time
                  vp(l,2) = y/simulated_time
                  vp(l,3) = z/simulated_time

                  ijkp(l,1) = search(xp(l,1), qx)
                  ijkp(l,2) = search(xp(l,2), qy)
                  ijkp(l,3) = searchk(xp(l,3), qz)

                  mrat(l) = ion_amu/m_pu
                  beta_p(l) = pu_beta_p
                  tags(l) = pu_tag
                  Ni_tot = l
                  do m=1,3
                     vp1(l,m) = vp(l,m)
                     ! Add the kinetic energy of the particle
                     input_E = input_E + 
     x               0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
     x                     (beta*beta_p(l))
                  enddo                     
               enddo
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
      real new_micro   
      integer new_macro   
      real pu_beta_p
      integer ierr
      integer l,m
      real r
      real N

      real cx,cy,cz         !x,y,z coord of neutral cloud center
      call Neut_Center2(cx,cy,cz) ! cz is in local coordinates

      t = simulated_time
      sigma = t*dsigma

      pu_beta_p = 5*b_sw_thermal_H 

      N = neutral_N(t)
      new_micro = N * ionization_rate(t)*dt
      new_macro = nint(new_micro*pu_beta_p)

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
         real chex_rate
         real chex_num
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
         chex_rate = nn*sigma*vrel ! * beta_p(l)
         chex_num = dt*chex_rate
         ! Chex_prob is not exactly a probablility. Chex prob is given
         ! by the predicted number of charge exchange events this
         ! timestep divided by how many micro particles are represented
         ! by this macro particle. Then we "charge exchange" the whole
         ! macro particle if rand() < chex_prob.
         ! You'd think beta_p would make an appearance for the number of
         ! micro particles, but since it contributes to the expected
         ! number of charge exchanges also it ends up canceling out.
         chex_prob = chex_num ! / beta_p(l)

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
