      MODULE gutsp_dd

      USE global
      USE misc
      USE mpi
      USE boundary
      USE grid_interp
      USE iso_fortran_env, only: error_unit

      contains

c----------------------------------------------------------------------
      SUBROUTINE remove_ion(xp,vp,vp1,ion_l,separate)
c Removes particles from simulation that have gone out of bounds
c----------------------------------------------------------------------


      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      integer ion_l
      integer, optional :: separate

      do 5 m=1,3   !remove ion energy from total input energy
         input_E = input_E
     x             -0.5*(mion/mrat(ion_l))*(vp(ion_l,m)*km_to_m)**2
     x             / (beta*beta_p(ion_l))
 5    continue

      if(present(separate))then
          write(*,*) 'Separate: removing ion...',ion_l
      else
          write(*,*) 'removing ion...',ion_l
      endif

      do 10 l=ion_l,Ni_tot-1
c         m_arr(l) = m_arr(l+1)
         mrat(l) = mrat(l+1)
         tags(l) = tags(l+1)
         beta_p(l) = beta_p(l+1)
         do 7 m=1,8
            wght(l,m) = wght(l+1,m)
 7       continue
         do 10 m=1,3 
            xp(l,m) = xp(l+1,m)
            vp(l,m) = vp(l+1,m)
            vp1(l,m) = vp1(l+1,m)
            ijkp(l,m) = ijkp(l+1,m)
c            wquad(l,m) = wquad(l+1,m)
 10      continue


      Ni_tot = Ni_tot - 1

      return
      end SUBROUTINE remove_ion
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE check_min_den(np,xp,vp,vp1,up,bt)
c----------------------------------------------------------------------

      real np(nx,ny,nz),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     up(nx,ny,nz,3),
     x     bt(nx,ny,nz,3)

      integer Ni_tot_in
      integer Ni_out, Ni_in
      integer source, dest
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)
      integer npart,ipart
      integer cnt

      real, dimension(:,:), allocatable :: in_part
      real, dimension(:,:), allocatable :: out_part
      integer, dimension(:,:), allocatable :: in_ijkp
      integer, dimension(:,:), allocatable :: out_ijkp
      real, dimension(:,:), allocatable :: in_part_wght
      real, dimension(:,:), allocatable :: out_part_wght
      real, dimension(:), allocatable :: in_mass
      real, dimension(:), allocatable :: out_mass
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:), allocatable :: out_tags
      real, dimension(:), allocatable :: in_tags
      real, dimension(:), allocatable :: vsqrd_out

      real minden !minimum wake density 
      real den_part ! density of 1 particle per cell
      real ak
      real btot,a1,a2,womega,phi,deltat

      den_part = 1/(beta*dx**3)

      minden = nf_init/10.0
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
               if (Ni_tot+1 .ge. Ni_max) then
                   write(*,*) "Warning: Ni_max ions. Minden failed"
               endif
               if (((np(i,j,k) .le. minden) .or. 
     x              (np(i,j,k) .le. 2.0*den_part)) .and.
     x              (Ni_tot + 1 .lt. Ni_max)) then
                  npart = nint(minden/(np(i,j,k)))
                  write(*,*) npart, "dummy particles added"
                  do ipart = 1,npart 
                     l=Ni_tot + 1 !beginning array element for new borns    
                  
                     vp(l,1) = up(i,j,k,1)
                     vp(l,2) = up(i,j,k,2)
                     vp(l,3) = up(i,j,k,3)
                     xp(l,1) = qx(i) + (0.5-pad_ranf())*dx_grid(i)
                     xp(l,2) = qy(j) + (0.5-pad_ranf())*dy_grid(j)
                     xp(l,3) = qz(k) + (0.5-pad_ranf())*dz_grid(k)
                     
                     ii=0
 16                  continue
                     ii = ii + 1
                     if (xp(l,1) .gt. qx(ii)) go to 16 !find i on non-uniform 
                     ii = ii-1
                     ijkp(l,1)= ii
                     
                     
                     jj=0
 18                  continue
                     jj = jj + 1
                     if (xp(l,2) .gt. qy(jj)) go to 18 !find j on non-uniform 
                     jj = jj-1
                     ijkp(l,2)= jj
                     
                     
                     kk=0
 15                  continue
                     kk = kk + 1
                     if (xp(l,3) .gt. qz(kk)) go to 15 !find k on non-uniform 
                     kk = kk-1
                     ijkp(l,3)= kk



                     mrat(l) = 1.0
                     beta_p(l) = 1.0
                     tags(l) = dummy_particle_tag
                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            enddo
         enddo
      enddo
c ---------z exchange down---------------------------------------

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      where (xp(1:Ni_tot,3) .le. qz(2))
         in_bounds(1:Ni_tot)= .false.
         ijkp(1:Ni_tot,3) = nz
         xp(1:Ni_tot,3) = qz(nz)-(qz(2)-xp(1:Ni_tot,3))
      endwhere

      Ni_tot_in = count(in_bounds(1:Ni_tot))

      Ni_out = count(.not.in_bounds(1:Ni_tot))
      allocate(out_part(Ni_out,3))
      allocate(out_ijkp(Ni_out,3))
      allocate(out_part_wght(Ni_out,8))
      allocate(out_mass(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(out_tags(Ni_out))
      allocate(vsqrd_out(Ni_out))

      dest = down_proc
      source = up_proc
      
      call MPI_ISEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)


      call MPI_WAITALL(2, reqs, stats, ierr)

      allocate(in_part(Ni_in,3))
      allocate(in_ijkp(Ni_in,3))
      allocate(in_part_wght(Ni_in,8))
      allocate(in_mass(Ni_in))
      allocate(in_tags(Ni_in))


      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(xp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(xp, in_bounds,3)

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      xp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp, in_bounds, 3)

      vsqrd_out(:) = out_part(:,1)**2 + 
     x               out_part(:,2)**2 + out_part(:,3)**2
      vsqrd_out(:) = vsqrd_out*(km_to_m)**2

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      vp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp1(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp1, in_bounds, 3)

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      vp1(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_ijkp(1:Ni_out,m) = 
     x          pack(ijkp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      cnt = 1
      do l = 1,Ni_tot
        if (in_bounds(l)) then
           ijkp(cnt,:) = ijkp(l,:)
           cnt = cnt+1
        endif
      enddo


      call MPI_ISEND(out_ijkp, 3*Ni_out, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_ijkp, 3*Ni_in, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      ijkp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_ijkp(:,:)


      do m = 1,8
         out_part_wght(1:Ni_out,m) = 
     x          pack(wght(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(wght, in_bounds, 8)

      call MPI_ISEND(out_part_wght, 8*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part_wght, 8*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      wght(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part_wght(:,:)


      out_beta_p(1:Ni_out) = 
     x     pack(beta_p(1:Ni_tot), .not.in_bounds(1:Ni_tot))

      call pack_pd(beta_p, in_bounds, 1)
      
      call MPI_ISEND(out_beta_p, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      beta_p(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      out_tags(1:Ni_out) =
     x     pack(tags(1:Ni_tot), .not.in_bounds(1:Ni_tot))

      call pack_pd(tags, in_bounds, 1)

      call MPI_ISEND(out_tags, Ni_out, MPI_REAL, dest, tag,
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_tags, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)

      tags(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_tags(:)


      out_mass(1:Ni_out) = 
     x     pack(mrat(1:Ni_tot), .not.in_bounds(1:Ni_tot))

      call pack_pd(mrat, in_bounds, 1)


      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      mrat(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      Ni_tot = Ni_tot - Ni_out + Ni_in

      ! remove energy of outgoing particles
      input_E = input_E - sum(0.5*(mion/out_mass(:))*vsqrd_out(:)/
     x                       (beta*out_beta_p(:)))

      ! add energy of incomming particles
      do m = 1,3
         input_E = input_E + sum(0.5*(mion/mrat(Ni_tot_in+1:Ni_tot))*
     x               (vp(Ni_tot_in+1:Ni_tot,m)*km_to_m)**2 
     x                /(beta*beta_p(Ni_tot_in+1:Ni_tot)))
      enddo

      deallocate(out_part)
      deallocate(out_ijkp)
      deallocate(out_part_wght)
      deallocate(out_mass)
      deallocate(out_beta_p)
      deallocate(out_tags)
      deallocate(vsqrd_out)

      deallocate(in_part)
      deallocate(in_ijkp)
      deallocate(in_part_wght)
      deallocate(in_mass)
      deallocate(in_tags)

      return
      end SUBROUTINE check_min_den
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE extrapol_up(up,vp,vp1,np)
c This subroutine does the provisional extrapolation of the particle
c bulk flow velocity to time level n, and replaces up_n-3/2 
c with up_n-1/2
c----------------------------------------------------------------------
CVD$R VECTOR
c      include 'incurv.h'

      real up(nx,ny,nz,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     np(nx,ny,nz)
      
      real v_at_n(Ni_max,3)

      do 10 m=1,3
         do 10 l=1,Ni_tot
            v_at_n(l,m) = 1.5*vp(l,m) - 0.5*vp1(l,m)
 10      continue

      call update_up(v_at_n,np,up)

      return
      end SUBROUTINE extrapol_up
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE get_Ep(Ep,aj,np,up,btc,nu)
c----------------------------------------------------------------------
CVD$F VECTOR
c      include 'incurv.h'

      real Ep(Ni_max,3),
     x     aj(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
     x     btc(nx,ny,nz,3),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)
      real us(ny,nz,3)

      real ajc(nx,ny,nz,3),       !aj at cell center
     x     upc(nx,ny,nz,3),       !up at cell center
     x     ufc(nx,ny,nz,3)       !uf at cell center
c     x     gradPc(nx,ny,nz,3)     !gradP at cell center

      real aa(3), bb(3), cc(3)    !dummy vars for doing cross product
      real ntot                   !total plasma density
      real fnf,fnp                !fraction nf,np of total n
      real aj3(3),up3(3),uf3(3),  !vector field values at Ba position
     x     btc3(3),gradP3(3)
      real np_at_Ba               !particle density at particle pos
      real nf_at_Ba

      us = 0.0
      call face_to_center(aj,ajc,us)! only if upstream B is const

      do 10 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

         do 15 m=1,3 
            aj3(m) = ajc(i,j,k,m)*wght(l,1) + ajc(ip,j,k,m)*wght(l,2) 
     x          + ajc(i,j,kp,m)*wght(l,3) + ajc(ip,j,kp,m)*wght(l,4)
     x          + ajc(i,jp,k,m)*wght(l,5) + ajc(ip,jp,k,m)*wght(l,6)
     x          + ajc(i,jp,kp,m)*wght(l,7) + ajc(ip,jp,kp,m)*wght(l,8)

            up3(m) = up(i,j,k,m)*wght(l,1) + up(ip,j,k,m)*wght(l,2) 
     x          + up(i,j,kp,m)*wght(l,3) + up(ip,j,kp,m)*wght(l,4)
     x          + up(i,jp,k,m)*wght(l,5) + up(ip,jp,k,m)*wght(l,6)
     x          + up(i,jp,kp,m)*wght(l,7) + up(ip,jp,kp,m)*wght(l,8)

            btc3(m) = btc(i,j,k,m)*wght(l,1) 
     x               + btc(ip,j,k,m)*wght(l,2) 
     x               + btc(i,j,kp,m)*wght(l,3) 
     x               + btc(ip,j,kp,m)*wght(l,4)
     x               + btc(i,jp,k,m)*wght(l,5) 
     x               + btc(ip,jp,k,m)*wght(l,6)
     x               + btc(i,jp,kp,m)*wght(l,7) 
     x               + btc(ip,jp,kp,m)*wght(l,8)

 15         continue

         do 20 m=1,3
            aa(m) = aj3(m) - up3(m)
            bb(m) = btc3(m)                   

 20         continue

         cc(1) = aa(2)*bb(3) - aa(3)*bb(2)    !do cross product
         cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
         cc(3) = aa(1)*bb(2) - aa(2)*bb(1)

         do 30 m=1,3
            Ep(l,m) = cc(m) 
            Ep(l,m) = Ep(l,m)*mrat(l) !O_to_Ba

 30         continue


 10      continue

      return
      end SUBROUTINE get_Ep
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vplus_vminus(Ep,btc,vp,vplus,vminus)
c----------------------------------------------------------------------

      real Ep(Ni_max,3),
     x     btc(nx,ny,nz,3),   !bt at cell center
     x     vp(Ni_max,3),      !particle velocities at t level n-1/2
     x     vplus(Ni_max,3),
     x     vminus(Ni_max,3)

      real a1, a2, a3,       !coefficients for calculating vplus
     x     a_d,              !denominator for a1,a2,a3
     x     B2,dt2,           !B*B,dt*dt
     x     Bx,By,Bz,         !B for cross product call
     x     vminus_x_B(3),    !v- x B
     x     vminus_dot_B     !v- . B

      real btc3(3)
      

      do 10 m=1,3
         do 10 l=1,Ni_tot 
            vminus(l,m) = vp(l,m) + 0.5*dt*Ep(l,m)
 10         continue

      do 30 l=1,Ni_tot 

         i = ijkp(l,1)
         j = ijkp(l,2)
         k = ijkp(l,3)
   
         ip = i+1
         jp = j+1
         kp = k+1

         do 35 m=1,3
            btc3(m) = btc(i,j,k,m)*wght(l,1) 
     x           + btc(ip,j,k,m)*wght(l,2) 
     x           + btc(i,j,kp,m)*wght(l,3) 
     x           + btc(ip,j,kp,m)*wght(l,4)
     x           + btc(i,jp,k,m)*wght(l,5) 
     x           + btc(ip,jp,k,m)*wght(l,6)
     x           + btc(i,jp,kp,m)*wght(l,7) 
     x           + btc(ip,jp,kp,m)*wght(l,8)
            
 35      continue

         vminus_x_B(1) = vminus(l,2)*btc3(3)*mrat(l) -
     x                     vminus(l,3)*btc3(2)*mrat(l)
         vminus_x_B(2) = vminus(l,3)*btc3(1)*mrat(l) -
     x                     vminus(l,1)*btc3(3)*mrat(l)
         vminus_x_B(3) = vminus(l,1)*btc3(2)*mrat(l) -
     x                     vminus(l,2)*btc3(1)*mrat(l)

         vminus_dot_B = vminus(l,1)*btc3(1)*mrat(l) +
     x                     vminus(l,2)*btc3(2)*mrat(l) +
     x                     vminus(l,3)*btc3(3)*mrat(l)

         Bx = btc3(1)*mrat(l)
         By = btc3(2)*mrat(l)
         Bz = btc3(3)*mrat(l)
      
         B2 = Bx**2 + By**2 + Bz**2
         dt2 = dt**2

         a_d = 1 + (B2*dt2/4.0)
         a1 = (1 - (B2*dt2/4.0)) / a_d
         a2 = dt / a_d
         a3 = 0.5*dt2 / a_d

         do 40 m=1,3
            vplus(l,m) = a1*vminus(l,m) + a2*vminus_x_B(m) + 
     x           a3*vminus_dot_B*btc3(m)*mrat(l)
 40      continue

 30   continue


      return
      end SUBROUTINE get_vplus_vminus
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE improve_up(vp1,vplus,vminus,up,np)
c The routine calculates v at time level n, and the associated bulk
c flow velocity up using the v+, v- technique.  The new up at
c time level n replaces the provisional extrapolation for up.
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vp1(Ni_max,3),    !particle velocities at t level n
     x     vplus(Ni_max,3),
     x     vminus(Ni_max,3),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz)

      do 10 m=1,3
         do 10 l = 1,Ni_tot
            vp1(l,m) = 0.5*(vplus(l,m) + vminus(l,m))
c            write(*,*) 'vp1....',m,vp1(l,m)
 10         continue

      call update_up(vp1,np,up)

      return
      end SUBROUTINE improve_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vp_final(Ep,vp,vp1,vplus)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real Ep(Ni_max,3),
     x     vp(Ni_max,3),    !particle velocities at t level n+1/2
     x     vp1(Ni_max,3),   !particle velocity at t level n-1/2
     x     vplus(Ni_max,3)

      do 10 m=1,3
         do 10 l = 1,Ni_tot
            vp1(l,m) = vp(l,m)  !to be used in extrapol_up for n-3/2
            vp(l,m) = vplus(l,m) + 0.5*dt*Ep(l,m)

 10         continue

      return
      end SUBROUTINE get_vp_final
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE exchange_ion_in(xp,vp,vp1,xp_buf,vp_buf) 
c----------------------------------------------------------------------
c Exchange ions from the main domain into the inflow buffer

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)
      real xp_buf(Ni_max_buf,3)
      real vp_buf(Ni_max_buf,3)

      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
      real, dimension(:), allocatable :: out_mrat
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:), allocatable :: out_tags
      integer, dimension(:,:), allocatable :: out_ijkp
      integer Ni_tot_in, Ni_out
      integer :: cnt


      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.
      
      where (xp(1:Ni_tot,1) .gt. qx(nx))
         in_bounds(1:Ni_tot) = .false.
      endwhere
      
      Ni_tot_in = count(in_bounds)
      Ni_out = count(.not.in_bounds(1:Ni_tot))
      
      allocate(out_xp(Ni_out,3))
      allocate(out_vp(Ni_out,3))
      allocate(out_mrat(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(out_tags(Ni_out))
      
      do m = 1,3 
         
         out_xp(1:Ni_out,m) = pack(xp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         out_vp(1:Ni_out,m) = pack(vp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
      enddo
 
      call pack_pd(xp, in_bounds, 3)
      call pack_pd(vp, in_bounds, 3)
      call pack_pd(vp1, in_bounds, 3)      

      cnt = 1
      do l = 1,Ni_tot
        if (in_bounds(l)) then
           ijkp(cnt,:) = ijkp(l,:)
           cnt = cnt+1
        endif
      enddo

     
      call pack_pd(wght, in_bounds, 8)
      
      out_mrat(1:Ni_out) = pack(mrat(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      out_beta_p(1:Ni_out) = pack(beta_p(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      out_tags(1:Ni_out) = pack(tags(1:Ni_tot),
     x     .not.in_bounds(1:Ni_tot))
      
      call pack_pd(mrat, in_bounds, 1)
      call pack_pd(beta_p, in_bounds, 1)
      call pack_pd(tags, in_bounds, 1)
            
      do m = 1,3
         xp_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out,m) = out_xp(:,m)
         vp_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out,m) = out_vp(:,m)
      enddo
      
      mrat_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out) = out_mrat(:)
      beta_p_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out) = out_beta_p(:)
      tags_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out) = out_tags(:)
      
c     remove energy

      do l = 1,Ni_out  
         do m=1,3
            input_E = input_E - 
     x           0.5*(mion/out_mrat(l))*(out_vp(l,m)*km_to_m)**2 /
     x           (beta*out_beta_p(l))
         enddo         
      enddo

      Ni_tot_buf = Ni_tot_buf + Ni_out
      Ni_tot = count(in_bounds)
      
      
      deallocate(out_xp)
      deallocate(out_vp)
      deallocate(out_mrat)
      deallocate(out_beta_p)
      deallocate(out_tags)
      

      return
      end SUBROUTINE exchange_ion_in
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE exchange_ion_out(xp,vp,vp1,xp_buf,vp_buf,
     x     E,Bt,xp_out_buf,vp_out_buf,E_out_buf,
     x     B_out_buf,mrat_out_buf, save_unit) 
c----------------------------------------------------------------------
c Exchange ions from the main domain to the outflow buffer (delete them)

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)
      real xp_buf(Ni_max_buf,3)
      real vp_buf(Ni_max_buf,3)
      real E(nx,ny,nz,3)
      real Bt(nx,ny,nz,3)
      real xp_out_buf(Ni_max_buf,3)
      real vp_out_buf(Ni_max_buf,3)
      real E_out_buf(Ni_max_buf,3)
      real B_out_buf(Ni_max_buf,3)
      real mrat_out_buf(Ni_max_buf)

      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
      real, dimension(:), allocatable :: out_mrat
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:), allocatable :: out_tags
      real, dimension(:,:), allocatable :: out_E
      real, dimension(:,:), allocatable :: out_B
      integer, dimension(:,:), allocatable :: out_ijkp
      integer Ni_tot_in, Ni_out
      integer :: cnt
      integer :: save_unit


      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      
      where(xp(1:Ni_tot,1) .le. qx(1))
           in_bounds(1:Ni_tot) = .false.
      endwhere         


      Ni_tot_in = count(in_bounds)
      Ni_out = count(.not.in_bounds(1:Ni_tot))
      
      allocate(out_xp(Ni_out,3))
      allocate(out_vp(Ni_out,3))
      allocate(out_E(Ni_out,3))
      allocate(out_B(Ni_out,3))
      allocate(out_ijkp(Ni_out,3))
      allocate(out_mrat(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(out_tags(Ni_out))
      

      do m = 1,3 
         out_xp(1:Ni_out,m) = pack(xp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         out_vp(1:Ni_out,m) = pack(vp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         out_ijkp(1:Ni_out,m) = pack(ijkp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         
      enddo

      call pack_pd(xp, in_bounds, 3)
      call pack_pd(vp, in_bounds, 3)
      call pack_pd(vp1, in_bounds, 3)

      cnt = 1
      do l = 1,Ni_tot
        if (in_bounds(l)) then
           ijkp(cnt,:) = ijkp(l,:)
           cnt = cnt+1
        endif
      enddo
      
      do l = 1,Ni_out 
         out_E(l,:) = E(out_ijkp(l,1),out_ijkp(l,2),out_ijkp(l,3),:)
         out_B(l,:) = Bt(out_ijkp(l,1),out_ijkp(l,2),out_ijkp(l,3),:)
      enddo

 
      call pack_pd(wght, in_bounds, 8)
     
      out_mrat(1:Ni_out) = pack(mrat(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      out_beta_p(1:Ni_out) = pack(beta_p(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      out_tags(1:Ni_out) = pack(tags(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))

      ! Save heavy ions that are leaving
      write(save_unit) pack(out_mrat,   out_mrat<0.1)
      write(save_unit) pack(out_beta_p, out_mrat<0.1)
      write(save_unit) pack(out_tags,   out_mrat<0.1)
      
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      do m = 1,3
         xp_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out,m) = 
     x        out_xp(1:Ni_out,m)
         vp_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out,m) = 
     x        out_vp(1:Ni_out,m)
         E_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out,m) = 
     x        out_E(1:Ni_out,m)
         B_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out,m) = 
     x        out_B(1:Ni_out,m)
      enddo
      
      mrat_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out) = 
     x     out_mrat(:)
      
      call pack_pd(mrat, in_bounds, 1)
      call pack_pd(beta_p, in_bounds, 1)
      call pack_pd(tags, in_bounds, 1)
            
c     remove energy

      do l = 1,Ni_out  
         do m=1,3
            input_E = input_E - 
     x           0.5*(mion/out_mrat(l))*(out_vp(l,m)*km_to_m)**2 /
     x           (beta*out_beta_p(l))
         enddo         
      enddo
      
      deallocate(out_xp)
      deallocate(out_vp)
      deallocate(out_E)
      deallocate(out_B)
      deallocate(out_ijkp)
      deallocate(out_mrat)
      deallocate(out_beta_p)
      deallocate(out_tags)
      
      Ni_tot = count(in_bounds)

         
      return
      end SUBROUTINE exchange_ion_out
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE move_ion_half(xp,vp,vp1,Ep)
c----------------------------------------------------------------------

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)
      real Ep(Ni_max,3)
      integer zindx(1)
      integer Ni_tot_in
      integer Ni_out, Ni_in
      integer source, dest
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)
      integer stat(MPI_STATUS_SIZE)

      real dth                 !half time step

      real, dimension(:,:), allocatable :: in_part
      real, dimension(:,:), allocatable :: out_part
      integer, dimension(:,:), allocatable :: in_ijkp
      integer, dimension(:,:), allocatable :: out_ijkp
      real, dimension(:,:), allocatable :: in_part_wght
      real, dimension(:,:), allocatable :: out_part_wght
      real, dimension(:), allocatable :: in_mass
      real, dimension(:), allocatable :: out_mass
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:), allocatable :: out_tags
      real, dimension(:), allocatable :: in_tags
      real, dimension(:), allocatable :: vsqrd_out
      integer :: cnt


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      dth = dt/2.

      do 10 l=1,Ni_tot                   !make 1/2 time step advance

         xp0 = xp(l,1)
         vp0 = vp(l,1)
         xp(l,1) = xp(l,1) + dth*vp(l,1)

 


         xp(l,2) = xp(l,2) + dth*vp(l,2)

         xp(l,3) = xp(l,3) + dth*vp(l,3)
 10      continue



      where (xp(:,2) .gt. qy(ny-1))
         xp(:,2) = qy(1) + ( xp(:,2) - qy(ny-1) )
      endwhere

      where (xp(:,2) .le. qy(1)) 
         xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
      endwhere


c -------------------z exchange, up-----------------------------

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      where (xp(1:Ni_tot,3) .gt. qz(nz))
         in_bounds(1:Ni_tot)= .false.
         ijkp(1:Ni_tot,3) = 2
         xp(1:Ni_tot,3) = qz(2)+(xp(1:Ni_tot,3)-qz(nz))
      endwhere

      Ni_tot_in = count(in_bounds(1:Ni_tot))

      Ni_out = count(.not.in_bounds(1:Ni_tot))

      allocate(out_part(Ni_out,3))
      allocate(out_ijkp(Ni_out,3))
      allocate(out_part_wght(Ni_out,8))
      allocate(out_mass(Ni_out))
      allocate(vsqrd_out(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(out_tags(Ni_out))

      dest = up_proc
      source = down_proc

      call MPI_SEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, ierr)
      call MPI_RECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, stat, ierr)

      if(Ni_tot_in + Ni_in .gt. Ni_max) then
          write(error_unit,*) 'Not enough space to pass particles in1',
     x                    Ni_tot_in, Ni_in, my_rank
          stop 1
      endif
      
      allocate(in_part(Ni_in,3))
      allocate(in_ijkp(Ni_in,3))
      allocate(in_part_wght(Ni_in,8))
      allocate(in_mass(Ni_in))
      allocate(in_tags(Ni_in))


      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      do m = 1,3
         out_part(1:Ni_out,m) = 
     x        pack(xp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(xp, in_bounds, 3)
                

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

            
      xp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp, in_bounds, 3)

      vsqrd_out(:) = out_part(:,1)**2 + 
     x               out_part(:,2)**2 + out_part(:,3)**2
      vsqrd_out(:) = vsqrd_out*(km_to_m)**2


      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      
      vp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp1(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp1, in_bounds, 3)

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

       vp1(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_ijkp(1:Ni_out,m) = 
     x          pack(ijkp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo


      cnt = 1
      do l = 1,Ni_tot
        if (in_bounds(l)) then
           ijkp(cnt,:) = ijkp(l,:)
           cnt = cnt+1
        endif
      enddo


      call MPI_ISEND(out_ijkp, 3*Ni_out, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_ijkp, 3*Ni_in, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      ijkp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_ijkp(:,:)


      do m = 1,8
         out_part_wght(1:Ni_out,m) = 
     x          pack(wght(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo


      call pack_pd(wght, in_bounds, 8)

      call MPI_ISEND(out_part_wght, 8*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part_wght, 8*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      
      wght(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part_wght(:,:)



      out_beta_p(1:Ni_out) = 
     x     pack(beta_p(1:Ni_tot), .not.in_bounds(1:Ni_tot))

      call pack_pd(beta_p, in_bounds, 1)

      call MPI_ISEND(out_beta_p, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      
      beta_p(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      if( my_rank .eq. procnum-1 ) then
        out_tags(1:Ni_out) = 2
      else
        out_tags(1:Ni_out) =
     x     pack(tags(1:Ni_tot), .not.in_bounds(1:Ni_tot))
      endif

      call pack_pd(tags, in_bounds, 1)

      call MPI_ISEND(out_tags, Ni_out, MPI_REAL, dest, tag,
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_tags, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)

      tags(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_tags(:)

      out_mass(1:Ni_out) = 
     x     pack(mrat(1:Ni_tot), .not.in_bounds(1:Ni_tot))

      call pack_pd(mrat, in_bounds, 1)

      ! remove energy of outgoing particles
      input_E = input_E - sum(0.5*(mion/out_mass(:))*vsqrd_out(:)/
     x                        (beta*out_beta_p(:)))


      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      

      mrat(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      Ni_tot = Ni_tot - Ni_out + Ni_in

      ! add energy of incoming particles
      do m = 1,3
         input_E = input_E + sum(0.5*(mion/mrat(Ni_tot_in+1:Ni_tot))*
     x               (vp(Ni_tot_in+1:Ni_tot,m)*km_to_m)**2 /
     x                (beta*beta_p(Ni_tot_in+1:Ni_tot)))
      enddo
      
      deallocate(out_part)
      deallocate(out_ijkp)
      deallocate(out_part_wght)
      deallocate(out_mass)
      deallocate(out_beta_p)
      deallocate(out_tags)
      deallocate(vsqrd_out)

      deallocate(in_part)
      deallocate(in_ijkp)
      deallocate(in_part_wght)
      deallocate(in_mass)
      deallocate(in_tags)


c ---------z exchange, down---------------------------------------

      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      where (xp(1:Ni_tot,3) .le. qz(2))
         in_bounds(1:Ni_tot)= .false.
         ijkp(1:Ni_tot,3) = nz
         xp(1:Ni_tot,3) = qz(nz)-(qz(2)-xp(1:Ni_tot,3))
      endwhere

      Ni_tot_in = count(in_bounds(1:Ni_tot))

      Ni_out = count(.not.in_bounds(1:Ni_tot))


      allocate(out_part(Ni_out,3))
      allocate(out_ijkp(Ni_out,3))
      allocate(out_part_wght(Ni_out,8))
      allocate(out_mass(Ni_out))
      allocate(vsqrd_out(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(out_tags(Ni_out))



      dest = down_proc
      source = up_proc
      
      call MPI_ISEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)


      call MPI_WAITALL(2, reqs, stats, ierr)

      if(Ni_tot_in + Ni_in .gt. Ni_max) then
          write(error_unit,*) 'Not enough space to pass particles in2',
     x                                                      my_rank
          stop 1
      endif


      allocate(in_part(Ni_in,3))
      allocate(in_ijkp(Ni_in,3))
      allocate(in_part_wght(Ni_in,8))
      allocate(in_mass(Ni_in))
      allocate(in_tags(Ni_in))


      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(xp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(xp, in_bounds, 3)

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      xp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp, in_bounds, 3)

      vsqrd_out(:) = out_part(:,1)**2 + 
     x               out_part(:,2)**2 + out_part(:,3)**2
      vsqrd_out(:) = vsqrd_out*(km_to_m)**2

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      vp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp1(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp1, in_bounds, 3)

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      vp1(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_ijkp(1:Ni_out,m) = 
     x          pack(ijkp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo


      cnt = 1
      do l = 1,Ni_tot
        if (in_bounds(l)) then
           ijkp(cnt,:) = ijkp(l,:)
           cnt = cnt+1
        endif
      enddo


      call MPI_ISEND(out_ijkp, 3*Ni_out, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_ijkp, 3*Ni_in, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      ijkp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_ijkp(:,:)


      do m = 1,8
         out_part_wght(1:Ni_out,m) = 
     x          pack(wght(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
      enddo

      call pack_pd(wght, in_bounds, 8)

      call MPI_ISEND(out_part_wght, 8*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part_wght, 8*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      wght(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part_wght(:,:)



      out_beta_p(1:Ni_out) = 
     x     pack(beta_p(1:Ni_tot), .not.in_bounds(1:Ni_tot))


      call pack_pd(beta_p, in_bounds, 1)

      call MPI_ISEND(out_beta_p, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      beta_p(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      if( my_rank .eq. 0 ) then
        out_tags(1:Ni_out) = 3
      else
        out_tags(1:Ni_out) =
     x     pack(tags(1:Ni_tot), .not.in_bounds(1:Ni_tot))
      endif

      call pack_pd(tags, in_bounds, 1)

      call MPI_ISEND(out_tags, Ni_out, MPI_REAL, dest, tag,
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_tags, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)

      tags(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_tags(:)


      out_mass(1:Ni_out) = 
     x     pack(mrat(1:Ni_tot), .not.in_bounds(1:Ni_tot))

      call pack_pd(mrat, in_bounds, 1)

      ! remove energy of outgoing particles
      input_E = input_E - sum(0.5*(mion/out_mass(:))*vsqrd_out(:)/
     x                        (beta*out_beta_p(:)))

      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      mrat(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      Ni_tot = Ni_tot - Ni_out + Ni_in

      ! add energy of incoming particles
      do m = 1,3
         input_E = input_E + sum(0.5*(mion/mrat(Ni_tot_in+1:Ni_tot))*
     x               (vp(Ni_tot_in+1:Ni_tot,m)*km_to_m)**2 /
     x                (beta*beta_p(Ni_tot_in+1:Ni_tot)))
      enddo

      deallocate(out_part)
      deallocate(out_ijkp)
      deallocate(out_part_wght)
      deallocate(out_mass)
      deallocate(out_beta_p)
      deallocate(out_tags)
      deallocate(vsqrd_out)

      deallocate(in_part)
      deallocate(in_ijkp)
      deallocate(in_part_wght)
      deallocate(in_mass)
      deallocate(in_tags)

      return
      end SUBROUTINE move_ion_half
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_interp_weights_2(xp)
c Weights are used for trilinear interpolation to/from main cell
c centers to particle positions.  For each particle there are 8
c grid points associated with the interpolation.  These 8 points
c are determined by the location of the particle within the main
c cell.  There are 8 sets of 8 grid points for each cell.
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3)
      real x1,x2,y1,y2,z1,z2,vol

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c      where((xp(:,1) .le. qx(ijkp(:,1))) .and. 
c     x    (xp(:,2) .le. qy(ijkp(:,2))) .and.
c     x    (xp(:,3) .le. qz(ijkp(:,3)))) 
c         wquad(:,1) = -1             
c         wquad(:,2) = -1
c         wquad(:,3) = -1  
c         wght(:,1) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,2) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,3) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,4) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,5) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,6) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,7) = abs(xp(:,1)-qx(ijkp(:,1)))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))
c         wght(:,8) = abs(xp(:,1)-qx(ijkp(:,1)-1))*
c     x               abs(xp(:,2)-qy(ijkp(:,2)-1))*
c     x               abs(xp(:,3)-qz(ijkp(:,3)-1))/
c     x               (dx*dy*(qz(ijkp(:,3))-qz(ijkp(:,3)-1)))

c      endwhere


      do 10 l=1,Ni_tot

         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,2) .gt. ny) .or. 
     x       (ijkp(l,3) .gt. nz) .or. (ijkp(l,1) .lt. 1) .or. 
     x       (ijkp(l,2) .lt. 1) .or. (ijkp(l,3) .lt. 2)) then
            write(*,*) 'Out of bounds...',l,my_rank,
     x           ijkp(l,:),xp(l,:)
c            call remove_ion(xp,vp,vp1,l)
            endif
   

c 111111111111111111111111111111111111111111111111111111111111111111111


      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = -1
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 222222222222222222222222222222222222222222222222222222222222222222222

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = -1
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 333333333333333333333333333333333333333333333333333333333333333333333

      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = -1
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol


      endif

c 444444444444444444444444444444444444444444444444444444444444444444444

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .le. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = -1
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j-1))
         y2=abs(xp(l,2)-qy(j))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 555555555555555555555555555555555555555555555555555555555555555555555

      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = 0
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 666666666666666666666666666666666666666666666666666666666666666666666

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .le. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = 0
c         wquad(l,3) = -1
         vol = dx*dy*(qz(k)-qz(k-1))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k-1))
         z2=abs(xp(l,3)-qz(k))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 777777777777777777777777777777777777777777777777777777777777777777777

      if ((xp(l,1) .le. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = -1
c         wquad(l,2) = 0
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i-1))
         x2=abs(xp(l,1)-qx(i))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c 888888888888888888888888888888888888888888888888888888888888888888888

      if ((xp(l,1) .gt. qx(i)) .and. (xp(l,2) .gt. qy(j)) .and. 
     x    (xp(l,3) .gt. qz(k))) then
c         wquad(l,1) = 0
c         wquad(l,2) = 0
c         wquad(l,3) = 0
         vol = dx*dy*(qz(k+1)-qz(k))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2/vol
         wght(l,2) = x1*y2*z2/vol
         wght(l,3) = x2*y2*z1/vol
         wght(l,4) = x1*y2*z1/vol
         wght(l,5) = x2*y1*z2/vol
         wght(l,6) = x1*y1*z2/vol
         wght(l,7) = x2*y1*z1/vol
         wght(l,8) = x1*y1*z1/vol

      endif

c      if ((sum(wght(l,:)) .gt. 1.01) .or. (sum(wght(l,:)).lt.0.9)) then
c         write(*,*) 'wght error...',sum(wght(l,:))
c      endif

      wght(l,:) = wght(l,:)/beta_p(l)

c      if (sum(wght(l,:)) .gt. 1.01) then
c         write(*,*) 'wght dense...',sum(wght(l,:)),beta_p(l)
c      endif



 10   continue

      
      
      return
      end SUBROUTINE get_interp_weights_2
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_interp_weights(xp)
c Weights are used for trilinear interpolation to/from main cell
c centers to particle positions.  For each particle there are 8
c grid points associated with the interpolation.  These 8 points
c are determined by the location of the particle within the main
c cell.  There are 8 sets of 8 grid points for each cell.
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real xp(Ni_max,3)
      real x1,x2,y1,y2,z1,z2,vol

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do 10 l=1,Ni_tot

c         i = floor(xp(l,1)/dx) 
c         ijkp(l,1) = i

         i=0
 16      continue
         i = i + 1
         if (xp(l,1) .gt. qx(i)) go to 16 !find i on non-uniform 
         i = i-1
         ijkp(l,1)= i


         j=0
 18      continue
         j = j + 1
         if (xp(l,2) .gt. qy(j)) go to 18 !find j on non-uniform 
         j = j-1
         ijkp(l,2)= j

c         j = floor(xp(l,2)/dy) 
c         ijkp(l,2) = j

         k=2
 15      continue
         k = k + 1
         if (xp(l,3) .gt. qz(k)) go to 15  !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k

c randomize particles within cell

c         if (pad_ranf() .lt. 0.05) then
c         xp(l,1) = qx(i) + pad_ranf()*dx_grid(i)
c         xp(l,2) = qy(j) + pad_ranf()*dy_grid(j)
c         xp(l,3) = qz(k) + pad_ranf()*dz_grid(k)
c         endif

c         vol = 1.0/(dx*dy*(qz(k+1)-qz(k)))
         vol = 1.0/((qx(i+1)-qx(i))*(qy(j+1)-qy(j))*(qz(k+1)-qz(k)))
         x1=abs(xp(l,1)-qx(i))
         x2=abs(xp(l,1)-qx(i+1))
         y1=abs(xp(l,2)-qy(j))
         y2=abs(xp(l,2)-qy(j+1))
         z1=abs(xp(l,3)-qz(k))
         z2=abs(xp(l,3)-qz(k+1))
         wght(l,1) = x2*y2*z2*vol
         wght(l,2) = x1*y2*z2*vol
         wght(l,3) = x2*y2*z1*vol
         wght(l,4) = x1*y2*z1*vol
         wght(l,5) = x2*y1*z2*vol
         wght(l,6) = x1*y1*z2*vol
         wght(l,7) = x2*y1*z1*vol
         wght(l,8) = x1*y1*z1*vol

c         wght(l,:) = wght(l,:)/beta_p(l)  !scale for non equal particle weights
c         if (beta_p(l) .ne. 1.0) then
c            write(*,*) 'beta_p....',beta_p(l)
c         endif

 10   continue


      return
      end SUBROUTINE get_interp_weights
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE update_np(xp, vp, vp1, np)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real us(ny,nz)

      real volb
      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               np(i,j,k)=0.0
 10            continue

      do l=1,Ni_tot

         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

         if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1 .or.
     x       ip .gt. nx .or. jp .gt. ny .or. kp .gt. nz) then
            call remove_ion(xp,vp,vp1,l)
            cycle
         endif


         volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))

         np(i,j,k) = np(i,j,k) + wght(l,1)*volb
         np(ip,j,k) = np(ip,j,k) + wght(l,2)*volb
         np(i,j,kp) = np(i,j,kp) + wght(l,3)*volb
         np(ip,j,kp) = np(ip,j,kp) + wght(l,4)*volb
         np(i,jp,k) = np(i,jp,k) + wght(l,5)*volb
         np(ip,jp,k) = np(ip,jp,k) + wght(l,6)*volb
         np(i,jp,kp) = np(i,jp,kp) + wght(l,7)*volb
         np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)*volb

      enddo

c use for periodic boundary conditions
         np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)

      us = nf_init
      call boundary_scalar(np, us)

      return
      end SUBROUTINE update_np
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE separate_np(np,mask)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real us(ny,nz)
      logical mask(Ni_max)

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      real volb

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               np(i,j,k)=0.0
 10         continue
               

      do 20 l=1,Ni_tot
         
         if (mask(l)) then 
            

            i=ijkp(l,1)
            j=ijkp(l,2)
            k=ijkp(l,3)
            
            ip = i+1
            jp = j+1
            kp = k+1

         if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1 .or.
     x       ip .gt. nx .or. jp .gt. ny .or. kp .gt. nz) then
            !Remove ion with the separate flag
            call remove_ion(xp,vp,vp1,l,1) 
            cycle
         endif
            
            volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))
            
            np(i,j,k) = np(i,j,k) + wght(l,1)*volb
            np(ip,j,k) = np(ip,j,k) + wght(l,2)*volb
            np(i,j,kp) = np(i,j,kp) + wght(l,3)*volb
            np(ip,j,kp) = np(ip,j,kp) + wght(l,4)*volb
            np(i,jp,k) = np(i,jp,k) + wght(l,5)*volb
            np(ip,jp,k) = np(ip,jp,k) + wght(l,6)*volb
            np(i,jp,kp) = np(i,jp,kp) + wght(l,7)*volb
            np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)*volb

         endif
         
 20   continue
      
c     use for periodic boundary conditions
      np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)
      
      us = nf_init
      call boundary_scalar(np, us)

      call update_np_boundary(np)

      return
      end SUBROUTINE separate_np
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE update_up(vp,np,up)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vp(Ni_max,3),
     x     np(nx,ny,nz),
     x     up(nx,ny,nz,3)

      real volb,nvolb      !np times vol times beta

      real cnt(nx,ny,nz)

      integer dest, source
      real out_buf_z(nx,ny,3)
      real in_buf_z(nx,ny,3)
      integer cnt_buf_z
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      cnt_buf_z = nx*ny*3


      do 10 m=1,3          !clear out temp variable ct
         do 10 i=1,nx
            do 10 j=1,ny
               do 10 k=1,nz
                  up(i,j,k,m)=0.0
                  ct(i,j,k,m)=0.0
 10            continue
               
               cnt(:,:,:) = 0.0
    
      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)/beta_p(l) 
         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)/beta_p(l) 
         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)/beta_p(l) 
         
         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)
         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)/beta_p(l) 
         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)/beta_p(l) 
         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)/beta_p(l) 
         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)/beta_p(l) 
         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)/beta_p(l) 
         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l) 
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)/beta_p(l) 
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)/beta_p(l) 
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)/beta_p(l) 
         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)/beta_p(l) 
         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)/beta_p(l) 
         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)/beta_p(l) 
         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)/beta_p(l) 
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)/beta_p(l) 
         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)/beta_p(l) 
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)/beta_p(l) 
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)/beta_p(l) 
         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)/beta_p(l) 
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)/beta_p(l) 
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)/beta_p(l) 

 20   continue

c use for periodic boundary conditions
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      where (cnt(:,:,:) .gt. 0.0)

         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere

      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      up = ct

      call boundaries(up, vsw_us)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end SUBROUTINE update_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE update_np_boundary(np)
c----------------------------------------------------------------------

      real np(nx,ny,nz)
      real us(ny,nz)

      integer dest, source
      real out_buf_z(nx,ny)
      real in_buf_z(nx,ny)
      integer cnt_buf_z
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)

      cnt_buf_z = nx*ny

      out_buf_z(:,:) = np(:,:,nz)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      np(:,:,2) = np(:,:,2) + in_buf_z
 
      us = nf_init
      call boundary_scalar(np, us)
     
      return 
      end SUBROUTINE update_np_boundary
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE pack_pd(arr,larr,sz)
c replaces the f90 pack function
c----------------------------------------------------------------------
      
      integer*4 sz
      real arr(Ni_max,sz)
      logical larr(Ni_max)
      integer*4 cnt

      cnt = 1
      do l = 1,Ni_tot
        if (larr(l)) then
           arr(cnt,:) = arr(l,:)
           cnt = cnt+1
        endif
      enddo

      return
      end SUBROUTINE pack_pd
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE get_temperature(xp,vp,np,temp_p)
c----------------------------------------------------------------------
      
      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     np(nx,ny,nz),
     x     up_ave(nx,ny,nz,3),
     x     up2(nx,ny,nz,3),
     x     temp_p(nx,ny,nz)

      integer dest, source
      real out_buf_z(nx,ny,3)
      real in_buf_z(nx,ny,3)
      integer cnt_buf_z
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)

      real cnt(nx,ny,nz)

      real mvp(Ni_max,3)

      real volb,nvolb      !np times vol times beta

      cnt_buf_z = nx*ny*3

      cnt(:,:,:) = 0.0

      up2(:,:,:,:) = 0.0
      up_ave(:,:,:,:) = 0.0
      ct(:,:,:,:) = 0.0

      do m = 1,3 
         mvp(1:Ni_tot,m) = vp(1:Ni_tot,m)/sqrt(mrat(1:Ni_tot))
      enddo


      do 20 l=1,Ni_tot


         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         ip = i+1
         jp = j+1
         kp = k+1
      
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/beta_p(l)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/beta_p(l)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/beta_p(l)

         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)         
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/beta_p(l)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/beta_p(l)

         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/beta_p(l)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/beta_p(l)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/beta_p(l)

         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/beta_p(l)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/beta_p(l)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/beta_p(l)

         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/beta_p(l)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/beta_p(l)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/beta_p(l)

         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/beta_p(l)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/beta_p(l)

         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/beta_p(l)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/beta_p(l)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/beta_p(l)

         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1)+mvp(l,1)**2*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2)+mvp(l,2)**2*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3)+mvp(l,3)**2*wght(l,8)/beta_p(l)


 20   continue

c use for periodic boundary conditions
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)


      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)
      
      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1

               up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue


      call periodic(up2)


      ct(:,:,:,:) = 0.0
      cnt(:,:,:) = 0.0


      do 40 l=1,Ni_tot

         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         ip = i+1
         jp = j+1
         kp = k+1
      
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/beta_p(l)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/beta_p(l)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/beta_p(l)

         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/beta_p(l)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/beta_p(l)

         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/beta_p(l)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/beta_p(l)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/beta_p(l)

         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/beta_p(l)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/beta_p(l)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/beta_p(l)

         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/beta_p(l)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/beta_p(l)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/beta_p(l)

         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/beta_p(l)

         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/beta_p(l)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/beta_p(l)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/beta_p(l)

         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/beta_p(l)


 40   continue

c use for periodic boundary conditions
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)

      do 50 i=1,nx-1      !interpolate back to contravarient positions
         do 50 j=1,ny-1
            do 50 k=1,nz-1

               up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 50            continue


      call periodic(up_ave)


      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               temp_p(i,j,k) = (1./3.)*1e6*mproton*(up2(i,j,k,1)+
     x                              up2(i,j,k,2) + 
     x                              up2(i,j,k,3) - up_ave(i,j,k,1)**2 -
     x                              up_ave(i,j,k,2)**2 -
     x                              up_ave(i,j,k,3)**2)
            enddo
         enddo
      enddo
      

      call periodic_scalar(temp_p)


      return
      end SUBROUTINE get_temperature
c----------------------------------------------------------------------

      SUBROUTINE separate_temperature(xp,vp,np,temp_p,mask)
c----------------------------------------------------------------------
      
      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     np(nx,ny,nz),
     x     up_ave(nx,ny,nz,3),
     x     up2(nx,ny,nz,3),
     x     temp_p(nx,ny,nz)
      logical mask(Ni_max)

      integer dest, source
      real out_buf_z(nx,ny,3)
      real in_buf_z(nx,ny,3)
      integer cnt_buf_z
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)

      real cnt(nx,ny,nz)

      real mvp(Ni_max,3)

      real volb,nvolb      !np times vol times beta

      cnt_buf_z = nx*ny*3

      cnt(:,:,:) = 0.0

      up2(:,:,:,:) = 0.0
      up_ave(:,:,:,:) = 0.0
      ct(:,:,:,:) = 0.0


      do m = 1,3 
         mvp(1:Ni_tot,m) = vp(1:Ni_tot,m)/sqrt(mrat(1:Ni_tot))
      enddo


      do 20 l=1,Ni_tot
      if(mask(l)) then


         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         ip = i+1
         jp = j+1
         kp = k+1
      
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/beta_p(l)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/beta_p(l)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/beta_p(l)

         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)         
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/beta_p(l)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/beta_p(l)

         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/beta_p(l)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/beta_p(l)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/beta_p(l)

         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/beta_p(l)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/beta_p(l)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/beta_p(l)

         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/beta_p(l)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/beta_p(l)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/beta_p(l)

         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/beta_p(l)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/beta_p(l)

         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/beta_p(l)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/beta_p(l)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/beta_p(l)

         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1)+mvp(l,1)**2*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2)+mvp(l,2)**2*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3)+mvp(l,3)**2*wght(l,8)/beta_p(l)


      endif
 20   continue

c use for periodic boundary conditions
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)
      
      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1

               up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue


      call periodic(up2)


      ct(:,:,:,:) = 0.0
      cnt(:,:,:) = 0.0


      do 40 l=1,Ni_tot
      if(mask(l)) then

         i=ijkp(l,1)
         j=ijkp(l,2)
         k=ijkp(l,3)

         ip = i+1
         jp = j+1
         kp = k+1
      
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/beta_p(l)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/beta_p(l)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/beta_p(l)

         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/beta_p(l)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/beta_p(l)

         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/beta_p(l)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/beta_p(l)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/beta_p(l)

         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/beta_p(l)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/beta_p(l)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/beta_p(l)

         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/beta_p(l)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/beta_p(l)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/beta_p(l)

         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/beta_p(l)

         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/beta_p(l)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/beta_p(l)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/beta_p(l)

         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/beta_p(l)


      endif
 40   continue

c use for periodic boundary conditions
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)

      do 50 i=1,nx-1      !interpolate back to contravarient positions
         do 50 j=1,ny-1
            do 50 k=1,nz-1

               up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 50            continue


      call periodic(up_ave)


      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               temp_p(i,j,k) = (1./3.)*1e6*mproton*(up2(i,j,k,1)+
     x                              up2(i,j,k,2) + 
     x                              up2(i,j,k,3) - up_ave(i,j,k,1)**2 -
     x                              up_ave(i,j,k,2)**2 -
     x                              up_ave(i,j,k,3)**2)
            enddo
         enddo
      enddo
      

      call periodic_scalar(temp_p)


      return
      end SUBROUTINE separate_temperature


      end MODULE gutsp_dd


