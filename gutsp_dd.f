      MODULE gutsp_dd

      USE global
      USE misc
      USE mpi
      USE boundary
      USE grid_interp

      contains

c----------------------------------------------------------------------
      SUBROUTINE remove_ion(xp,vp,vp1,ion_l)
c Removes particles from simulation that have gone out of bounds
c----------------------------------------------------------------------
CVD$R VECTOR

c      include 'incurv.h'

      real xp(Ni_max,3)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)

      do 5 m=1,3   !remove ion energy from total input energy
         input_E = input_E-0.5*(mion/mrat(ion_l))*(vp(ion_l,m)*km_to_m)**2 /
     x             (beta*beta_p(ion_l))
 5    continue
      write(*,*) 'removing ion...',ion_l

      do 10 l=ion_l,Ni_tot-1
c         m_arr(l) = m_arr(l+1)
         mrat(l) = mrat(l+1)
         do 10 m=1,3 
            xp(l,m) = xp(l+1,m)
            vp(l,m) = vp(l+1,m)
            vp1(l,m) = vp1(l+1,m)
            ijkp(l,m) = ijkp(l+1,m)
c            wquad(l,m) = wquad(l+1,m)
 10      continue



      do 20 m=1,8
         do 20 l=ion_l,Ni_tot-1
            wght(l,m) = wght(l+1,m)
 20         continue

      Ni_tot = Ni_tot - 1

      return
      end SUBROUTINE remove_ion
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE check_min_den(np,xp,vp,vp1,up,bt)
c----------------------------------------------------------------------
c      include 'incurv.h'

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
      real, dimension(:), allocatable :: vsqrd_out

      real minden !minimum wake density 
      real den_part ! density of 1 particle per cell
      real ak
      real btot,a1,a2,womega,phi,deltat

      den_part = 1/(beta*dx**3)

      minden = nf_init/5.
c      minden = 2.0*den_part
      do i = 2,nx-1
         do j = 2,ny-1
            do k = 2,nz-1
c               ak = PI/dx
c               btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + 
c     x              bt(i,j,k,3)**2)
c               a1 = ak**2*Btot/(alpha*(np(i,j,k)))
c               a2 = (ak*Btot)**2/(alpha*(np(i,j,k)))
c               womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
c               phi = womega/ak
c               deltat = 0.1*dx/phi
c               if (deltat .le. dtsub) then
c                  write(*,*) 'Time stepping error...'
c               endif
               if ((np(i,j,k) .le. minden) .or. 
     x              (np(i,j,k) .le. 2.0*den_part)) then
                  npart = nint(minden/(np(i,j,k)))
                  do ipart = 1,npart 
c                     write(*,*) 'min den...',np(i,j,k),min_den,den_part,
c     x                                  npart,ipart
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
                     
c     j = floor(xp(l,2)/dy) 
c     ijkp(l,2) = j
                     
                     kk=0
 15                  continue
                     kk = kk + 1
                     if (xp(l,3) .gt. qz(kk)) go to 15 !find k on non-uniform 
                     kk = kk-1
                     ijkp(l,3)= kk



c                     ijkp(l,1) = floor(xp(l,1)/dx) !particle grid location index
c                     ijkp(l,2) = floor(xp(l,2)/dy)
                     
c                     kk=1
c                     do 100 while(xp(l,3) .gt. qz(kk)) !find k on non-uniform 
c                        ijkp(l,3) = kk !grid
c                        kk=kk+1
c 100                 continue
c                     kk=ijkp(l,3)
c                     if (xp(l,3) .gt. (qz(kk)+(dz_grid(kk)/2))) then
c                        ijkp(l,3) = kk+1
c                     endif


                     mrat(l) = 1.0
c                     m_arr(l) = mproton
                     beta_p(l) = 1.0
                     Ni_tot = Ni_tot + 1
                  enddo
               endif
            enddo
         enddo
      enddo
c ---------z exchange down---------------------------------------

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      where (xp(:,3) .le. qz(1)) 
c         ijkp(:,3) = nz
c         wquad(:,3) = -1.0
c         xp(:,3) = qz(nz) - (qz(1) - xp(:,3))
c      endwhere

      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      where (xp(1:Ni_tot,3) .le. qz(2))
         in_bounds(1:Ni_tot)= .false.
         ijkp(1:Ni_tot,3) = nz
c         wquad(1:Ni_tot,3) = -1.0
         xp(1:Ni_tot,3) = qz(nz)-(qz(2)-xp(1:Ni_tot,3))
      endwhere

      Ni_tot_in = count(in_bounds(1:Ni_tot))

      Ni_out = count(.not.in_bounds(1:Ni_tot))
      allocate(out_part(Ni_out,3))
      allocate(out_ijkp(Ni_out,3))
      allocate(out_part_wght(Ni_out,8))
      allocate(out_mass(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(vsqrd_out(Ni_out))

      dest = nbrs(n_down)
      source = nbrs(n_up)
      
      call MPI_ISEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)


      call MPI_WAITALL(2, reqs, stats, ierr)
c      write(*,*) 'down exchange...',Ni_in,Ni_out,my_rank      

      allocate(in_part(Ni_in,3))
      allocate(in_ijkp(Ni_in,3))
      allocate(in_part_wght(Ni_in,8))
      allocate(in_mass(Ni_in))


      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(xp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c         vp(1:Ni_tot_in,m) = pack(vp(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c         vp1(1:Ni_tot_in,m) = pack(vp1(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c         ijkp(1:Ni_tot_in,m) = 
c     x          pack(ijkp(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo

c     call pack_pd_3(ijkp, in_bounds)

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


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x          pack(wquad(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         wquad(1:Ni_tot_in,m) = 
c     x          pack(wquad(1:Ni_tot,m), in_bounds(1:Ni_tot))
c      enddo

c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)
      
c      wquad(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,8
         out_part_wght(1:Ni_out,m) = 
     x          pack(wght(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         wght(1:Ni_tot_in,m) = 
c    x          pack(wght(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c      beta_p(1:Ni_tot_in) = pack(beta_p(1:Ni_tot), in_bounds(1:Ni_tot))

      call pack_pd(beta_p, in_bounds, 1)
      
      call MPI_ISEND(out_beta_p, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      beta_p(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


c      out_mass(1:Ni_out) = 
c     x     pack(m_arr(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      m_arr(1:Ni_tot_in) = pack(m_arr(1:Ni_tot), in_bounds(1:Ni_tot))

c      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      m_arr(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      out_mass(1:Ni_out) = 
     x     pack(mrat(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      mrat(1:Ni_tot_in) = pack(mrat(1:Ni_tot), in_bounds(1:Ni_tot))

      call pack_pd(mrat, in_bounds, 1)

      ! remove energy of outgoing particles
      input_E = input_E - sum(0.5*(mion/out_mass(:))*vsqrd_out(:)/
     x                       (beta*out_beta_p(:)))

      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
      mrat(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      Ni_tot = Ni_tot - Ni_out + Ni_in

      ! add energy back in
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
      deallocate(vsqrd_out)

      deallocate(in_part)
      deallocate(in_ijkp)
      deallocate(in_part_wght)
      deallocate(in_mass)

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

      call face_to_center(aj,ajc)
      call face_to_center(up,upc)
c      call face_to_center(uf,ufc)
c      call face_to_center(gradP,gradPc)

      do 10 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

            
c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         np_at_Ba = np(i,j,k)*wght(l,1) + np(ip,j,k)*wght(l,2) + 
c     x              np(i,j,kp)*wght(l,3) + np(ip,j,kp)*wght(l,4) + 
c     x              np(i,jp,k)*wght(l,5) + np(ip,jp,k)*wght(l,6) +
c     x              np(i,jp,kp)*wght(l,7) + np(ip,jp,kp)*wght(l,8)

c         nf_at_Ba = nf(i,j,k)*wght(l,1) + nf(ip,j,k)*wght(l,2) + 
c     x              nf(i,j,kp)*wght(l,3) + nf(ip,j,kp)*wght(l,4) + 
c     x              nf(i,jp,k)*wght(l,5) + nf(ip,jp,k)*wght(l,6) +
c     x              nf(i,jp,kp)*wght(l,7) + nf(ip,jp,kp)*wght(l,8)

c         ntot = nf_at_Ba + np_at_Ba
c         fnf = nf_at_Ba/ntot
c         fnp = np_at_Ba/ntot

         do 15 m=1,3 
            aj3(m) = ajc(i,j,k,m)*wght(l,1) + ajc(ip,j,k,m)*wght(l,2) 
     x          + ajc(i,j,kp,m)*wght(l,3) + ajc(ip,j,kp,m)*wght(l,4)
     x          + ajc(i,jp,k,m)*wght(l,5) + ajc(ip,jp,k,m)*wght(l,6)
     x          + ajc(i,jp,kp,m)*wght(l,7) + ajc(ip,jp,kp,m)*wght(l,8)

            up3(m) = upc(i,j,k,m)*wght(l,1) + upc(ip,j,k,m)*wght(l,2) 
     x          + upc(i,j,kp,m)*wght(l,3) + upc(ip,j,kp,m)*wght(l,4)
     x          + upc(i,jp,k,m)*wght(l,5) + upc(ip,jp,k,m)*wght(l,6)
     x          + upc(i,jp,kp,m)*wght(l,7) + upc(ip,jp,kp,m)*wght(l,8)

c            uf3(m) = ufc(i,j,k,m)*wght(l,1) + ufc(ip,j,k,m)*wght(l,2) 
c     x          + ufc(i,j,kp,m)*wght(l,3) + ufc(ip,j,kp,m)*wght(l,4)
c     x          + ufc(i,jp,k,m)*wght(l,5) + ufc(ip,jp,k,m)*wght(l,6)
c     x          + ufc(i,jp,kp,m)*wght(l,7) + ufc(ip,jp,kp,m)*wght(l,8)

            btc3(m) = btc(i,j,k,m)*wght(l,1) 
     x               + btc(ip,j,k,m)*wght(l,2) 
     x               + btc(i,j,kp,m)*wght(l,3) 
     x               + btc(ip,j,kp,m)*wght(l,4)
     x               + btc(i,jp,k,m)*wght(l,5) 
     x               + btc(ip,jp,k,m)*wght(l,6)
     x               + btc(i,jp,kp,m)*wght(l,7) 
     x               + btc(ip,jp,kp,m)*wght(l,8)

c            gradP3(m) = gradPc(i,j,k,m)*wght(l,1) 
c     x                + gradPc(ip,j,k,m)*wght(l,2) 
c     x                + gradPc(i,j,kp,m)*wght(l,3) 
c     x                + gradPc(ip,j,kp,m)*wght(l,4)
c     x                + gradPc(i,jp,k,m)*wght(l,5) 
c     x                + gradPc(ip,jp,k,m)*wght(l,6)
c     x                + gradPc(i,jp,kp,m)*wght(l,7) 
c     x                + gradPc(ip,jp,kp,m)*wght(l,8)

 15         continue

         do 20 m=1,3
c            aa(m) = aj3(m) - fnp*up3(m)  
c     x                     - fnf*uf3(m)
            aa(m) = aj3(m) - up3(m)
            bb(m) = btc3(m)                   

 20         continue

         cc(1) = aa(2)*bb(3) - aa(3)*bb(2)    !do cross product
         cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
         cc(3) = aa(1)*bb(2) - aa(2)*bb(1)

         do 30 m=1,3
            Ep(l,m) = cc(m) 
c     x                + nu(i,j,k)*fnf*(uf3(m)-up3(m)) 
c     x                - gradP3(m) 
c     x                + nuei*aj3(m) 
c                     + etar(i,j,k,m)*aj3(m)
            Ep(l,m) = Ep(l,m)*mrat(l) !O_to_Ba

 30         continue


 10      continue

      return
      end SUBROUTINE get_Ep
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vplus_vminus(Ep,btc,vp,vplus,vminus)
c----------------------------------------------------------------------

c      include 'incurv.h'

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

c      do 20 l=1,Ni_tot

c         i = ijkp(l,1)+wquad(l,1)
c         j = ijkp(l,2)+wquad(l,2)
c         k = ijkp(l,3)+wquad(l,3)
         
c         ip=i+1
c         jp=j+1
c         kp=k+1

c         if (ip .ge. nx) ip = 2      !periodic boundary conditions
c         if (jp .ge. ny) jp = 2      !periodic boundary conditions
c         if (kp .ge. nz) kp = 2      !periodic boundary conditions

c         do 25 m=1,3    !interpolate to particle position
c            btc3(l,m) = btc(i,j,k,m)*wght(l,1) 
c     x               + btc(ip,j,k,m)*wght(l,2) 
c     x               + btc(i,j,kp,m)*wght(l,3) 
c     x               + btc(ip,j,kp,m)*wght(l,4)
c     x               + btc(i,jp,k,m)*wght(l,5) 
c     x               + btc(ip,jp,k,m)*wght(l,6)
c     x               + btc(i,jp,kp,m)*wght(l,7) 
c     x               + btc(ip,jp,kp,m)*wght(l,8)
c 25         continue

c         vminus_x_B(l,1) = vminus(l,2)*btc3(l,3)*mrat(l) - !O_to_Ba - 
c     x                     vminus(l,3)*btc3(l,2)*mrat(l)   !O_to_Ba
c         vminus_x_B(l,2) = vminus(l,3)*btc3(l,1)*mrat(l) - !O_to_Ba - 
c     x                     vminus(l,1)*btc3(l,3)*mrat(l)   !O_to_Ba
c         vminus_x_B(l,3) = vminus(l,1)*btc3(l,2)*mrat(l) - !O_to_Ba -
c     x                     vminus(l,2)*btc3(l,1)*mrat(l)   !O_to_Ba

c         vminus_dot_B(l) = vminus(l,1)*btc3(l,1)*mrat(l) + !O_to_Ba +
c     x                     vminus(l,2)*btc3(l,2)*mrat(l) + !O_to_Ba +
c     x                     vminus(l,3)*btc3(l,3)*mrat(l)   !O_to_Ba

c 20   continue
   
      do 30 l=1,Ni_tot 

         i = ijkp(l,1)!+wquad(l,1)
         j = ijkp(l,2)!+wquad(l,2)
         k = ijkp(l,3)!+wquad(l,3)
   
         ip = i+1
         jp = j+1
         kp = k+1


c         if (ip .ge. nx) ip = 2 !periodic boundary conditions
c         if (jp .ge. ny) jp = 2 !periodic boundary conditions
c         if (kp .ge. nz) kp = 2 !periodic boundary conditions

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

         vminus_x_B(1) = vminus(l,2)*btc3(3)*mrat(l) - !O_to_Ba - 
     x                     vminus(l,3)*btc3(2)*mrat(l)   !O_to_Ba
         vminus_x_B(2) = vminus(l,3)*btc3(1)*mrat(l) - !O_to_Ba - 
     x                     vminus(l,1)*btc3(3)*mrat(l)   !O_to_Ba
         vminus_x_B(3) = vminus(l,1)*btc3(2)*mrat(l) - !O_to_Ba -
     x                     vminus(l,2)*btc3(1)*mrat(l)   !O_to_Ba

         vminus_dot_B = vminus(l,1)*btc3(1)*mrat(l) + !O_to_Ba +
     x                     vminus(l,2)*btc3(2)*mrat(l) + !O_to_Ba +
     x                     vminus(l,3)*btc3(3)*mrat(l)   !O_to_Ba

         Bx = btc3(1)*mrat(l) !O_to_Ba
         By = btc3(2)*mrat(l) !O_to_Ba
         Bz = btc3(3)*mrat(l) !O_to_Ba
      
         B2 = Bx**2 + By**2 + Bz**2
         dt2 = dt**2

         a_d = 1 + (B2*dt2/4.0)
         a1 = (1 - (B2*dt2/4.0)) / a_d
         a2 = dt / a_d
         a3 = 0.5*dt2 / a_d

         do 40 m=1,3
            vplus(l,m) = a1*vminus(l,m) + a2*vminus_x_B(m) + 
     x           a3*vminus_dot_B*btc3(m)*mrat(l) !O_to_Ba
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
      SUBROUTINE exchange_ion_in(xp,vp,vp1,input_p,xp_buf,vp_buf) 
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     input_p(3)
      real xp_buf(Ni_max_buf,3)
      real vp_buf(Ni_max_buf,3)
      real E(nx,ny,nz,3)
      real Bt(nx,ny,nz,3)
      real xp_out_buf(Ni_max_buf,3)
      real vp_out_buf(Ni_max_buf,3)
      real E_out_buf(Ni_max_buf,3)
      real B_out_buf(Ni_max_buf,3)
      real mrat_out_buf(Ni_max_buf)
c      real m_arr_out_buf(Ni_max_buf)

      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
c      real, dimension(:), allocatable :: out_m_arr
      real, dimension(:), allocatable :: out_mrat
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:,:), allocatable :: out_E
      real, dimension(:,:), allocatable :: out_B
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
c      allocate(out_m_arr(Ni_out))
      allocate(out_mrat(Ni_out))
      allocate(out_beta_p(Ni_out))
      
      do m = 1,3 
         
         out_xp(1:Ni_out,m) = pack(xp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         out_vp(1:Ni_out,m) = pack(vp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
c         xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), 
c     x        in_bounds(1:Ni_tot))
c         vp(1:Ni_tot_in,m) = pack(vp(1:Ni_tot,m), 
c     x        in_bounds(1:Ni_tot))
c         vp1(1:Ni_tot_in,m) = pack(vp1(1:Ni_tot,m), 
c     x        in_bounds(1:Ni_tot))
c         ijkp(1:Ni_tot_in,m) = pack(ijkp(1:Ni_tot,m),
c     x        in_bounds(1:Ni_tot))
cc         wquad(1:Ni_tot_in,m) = pack(wquad(1:Ni_tot,m),
cc     x        in_bounds(1:Ni_tot))
      enddo
 
      call pack_pd(xp, in_bounds, 3)
      call pack_pd(vp, in_bounds, 3)
      call pack_pd(vp1, in_bounds, 3)      
c      call pack_pd_3(ijkp, in_bounds)

      cnt = 1
      do l = 1,Ni_tot
        if (in_bounds(l)) then
           ijkp(cnt,:) = ijkp(l,:)
           cnt = cnt+1
        endif
      enddo

     
c      do m = 1,8
c         wght(1:Ni_tot_in,m) = pack(wght(1:Ni_tot,m),
c     x        in_bounds(1:Ni_tot))
c      enddo
 
      call pack_pd(wght, in_bounds, 8)
      
c      out_m_arr(1:Ni_out) = pack(m_arr(1:Ni_tot), 
c     x     .not.in_bounds(1:Ni_tot))
      out_mrat(1:Ni_out) = pack(mrat(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      out_beta_p(1:Ni_out) = pack(beta_p(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      
cc      m_arr(1:Ni_tot_in) = pack(m_arr(1:Ni_tot), 
cc     x     in_bounds(1:Ni_tot))
c      mrat(1:Ni_tot_in) = pack(mrat(1:Ni_tot), 
c     x     in_bounds(1:Ni_tot))
c      beta_p(1:Ni_tot_in) = pack(beta_p(1:Ni_tot), 
c     x     in_bounds(1:Ni_tot))
      
      call pack_pd(mrat, in_bounds, 1)
      call pack_pd(beta_p, in_bounds, 1)
            
      do m = 1,3
         xp_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out,m) = out_xp(:,m)
         vp_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out,m) = out_vp(:,m)
      enddo
      
c      m_arr_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out) = out_m_arr(:)
      mrat_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out) = out_mrat(:)
      beta_p_buf(Ni_tot_buf+1:Ni_tot_buf+Ni_out) = out_beta_p(:)
      
c     remove energy

      do l = 1,Ni_out  
         do m=1,3
            input_E = input_E - 
     x           0.5*(mion/out_mrat(l))*(out_vp(l,m)*km_to_m)**2 /
     x           (beta*out_beta_p(l))
c     input_p(m) = input_p(m) - (mion/mrat(l))*vp(l,m) / 
c     x                    (beta*beta_p(l))
         enddo         
      enddo

c      do l = Ni_tot_in+1,Ni_tot  
c         do m=1,3
c            input_E = input_E - 
c     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
c     x               (beta*beta_p(l))
c            input_p(m) = input_p(m) - (mion/mrat(l))*vp(l,m) / 
c     x                    (beta*beta_p(l))
c         enddo         
c      enddo
      
      Ni_tot_buf = Ni_tot_buf + Ni_out
      Ni_tot = count(in_bounds)
      
      
      deallocate(out_xp)
      deallocate(out_vp)
      deallocate(out_mrat)
      deallocate(out_beta_p)
      

      return
      end SUBROUTINE exchange_ion_in
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE exchange_ion_out(xp,vp,vp1,input_p,xp_buf,vp_buf,
     x     E,Bt,xp_out_buf,vp_out_buf,E_out_buf,
     x     B_out_buf,mrat_out_buf) 
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     input_p(3)
      real xp_buf(Ni_max_buf,3)
      real vp_buf(Ni_max_buf,3)
      real E(nx,ny,nz,3)
      real Bt(nx,ny,nz,3)
      real xp_out_buf(Ni_max_buf,3)
      real vp_out_buf(Ni_max_buf,3)
      real E_out_buf(Ni_max_buf,3)
      real B_out_buf(Ni_max_buf,3)
      real mrat_out_buf(Ni_max_buf)
c      real m_arr_out_buf(Ni_max_buf)

      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
c      real, dimension(:), allocatable :: out_m_arr
      real, dimension(:), allocatable :: out_mrat
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:,:), allocatable :: out_E
      real, dimension(:,:), allocatable :: out_B
      integer, dimension(:,:), allocatable :: out_ijkp
      integer Ni_tot_in, Ni_out
      integer :: cnt


      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.
      
c     where(xp(Ni_tot_sw+1:Ni_tot,1) .le. qx(1)) 
c     x            in_bounds(Ni_tot_sw+1:Ni_tot) = .false.
      
      where(xp(1:Ni_tot,1) .le. qx(1))
     X     in_bounds(1:Ni_tot) = .false.
c      endwhere         
      Ni_tot_in = count(in_bounds)
      Ni_out = count(.not.in_bounds(1:Ni_tot))
      
      allocate(out_xp(Ni_out,3))
      allocate(out_vp(Ni_out,3))
      allocate(out_E(Ni_out,3))
      allocate(out_B(Ni_out,3))
      allocate(out_ijkp(Ni_out,3))
c      allocate(out_m_arr(Ni_out))
      allocate(out_mrat(Ni_out))
      allocate(out_beta_p(Ni_out))
      

      do m = 1,3 
         out_xp(1:Ni_out,m) = pack(xp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         out_vp(1:Ni_out,m) = pack(vp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         out_ijkp(1:Ni_out,m) = pack(ijkp(1:Ni_tot,m), 
     x        .not.in_bounds(1:Ni_tot))
         
c         xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), in_bounds(1:Ni_tot))
c         vp(1:Ni_tot_in,m) = pack(vp(1:Ni_tot,m), in_bounds(1:Ni_tot))
c         vp1(1:Ni_tot_in,m) = pack(vp1(1:Ni_tot,m), in_bounds(1:Ni_tot))
c         ijkp(1:Ni_tot_in,m)=pack(ijkp(1:Ni_tot,m), in_bounds(1:Ni_tot))
cc        wquad(1:Ni_tot_in,m)=pack(wquad(1:Ni_tot,m),in_bounds(1:Ni_tot))
      enddo

      call pack_pd(xp, in_bounds, 3)
      call pack_pd(vp, in_bounds, 3)
      call pack_pd(vp1, in_bounds, 3)
c      call pack_pd_3(ijkp, in_bounds)

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


c      do m = 1,8
c         wght(1:Ni_tot_in,m)=pack(wght(1:Ni_tot,m), in_bounds(1:Ni_tot))
c      enddo
 
      call pack_pd(wght, in_bounds, 8)
     
c      out_m_arr(1:Ni_out) = pack(m_arr(1:Ni_tot), 
c     x     .not.in_bounds(1:Ni_tot))
      out_mrat(1:Ni_out) = pack(mrat(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      out_beta_p(1:Ni_out) = pack(beta_p(1:Ni_tot), 
     x     .not.in_bounds(1:Ni_tot))
      
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
      

c      m_arr_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out) = 
c     x     out_m_arr(:)
      
      mrat_out_buf(Ni_tot_out_buf+1:Ni_tot_out_buf+Ni_out) = 
     x     out_mrat(:)
      
cc      m_arr(1:Ni_tot_in) = pack(m_arr(1:Ni_tot), in_bounds(1:Ni_tot))
c      mrat(1:Ni_tot_in) = pack(mrat(1:Ni_tot), in_bounds(1:Ni_tot))
c      beta_p(1:Ni_tot_in) = pack(beta_p(1:Ni_tot), in_bounds(1:Ni_tot))

      call pack_pd(mrat, in_bounds, 1)
      call pack_pd(beta_p, in_bounds, 1)
            
c     remove energy

      do l = 1,Ni_out  
         do m=1,3
            input_E = input_E - 
     x           0.5*(mion/out_mrat(l))*(out_vp(l,m)*km_to_m)**2 /
     x           (beta*out_beta_p(l))
c     input_p(m) = input_p(m) - (mion/mrat(l))*vp(l,m) / 
c     x                    (beta*beta_p(l))
         enddo         
      enddo
      
c      do l = Ni_tot_in+1,Ni_tot  
c         do m=1,3
c            input_E = input_E - 
c     x           0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /
c     x           (beta*beta_p(l))
c            input_p(m) = input_p(m) - (mion/mrat(l))*vp(l,m) / 
c     x                   (beta*beta_p(l))
c         enddo         
c      enddo
      
c      Ni_tot_out_buf = Ni_tot_out_buf + Ni_out

c      if (Ni_tot_out_buf .ge. Ni_max_buf) then 
c         write(*,*) 'Error....out buf too small...',Ni_tot_out_buf,
c     x        Ni_max_buf
c      endif
      

      deallocate(out_xp)
      deallocate(out_vp)
      deallocate(out_E)
      deallocate(out_B)
      deallocate(out_ijkp)
c      deallocate(out_m_arr)
      deallocate(out_mrat)
      deallocate(out_beta_p)
      
      Ni_tot = count(in_bounds)

         
      return
      end SUBROUTINE exchange_ion_out
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE move_ion_half(xp,vp,vp1,input_p,Ep)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3),
     x     input_p(3)
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
      real, dimension(:), allocatable :: vsqrd_out
      integer :: cnt


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      dth = dt/2.

      do 10 l=1,Ni_tot                   !make 1/2 time step advance

         xp0 = xp(l,1)
         vp0 = vp(l,1)
         xp(l,1) = xp(l,1) + dth*vp(l,1)
c         ijkp(l,1) = nint(xp(l,1)/dx) 


c         if ((xp(l,1) .gt. qx(nx)) .or. (xp(l,1) .lt. qx(1))) then
c            write(*,*) 'WARNING...part OB...',xp(l,:),nint(xp(l,:)/dx)
c         endif
c         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,1) .lt. 1)) then
c            write(*,*) 'WARNING...part OB...',ijkp(l,:),xp(l,1),vp(l,1),
c     x                  xp0,vp0,vp1(l,1),Ep(l,:)
c            stop
c         endif

 


         xp(l,2) = xp(l,2) + dth*vp(l,2)
c         ijkp(l,2) = nint(xp(l,2)/dy) 

         xp(l,3) = xp(l,3) + dth*vp(l,3)
c         ijkp(l,3) = nint(xp(l,3)/delz)

c         k=1
c         do 15 while((xp(l,3) .gt. qz(k)) .and. (k .le. nz))  !find k
c            ijkp(l,3) = k                 !grid
c            k=k+1
c 15      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif

c         if (abs(xp(l,3)-qz(ijkp(l,3))) .gt. delz) then
c            write(*,*) 'k index error...',xp(l,3),ijkp(l,3)
c            stop
c         endif

c         if ((ijkp(l,3) .gt. nz) .or. (ijkp(l,3) .lt. 2)) then
c            write(*,*) 'WARNING...part OB...',ijkp(l,3)
c         endif

c         if ((ijkp(l,1) .gt. nx) .or. (ijkp(l,1) .lt. 1)) then
c            write(*,*) 'WARNING...part OB...',ijkp(l,:)
c         endif


 10      continue


c      do l = 1,Ni_tot
c         zindx(:) = minloc(abs(xp(l,3) - qz(:)))      
c         ijkp(l,3) = zindx(1)
cc         write(*,*) 'zindx...',ijkp(l,3)
c      enddo

c      where (xp(:,1) .ge. qx(nx))
c         ijkp(:,1) = 1
c         wquad(:,1) = 0.0
c         xp(:,1) = qx(1) + ( xp(:,1) - qx(nx) )
c      endwhere

c      do l = Ni_tot+1,Ni_tot+dNi_sw
cc      do l = 1,Ni_tot_sw
cc         if (xp(l,1) .le. qx(1)) then
cc            write(*,*) 'x boundary 1...',l,ijkp(l,1),ijkp(l,2),ijkp(l,3)
cc            do m=1,3
cc               vp1(l,m) = vp(l,m)
cc               input_E = input_E - 
cc     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
cc               input_p(m) = input_p(m) - m_arr(l)*vp(l,m) / beta
cc            enddo
c            ijkp(l,1) = nx
c            wquad(l,1) = -1.0
c            xp(l,1) = qx(nx) - pad_ranf()*dx/2 !(qx(1) - xp(l,1))
c            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny)-qy(1))
cc            xp(l,2) = pad_ranf()*qy(ny)
c            ijkp(l,2) = nint(xp(l,2)/dy)
c            xp(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
cc            xp(l,3) = pad_ranf()*qz(nz)
cc            ijkp(l,3) = nint(xp(l,3)/dz_grid(nz/2))
c            k=1
c            do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c               ijkp(l,3) = k    !grid
c               k=k+1
c 50         continue
c            k=ijkp(l,3)
c            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c               ijkp(l,3) = k+1
c            endif
c            phi = 2.0*pi*pad_ranf()
cc            vp(l,1) = -vsw
c            flg = 0
c            do 30 while (flg .eq. 0)
c               theta = pi*pad_ranf()
c               f = sin(theta)
c               rnd = pad_ranf()
c               if (f .ge. rnd) flg = 1
c 30         continue
c            flg = 0
c            do 40 while (flg .eq. 0)
c               v = (100*pad_ranf())
cc               f = (vth**2/exp(1.0))*v**2*exp(-(v)**2 / vth**2)
c               f = exp(-(v)**2 / 19.3**2)
c               rnd = pad_ranf()
c               if (f .ge. rnd) then 
c                  flg = 1
c                  vp(l,1) = -vsw + v*cos(phi)*sin(theta)
c                  vp(l,2) = v*sin(phi)*sin(theta)
c                  vp(l,3) = v*cos(theta)
c               endif
c 40         continue
c            m_arr(l) = mproton
c            mrat(l) = 1.0
c            do m=1,3
c               vp1(l,m) = vp(l,m)
c               input_E = input_E + 
c     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
c               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
c            enddo
cc            write(*,*) 'x boundary 2...',l,ijkp(l,1),ijkp(l,2),ijkp(l,3) 
cc         endif
c      enddo


c      Ni_tot = Ni_tot + dNi_sw

c      write(*,*) 'Ni_tot + dNi_sw....',Ni_tot

c      in_bounds(1:Ni_tot) = .true.
c      in_bounds(Ni_tot+1:) = .false.

cc      where(xp(Ni_tot_sw+1:Ni_tot,1) .le. qx(1)) 
cc     x            in_bounds(Ni_tot_sw+1:Ni_tot) = .false.

c      where(xp(1:Ni_tot,1) .le. qx(1))
c     x     in_bounds(1:Ni_tot) = .false.

c      Ni_tot_in = count(in_bounds)

cc      write(*,*) 'more outflow particles...',Ni_tot,Ni_tot_in,
cc     x       Ni_tot-Ni_tot_in,qx(1)
      
c      do m = 1,3 
c        xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), in_bounds(1:Ni_tot))
c        vp(1:Ni_tot_in,m) = pack(vp(1:Ni_tot,m), in_bounds(1:Ni_tot))
c        vp1(1:Ni_tot_in,m) = pack(vp1(1:Ni_tot,m), in_bounds(1:Ni_tot))
c        ijkp(1:Ni_tot_in,m)=pack(ijkp(1:Ni_tot,m), in_bounds(1:Ni_tot))
c        wquad(1:Ni_tot_in,m)=pack(wquad(1:Ni_tot,m),in_bounds(1:Ni_tot))
c      enddo

c      do m = 1,8
c        wght(1:Ni_tot_in,m)=pack(wght(1:Ni_tot,m), in_bounds(1:Ni_tot))
c      enddo

c      m_arr(1:Ni_tot_in) = pack(m_arr(1:Ni_tot), in_bounds(1:Ni_tot))
c      mrat(1:Ni_tot_in) = pack(mrat(1:Ni_tot), in_bounds(1:Ni_tot))

cc remove energy

c      do l = Ni_tot_in+1,Ni_tot  
c         do m=1,3
c            input_E = input_E - 
c     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
c            input_p(m) = input_p(m) - m_arr(l)*vp(l,m) / beta
c         enddo         
c      enddo

c      Ni_tot = count(in_bounds)

c      print *,'Ni_tot in_bounds...',Ni_tot

c      where (xp(:,1) .le. qx(1)) 
c         ijkp(:,1) = nx
c         wquad(:,1) = -1.0
c         xp(:,1) = qx(nx) - (qx(1) - xp(:,1))
c         xp(:,2) = pad_ranf()*qy(ny)
c         ijkp(:,2) = ninit(xp(:,2)/dy
c         vp(:,1) = -vsw
c         vp(:,2) = 0.0
c         vp(:,3) = 0.0
c      endwhere


c y periodic
c      where (xp(:,2) .ge. qy(ny))
c         ijkp(:,2) = 1
c         wquad(:,2) = 0.0
c         xp(:,2) = qy(1) + ( xp(:,2) - qy(ny) )
c      endwhere

c      where (xp(:,2) .le. qy(1)) 
c         ijkp(:,2) = ny
c         wquad(:,2) = -1.0
c         xp(:,2) = qy(ny) - (qy(1) - xp(:,2))
c      endwhere


      where (xp(:,2) .gt. qy(ny-1))
c         ijkp(:,2) = 1
c         wquad(:,2) = 0
         xp(:,2) = qy(1) + ( xp(:,2) - qy(ny-1) )
c         ijkp(:,2) = nint(xp(:,2)/dy)
      endwhere

      where (xp(:,2) .le. qy(1)) 
c         ijkp(:,2) = ny-1
c         wquad(:,2) = -1
         xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
c         ijkp(:,2) = nint(xp(:,2)/dy)
      endwhere


c -------------------z exchange, up-----------------------------

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      where (xp(1:Ni_tot,3) .gt. qz(nz))
         in_bounds(1:Ni_tot)= .false.
         ijkp(1:Ni_tot,3) = 2
c         wquad(1:Ni_tot,3) = 0.0
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

c      write(*,*) 'Ni_out up...',Ni_out, my_rank

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      dest = nbrs(n_up)
      source = nbrs(n_down)

c      call MPI_ISEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(Ni_in, 1, MPI_INTEGER, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

      call MPI_SEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, ierr)
      call MPI_RECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, stat, ierr)
      
      allocate(in_part(Ni_in,3))
      allocate(in_ijkp(Ni_in,3))
      allocate(in_part_wght(Ni_in,8))
      allocate(in_mass(Ni_in))


      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      do m = 1,3
         out_part(1:Ni_out,m) = 
     x        pack(xp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         write(*,*) 'finished out_part...',Ni_out,Ni_tot_in,Ni_tot
c         xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo

      call pack_pd(xp, in_bounds, 3)
c      write(*,*) 'finished xp pack...'
                

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

     


c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_SEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, stat, ierr)
            
      xp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         vp(1:Ni_tot_in,m) = pack(vp(1:Ni_tot,m), in_bounds(1:Ni_tot))
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

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_SEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, stat, ierr)
      
      vp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(vp1(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         vp1(1:Ni_tot_in,m) = pack(vp1(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo

      call pack_pd(vp1, in_bounds, 3)

      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_SEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, stat, ierr)
      
       vp1(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,3
         out_ijkp(1:Ni_out,m) = 
     x          pack(ijkp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         ijkp(1:Ni_tot_in,m) = 
c     x          pack(ijkp(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo

c      call pack_pd_3(ijkp, in_bounds)

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

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_SEND(out_ijkp, 3*Ni_out, MPI_INTEGER, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(in_ijkp, 3*Ni_in, MPI_INTEGER, source, tag,
c     x     cartcomm, stat, ierr)
      
      
      ijkp(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_ijkp(:,:)


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x          pack(wquad(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         wquad(1:Ni_tot_in,m) = 
c     x          pack(wquad(1:Ni_tot,m), in_bounds(1:Ni_tot))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
c      
c      call MPI_WAITALL(2, reqs, stats, ierr)

cc      call MPI_Barrier(MPI_COMM_WORLD,ierr)

cc      call MPI_SEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
cc     x     cartcomm, ierr)
cc      call MPI_RECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
cc     x     cartcomm, stat, ierr)
      
      
c      wquad(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,8
         out_part_wght(1:Ni_out,m) = 
     x          pack(wght(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         wght(1:Ni_tot_in,m) = 
c     x          pack(wght(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo


      call pack_pd(wght, in_bounds, 8)

      call MPI_ISEND(out_part_wght, 8*Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part_wght, 8*Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_SEND(out_part_wght, 8*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(in_part_wght, 8*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, stat, ierr)
      
      
      wght(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part_wght(:,:)

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      out_beta_p(1:Ni_out) = 
     x     pack(beta_p(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      beta_p(1:Ni_tot_in) = pack(beta_p(1:Ni_tot), in_bounds(1:Ni_tot))

      call pack_pd(beta_p, in_bounds, 1)

      call MPI_ISEND(out_beta_p, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      
      beta_p(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


c      out_mass(1:Ni_out) = 
c     x     pack(m_arr(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      m_arr(1:Ni_tot_in) = pack(m_arr(1:Ni_tot), in_bounds(1:Ni_tot))

c      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

      
c      m_arr(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      out_mass(1:Ni_out) = 
     x     pack(mrat(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      mrat(1:Ni_tot_in) = pack(mrat(1:Ni_tot), in_bounds(1:Ni_tot))

      call pack_pd(mrat, in_bounds, 1)

      ! remove energy of outgoing particles
      input_E = input_E - sum(0.5*(mion/out_mass(:))*vsqrd_out(:)/
     x                        (beta*out_beta_p(:)))


      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      
c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      call MPI_SEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(in_mass, Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, stat, ierr)
      

      mrat(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      Ni_tot = Ni_tot - Ni_out + Ni_in

      ! add energy back in
      do m = 1,3
         input_E = input_E + sum(0.5*(mion/mrat(Ni_tot_in+1:Ni_tot))*
     x               (vp(Ni_tot_in+1:Ni_tot,m)*km_to_m)**2 /
     x                (beta*beta_p(Ni_tot_in+1:Ni_tot)))
      enddo
      
c      do l = Ni_tot_in+1,Ni_tot  
c         do m=1,3
c            input_E = input_E + 
c     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
cc            input_p(m) = input_p(m) - m_arr(l)*vp(l,m) / beta
c         enddo         
c      enddo


      deallocate(out_part)
      deallocate(out_ijkp)
      deallocate(out_part_wght)
      deallocate(out_mass)
      deallocate(out_beta_p)
      deallocate(vsqrd_out)

      deallocate(in_part)
      deallocate(in_ijkp)
      deallocate(in_part_wght)
      deallocate(in_mass)


c ---------z exchange down---------------------------------------

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      where (xp(:,3) .le. qz(1)) 
c         ijkp(:,3) = nz
c         wquad(:,3) = -1.0
c         xp(:,3) = qz(nz) - (qz(1) - xp(:,3))
c      endwhere

      in_bounds(1:Ni_tot) = .true.
      in_bounds(Ni_tot+1:) = .false.

      where (xp(1:Ni_tot,3) .le. qz(2))
         in_bounds(1:Ni_tot)= .false.
         ijkp(1:Ni_tot,3) = nz
c         wquad(1:Ni_tot,3) = -1.0
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

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      dest = nbrs(n_down)
      source = nbrs(n_up)
      
      call MPI_ISEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, reqs(2), ierr)


      call MPI_WAITALL(2, reqs, stats, ierr)
c      write(*,*) 'down exchange...',Ni_in,Ni_out,my_rank      


      allocate(in_part(Ni_in,3))
      allocate(in_ijkp(Ni_in,3))
      allocate(in_part_wght(Ni_in,8))
      allocate(in_mass(Ni_in))

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      do m = 1,3
         out_part(1:Ni_out,m) = 
     x          pack(xp(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c         vp(1:Ni_tot_in,m) = pack(vp(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c         vp1(1:Ni_tot_in,m) = pack(vp1(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c         ijkp(1:Ni_tot_in,m) = 
c     x          pack(ijkp(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo

c      call pack_pd_3(ijkp, in_bounds)

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


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x          pack(wquad(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         wquad(1:Ni_tot_in,m) = 
c     x          pack(wquad(1:Ni_tot,m), in_bounds(1:Ni_tot))
c      enddo

c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)
      
c      wquad(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

      do m = 1,8
         out_part_wght(1:Ni_out,m) = 
     x          pack(wght(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         wght(1:Ni_tot_in,m) = 
c     x          pack(wght(1:Ni_tot,m), in_bounds(1:Ni_tot))
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
c      beta_p(1:Ni_tot_in) = pack(beta_p(1:Ni_tot), in_bounds(1:Ni_tot))


      call pack_pd(beta_p, in_bounds, 1)

      call MPI_ISEND(out_beta_p, Ni_out, MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

      beta_p(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


c      out_mass(1:Ni_out) = 
c     x     pack(m_arr(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      m_arr(1:Ni_tot_in) = pack(m_arr(1:Ni_tot), in_bounds(1:Ni_tot))

c      ! remove energy of outgoing particles
c      input_E = input_E - sum(0.5*out_mass(:)*vsqrd_out(:)/
c     x                        (beta*out_beta_p(:)))

c      call MPI_ISEND(out_mass, Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_mass, Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      m_arr(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


      out_mass(1:Ni_out) = 
     x     pack(mrat(1:Ni_tot), .not.in_bounds(1:Ni_tot))
c      mrat(1:Ni_tot_in) = pack(mrat(1:Ni_tot), in_bounds(1:Ni_tot))

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

      ! add energy back in
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
      deallocate(vsqrd_out)

      deallocate(in_part)
      deallocate(in_ijkp)
      deallocate(in_part_wght)
      deallocate(in_mass)

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
      SUBROUTINE update_np(np)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------
c      include 'incurv.h'

      real np(nx,ny,nz)

      real volb              !cell volume times beta
c      real recvbuf(nx*ny*nz)
c      integer count
c      count = nx*ny*nz

c      real sumnp,vol

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               np(i,j,k)=0.0
 10            continue

      do 20 l=1,Ni_tot

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         if (i .lt. 1) then
            i = 1
            write(*,*) 'i index error...',ijkp(l,1)!,wquad(l,1)
         endif
         if (j .lt. 1) then 
            j = 1
            write(*,*) 'j index error...',ijkp(l,2)!,wquad(l,2)
         endif
         if (k .lt. 1) then
            k = 1
            write(*,*) 'k index error...',ijkp(l,3)!,wquad(l,3)
         endif

         ip = i+1
         jp = j+1
         kp = k+1

         if (ip .gt. nx) write(*,*) 'ip index error...',ip !ip = 2      !periodic boundary conditions
         if (jp .gt. ny) write(*,*) 'jp index error...',jp !jp = 2      !periodic boundary conditions
         if (kp .gt. nz) write(*,*) 'kp index error...',kp !kp = 2      !periodic boundary conditions

c         volb = dx*dy*(qz(k+1)-qz(k))*beta
c         volb = dx*dy*dz_cell(k)*beta
         volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))

         np(i,j,k) = np(i,j,k) + wght(l,1)*volb
         np(ip,j,k) = np(ip,j,k) + wght(l,2)*volb
         np(i,j,kp) = np(i,j,kp) + wght(l,3)*volb
         np(ip,j,kp) = np(ip,j,kp) + wght(l,4)*volb
         np(i,jp,k) = np(i,jp,k) + wght(l,5)*volb
         np(ip,jp,k) = np(ip,jp,k) + wght(l,6)*volb
         np(i,jp,kp) = np(i,jp,kp) + wght(l,7)*volb
         np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)*volb

 20      continue

c use for periodic boundary conditions
c         np(nx-1,:,:) = np(nx-1,:,:)+np(1,:,:)
         np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)
c         np(:,:,nz-1) = np(:,:,nz-1)+np(:,:,1)

c         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c         out_buf_z(:,:) = np(:,:,nz)         
         
c         dest = nbrs(n_up)
c         source = nbrs(n_down)
c         call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
c     x        cartcomm, reqs(1), ierr)
c         call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
c     x        cartcomm, reqs(2), ierr)
         
c         call MPI_WAITALL(2, reqs, stats, ierr)
c         np(:,:,2) = np(:,:,2) + in_buf_z

c         call MPI_Barrier(MPI_COMM_WORLD,ierr)

c         call MPI_ALLREDUCE(np(:,:,:),recvbuf,count,
c     x        MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

c         np(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))

c         write(*,*) 'recvbuf...',recvbuf(nx*ny+1:nx*ny+10)
c         write(*,*) 'np........',np(1:10,1,2)
         
c         write(*,*) 'np ...',np(20,20,20)         
c         do i = 1,nx
c            do j = 1,ny
c               do k = 1,nz

c                  call MPI_ALLREDUCE(np(i,j,k),recvbuf,count,
c     x                 MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c                  np(i,j,k) = recvbuf
c               enddo
c            enddo
c         enddo
c         write(*,*) 'np1...',np(20,20,20)

      call periodic_scalar(np)

      return
      end SUBROUTINE update_np
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE separate_np(np,mr)
c Weight density to eight nearest grid points.
c----------------------------------------------------------------------
c      include 'incurv.h'

      real np(nx,ny,nz)
      real mr

      real volb              !cell volume times beta

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               np(i,j,k)=0.0
 10         continue
               

      do 20 l=1,Ni_tot
         
         if (mrat(l) .eq. mr) then 
            

            i=ijkp(l,1)!+wquad(l,1)
            j=ijkp(l,2)!+wquad(l,2)
            k=ijkp(l,3)!+wquad(l,3)
            
            ip = i+1
            jp = j+1
            kp = k+1
            
            volb = dx*dy*dz_cell(k)*beta*beta_p(l)
            
            np(i,j,k) = np(i,j,k) + wght(l,1)/volb
            np(ip,j,k) = np(ip,j,k) + wght(l,2)/volb
            np(i,j,kp) = np(i,j,kp) + wght(l,3)/volb
            np(ip,j,kp) = np(ip,j,kp) + wght(l,4)/volb
            np(i,jp,k) = np(i,jp,k) + wght(l,5)/volb
            np(ip,jp,k) = np(ip,jp,k) + wght(l,6)/volb
            np(i,jp,kp) = np(i,jp,kp) + wght(l,7)/volb
            np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)/volb

         endif
         
 20   continue
      
c     use for periodic boundary conditions
c     np(nx-1,:,:) = np(nx-1,:,:)+np(1,:,:)
      np(:,ny-1,:) = np(:,ny-1,:)+np(:,1,:)
c     np(:,:,nz-1) = np(:,:,nz-1)+np(:,:,1)
      
      call periodic_scalar(np)

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

c      real recvbuf(nx*ny*nz)
c      integer count
c      count = nx*ny*nz


     

      do 10 m=1,3          !clear out temp variable ct
         do 10 i=1,nx
            do 10 j=1,ny
               do 10 k=1,nz
                  up(i,j,k,m)=0.0
                  ct(i,j,k,m)=0.0
 10            continue
               
               cnt(:,:,:) = 0.0
               
c      where (ijkp(:,1) .eq. 1)
c         wquad(:,1) = 0.0
c      endwhere
     

    
      do 20 l=1,Ni_tot


         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         volb = dx*dy*(qz(kp)-qz(k))*beta
c         volb = dx*dy*dz_cell(k)*beta
      
c         np_at_Ba = np(i,j,k)*wght(l,1) + np(ip,j,k)*wght(l,2) + 
c     x              np(i,j,kp)*wght(l,3) + np(ip,j,kp)*wght(l,4) + 
c     x              np(i,jp,k)*wght(l,5) + np(ip,jp,k)*wght(l,6) +
c     x              np(i,jp,kp)*wght(l,7) + np(ip,jp,kp)*wght(l,8)

c         if (sum(wght(l,:)) .gt. 1.01) then 
c            write(*,*) 'wght error...',sum(wght(l,:)),ijkp(l,:),my_rank
c         endif


c         nvolb = 1.0
c         if (np(i,j,k) .gt. 0.0) then
c         nvolb = np(i,j,k)*volb
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)/beta_p(l) 
         ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)/beta_p(l) 
         ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)/beta_p(l) 
         
c         endif

c         if (np(ip,j,k) .gt. 0.0) then
c         nvolb = np(ip,j,k)*volb
         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)
         ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)/beta_p(l) 
         ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)/beta_p(l) 
c         endif

c         if (np(i,j,kp) .gt. 0.0) then
c         nvolb = np(i,j,kp)*volb
         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)/beta_p(l) 
         ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)/beta_p(l) 
         ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)/beta_p(l) 
c         endif

c         if (np(ip,j,kp) .gt. 0.0) then
c         nvolb = np(ip,j,kp)*volb
         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l) 
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)/beta_p(l) 
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)/beta_p(l) 
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)/beta_p(l) 
c         endif

c         if (np(i,jp,k) .gt. 0.0) then
c         nvolb = np(i,jp,k)*volb
         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)/beta_p(l) 
         ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)/beta_p(l) 
         ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)/beta_p(l) 
c         endif

c         if (np(ip,jp,k) .gt. 0.0) then
c         nvolb = np(ip,jp,k)*volb
         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)/beta_p(l) 
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)/beta_p(l) 
c         endif

c         if (np(i,jp,kp) .gt. 0.0) then
c         nvolb = np(i,jp,kp)*volb
         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)/beta_p(l) 
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)/beta_p(l) 
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)/beta_p(l) 
c         endif

c         if (np(ip,jp,kp) .gt. 0.0) then
c         nvolb = np(ip,jp,kp)*volb
         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)/beta_p(l) 
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)/beta_p(l) 
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)/beta_p(l) 
c         endif

 20   continue





c use for periodic boundary conditions
c      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)
c      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


c      do mrnk = 0,procnum-1
c         if (my_rank .eq. mrnk) then
c            do i = 1,nx
c               do j = 1,ny
c                  do k = 1,nz
c                     if (cnt(i,j,k) .gt. 0.0) then
cc                        write(*,*) 'cnt...',cnt(i,j,k),my_rank
c                        ct(i,j,k,1) = ct(i,j,k,1)/cnt(i,j,k)
c                        ct(i,j,k,2) = ct(i,j,k,2)/cnt(i,j,k)
c                        ct(i,j,k,3) = ct(i,j,k,3)/cnt(i,j,k)
c                     endif
c                  enddo
c               enddo
c            enddo
cc         endif
cc      enddo
cc      stop

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      where (cnt(:,:,:) .gt. 0.0)

         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere

      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
c      up(:,:,2,:) = (up(:,:,2,:) + in_buf_z)/2
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)


c      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
c      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
c      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
c      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

      do 30 i=1,nx-1      !interpolate back to contravarient positions
         do 30 j=1,ny-1
            do 30 k=1,nz-1
c               up(i,j,k,1) = ct(i,j,k,1)
c               up(i,j,k,2) = ct(i,j,k,2)
c               up(i,j,k,3) = ct(i,j,k,3)

               up(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
               up(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
               up(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))

 30            continue


      call periodic(up)


c must add density contribution from processor below at qz(nz) to 
c density at qz(2)


c      write(*,*) 'up...',up(20,20,21,1),up(20,20,21,2),up(20,20,21,3)

c      up(nx,:,:,1) = -vsw
c      up(nx,:,:,2) = 0.0
c      up(nx,:,:,3) = 0.0


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end SUBROUTINE update_up
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE update_np_boundary(np)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real np(nx,ny,nz)

      integer dest, source
      real out_buf_z(nx,ny)
      real in_buf_z(nx,ny)
      integer cnt_buf_z
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      cnt_buf_z = nx*ny

      out_buf_z(:,:) = np(:,:,nz)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      np(:,:,2) = np(:,:,2) + in_buf_z
 
      call periodic_scalar(np)
     
c      out_buf_z(:,:) = np(:,:,2)         

c      dest = nbrs(n_down)
c      source = nbrs(n_up)
c      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)

c      call MPI_WAITALL(2, reqs, stats, ierr)
c      np(:,:,nz) = np(:,:,nz) + in_buf_z


      return 
      end SUBROUTINE update_np_boundary
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_temperature(xp,vp,np,temp_p)
c----------------------------------------------------------------------
c      include 'incurv.h'
      
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


c      real recvbuf(nx*ny*nz)
c      integer count
c      count = nx*ny*nz

      cnt_buf_z = nx*ny*3

      cnt(:,:,:) = 0.0

c      do 10 m=1,3          !clear out temp variable ct
c         do 10 i=1,nx
c            do 10 j=1,ny
c               do 10 k=1,nz
c                  up2(i,j,k,m)=0.0
c                  ct(i,j,k,m)=0.0
c 10               continue
      up2(:,:,:,:) = 0.0
      up_ave(:,:,:,:) = 0.0
      ct(:,:,:,:) = 0.0


      do m = 1,3 
         mvp(:,m) = vp(:,m)/sqrt(mrat(:))
      enddo


c      do m = 1,3 
c         mvp(:,m) = vp(:,m)/mrat(:)
c      enddo


      do 20 l=1,Ni_tot


         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         volb = dx*dy*dz_cell(k)*beta
      
c         if (np(i,j,k) .gt. 1e15) then
c         nvolb = np(i,j,k)*volb
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/beta_p(l)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/beta_p(l)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/beta_p(l)
c         endif

c         if (np(ip,j,k) .gt. 1e15) then
c         nvolb = np(ip,j,k)*volb
         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)         
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/beta_p(l)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/beta_p(l)
c         endif

c         if (np(i,j,kp) .gt. 1e15) then
c         nvolb = np(i,j,kp)*volb
         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/beta_p(l)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/beta_p(l)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/beta_p(l)
c         endif

c        if (np(ip,j,kp) .gt. 1e15) then
c         nvolb = np(ip,j,kp)*volb
         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/beta_p(l)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/beta_p(l)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/beta_p(l)
c         endif

c         if (np(i,jp,k) .gt. 1e15) then
c         nvolb = np(i,jp,k)*volb
         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/beta_p(l)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/beta_p(l)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/beta_p(l)
c         endif

c         if (np(ip,jp,k) .gt. 1e15) then
c         nvolb = np(ip,jp,k)*volb
         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/beta_p(l)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/beta_p(l)
c        endif

c         if (np(i,jp,kp) .gt. 1e15) then
c         nvolb = np(i,jp,kp)*volb
         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/beta_p(l)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/beta_p(l)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/beta_p(l)
c         endif

c         if (np(ip,jp,kp) .gt. 1e15) then
c         nvolb = np(ip,jp,kp)*volb
         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1)+mvp(l,1)**2*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2)+mvp(l,2)**2*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3)+mvp(l,3)**2*wght(l,8)/beta_p(l)
c         endif


 20   continue

c use for periodic boundary conditions
c      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)
c      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)
      

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
c      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
c      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
c      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

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

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

c         volb = dx*dy*dz_cell(k)*beta
      
c         if (np(i,j,k) .gt. 1e15) then
c         nvolb = np(i,j,k)*volb
         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)/beta_p(l)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/beta_p(l)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/beta_p(l)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/beta_p(l)
c         endif

c         if (np(ip,j,k) .gt. 1e15) then
c         nvolb = np(ip,j,k)*volb
         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)/beta_p(l)
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/beta_p(l)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/beta_p(l)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/beta_p(l)
c         endif

c         if (np(i,j,kp) .gt. 1e15) then
c         nvolb = np(i,j,kp)*volb
         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)/beta_p(l)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/beta_p(l)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/beta_p(l)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/beta_p(l)
c         endif

c         if (np(ip,j,kp) .gt. 1e15) then
c         nvolb = np(ip,j,kp)*volb
         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)/beta_p(l)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/beta_p(l)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/beta_p(l)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/beta_p(l)
c         endif

c         if (np(i,jp,k) .gt. 1e15) then
c         nvolb = np(i,jp,k)*volb
         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)/beta_p(l)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/beta_p(l)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/beta_p(l)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/beta_p(l)
c         endif

c         if (np(ip,jp,k) .gt. 1e15) then
c         nvolb = np(ip,jp,k)*volb
         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)/beta_p(l)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/beta_p(l)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/beta_p(l)
c         endif

c         if (np(i,jp,kp) .gt. 1e15) then
c         nvolb = np(i,jp,kp)*volb
         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)/beta_p(l)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/beta_p(l)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/beta_p(l)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/beta_p(l)
c         endif

c         if (np(ip,jp,kp) .gt. 1e15) then
c         nvolb = np(ip,jp,kp)*volb
         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)/beta_p(l)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/beta_p(l)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/beta_p(l)
c         endif


 40   continue

c use for periodic boundary conditions
c      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)
c      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      ct(:,:,2,:) = (ct(:,:,2,:) + in_buf_z)/2

      call periodic(ct)

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)
         
c      call MPI_ALLREDUCE(ct(:,:,:,1),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,1) = reshape(recvbuf,(/nx,ny,nz/))
      
c      call MPI_ALLREDUCE(ct(:,:,:,2),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,2) = reshape(recvbuf,(/nx,ny,nz/))
      
c      call MPI_ALLREDUCE(ct(:,:,:,3),recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
c      ct(:,:,:,3) = reshape(recvbuf,(/nx,ny,nz/))
      

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


c----------------------------------------------------------------------
      SUBROUTINE separate_temp(vp,temp_p,mr)
c----------------------------------------------------------------------
c      include 'incurv.h'
      
      real vp(Ni_max,3)
      real temp_p(nx,ny,nz)
      real mr
      real up_ave(nx,ny,nz,3),
     x     up2(nx,ny,nz,3)


      integer dest, source
      real out_buf_z(nx,ny,3)
      real in_buf_z(nx,ny,3)
      integer cnt_buf_z
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)

      real cnt(nx,ny,nz)

      real mvp(Ni_max,3)

      cnt_buf_z = nx*ny*3

      cnt(:,:,:) = 0.0

      up2(:,:,:,:) = 0.0
      up_ave(:,:,:,:) = 0.0
      ct(:,:,:,:) = 0.0


      do m = 1,3 
         mvp(:,m) = vp(:,m)/sqrt(mrat(:))
      enddo


      do 20 l=1,Ni_tot

         if (mrat(l) .eq. mr) then

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)

         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)         
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)

         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)

         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)

         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)

         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)

         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)

         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)**2*wght(l,8)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)**2*wght(l,8)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)**2*wght(l,8)

         endif

 20   continue

c use for periodic boundary conditions
c      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)
c      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)


      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
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

         if (mrat(l) .eq. mr) then

         i=ijkp(l,1)!+wquad(l,1)
         j=ijkp(l,2)!+wquad(l,2)
         k=ijkp(l,3)!+wquad(l,3)

         ip = i+1
         jp = j+1
         kp = k+1

         cnt(i,j,k) = cnt(i,j,k) + wght(l,1)
         ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)
         ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)
         ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)

         cnt(ip,j,k) = cnt(ip,j,k) + wght(l,2)
         ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)
         ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)
         ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)

         cnt(i,j,kp) = cnt(i,j,kp) + wght(l,3)
         ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)
         ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)
         ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)

         cnt(ip,j,kp) = cnt(ip,j,kp) + wght(l,4)
         ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)
         ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)
         ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)

         cnt(i,jp,k) = cnt(i,jp,k) + wght(l,5)
         ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)
         ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)
         ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)

         cnt(ip,jp,k) = cnt(ip,jp,k) + wght(l,6)
         ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)
         ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)
         ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)

         cnt(i,jp,kp) = cnt(i,jp,kp) + wght(l,7)
         ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)
         ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)
         ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)

         cnt(ip,jp,kp) = cnt(ip,jp,kp) + wght(l,8)
         ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)
         ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)
         ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)

         endif

 40   continue

c use for periodic boundary conditions
c      ct(nx-1,:,:,:) = ct(nx-1,:,:,:)+ct(1,:,:,:)
      ct(:,ny-1,:,:) = ct(:,ny-1,:,:)+ct(:,1,:,:)
      cnt(:,ny-1,:) = cnt(:,ny-1,:)+cnt(:,1,:)
c      ct(:,:,nz-1,:) = ct(:,:,nz-1,:)+ct(:,:,1,:)


      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      where (cnt(:,:,:) .gt. 0.0)
         ct(:,:,:,1) = ct(:,:,:,1)/cnt(:,:,:)
         ct(:,:,:,2) = ct(:,:,:,2)/cnt(:,:,:)
         ct(:,:,:,3) = ct(:,:,:,3)/cnt(:,:,:)
      endwhere


      out_buf_z(:,:,:) = ct(:,:,nz,:)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
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
      end SUBROUTINE separate_temp
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


cc----------------------------------------------------------------------
c      SUBROUTINE get_np3(np,np3)
cc----------------------------------------------------------------------

c      include 'incurv.h'
c      real np(nx,ny,nz)
cc      real nf(nx,ny,nz)
c      real np3(nx,ny,nz,3)

c      real nfp(nx,ny,nz)

c      nfp = np

c      do i = 2,nx-1
c         do j = 2,ny-1
c            do k = 2,nz-1
c               np3(i,j,k,1) = 0.5*(nfp(i,j,k)+nfp(i+1,j,k))
c               np3(i,j,k,2) = 0.5*(nfp(i,j,k)+nfp(i,j+1,k))
c               np3(i,j,k,3) = 0.5*(nfp(i,j,k)+nfp(i,j,k+1))
c            enddo
c         enddo
c      enddo

c      call periodic(np3)

c      return
c      end SUBROUTINE get_np3
cc----------------------------------------------------------------------


      end MODULE gutsp_dd
















































