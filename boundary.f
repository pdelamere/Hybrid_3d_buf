      MODULE boundary

      USE global
c      USE chem_rates
     
      contains


c---------------------------------------------------------------------
      SUBROUTINE periodic(b)
c---------------------------------------------------------------------
CVD$F VECTOR
c      include 'incurv.h'

      real b(nx,ny,nz,3)

      integer stats(MPI_STATUS_SIZE,2), reqs(2)
      integer dest, source
      real out_buf_z(nx,ny,3)
      real in_buf_z(nx,ny,3)
      integer cnt_buf_z
      cnt_buf_z = nx*ny*3


c x direction

      do 10 j=1,ny
         do 10 k=1,nz
            do 10 m=1,3
c     b(1,j,k,m) = b(nx-1,j,k,m)
c     b(nx,j,k,m) = b(2,j,k,m)
               b(1,j,k,m) = b(2,j,k,m)
c     b(nx,j,k,m) = b(2,j,k,m)
c     b(nx-1,j,k,m) = b(nx-2,j,k,m)
               b(nx,j,k,m) = b(nx-1,j,k,m)
 10   continue
      
      
c     y direction
      
      do 20 i=1,nx
         do 20 k=1,nz
            do 20 m=1,3
               b(i,1,k,m) = b(i,ny-1,k,m)
               b(i,ny,k,m) = b(i,2,k,m)
 20         continue
            
c     z direction
            
c      do 30 i=1,nx
c         do 30 j=1,ny
c            do 30 m=1,3
c               b(i,j,1,m) = b(i,j,nz-1,m)
c               b(i,j,nz,m) = b(i,j,2,m)
c 30   continue
            
                        
      out_buf_z(:,:,:) = b(:,:,nz-1,:)         
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,1,:) = in_buf_z
      
      
      out_buf_z(:,:,:) = b(:,:,2,:)         
      
      dest = down_proc
      source = up_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,nz,:) = in_buf_z


      return
      end SUBROUTINE periodic
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      SUBROUTINE x_boundary_scalar(b, us)
c---------------------------------------------------------------------
c upstream is fixed
      b(1,:,:,:) = us
c downstream is extrap
      b(nx,:,:,:) = b(nx-1,:,:,:)
      end SUBROUTINE
c---------------------------------------------------------------------
      SUBROUTINE boundary_scalar(b, us)
c---------------------------------------------------------------------
      call x_boundary_scalar(b,us)
      call periodic_scalar(b)
      end SUBROUTINE
c---------------------------------------------------------------------

      SUBROUTINE periodic_scalar(b)
c---------------------------------------------------------------------
c      include 'incurv.h'

      real b(nx,ny,nz)

      integer stats(MPI_STATUS_SIZE,2), reqs(2)
      integer dest, source
      real out_buf_z(nx,ny)
      real in_buf_z(nx,ny)
      integer cnt_buf_z
      cnt_buf_z = nx*ny

c x surfaces      !periodic
      do 10 j=1,ny
         do 10 k=1,nz
c            do 10 m=1,3
               b(1,j,k) = b(2,j,k)       !tangential
               b(nx,j,k) = b(nx-1,j,k)
 10         continue


c y surfaces
       do 20 i=1,nx
          do 20 k=1,nz
c             do 20 m=1,3
                b(i,1,k) = b(i,ny-1,k)      !tangential
                b(i,ny,k) = b(i,2,k)
 20          continue

c z surfaces
c       do 30 i=1,nx
c          do 30 j=1,ny
cc             do 30 m=1,3
c                b(i,j,1) = b(i,j,nz-1)      !tangential
c                b(i,j,nz) = b(i,j,2)
c 30          continue

      out_buf_z(:,:) = b(:,:,nz-1)         

      dest = up_proc
      source = down_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,1) = in_buf_z
      
      out_buf_z(:,:) = b(:,:,2)         

      dest = down_proc
      source = up_proc
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,nz) = in_buf_z
      
      return
      end SUBROUTINE periodic_scalar
c---------------------------------------------------------------------



c---------------------------------------------------------------------
      SUBROUTINE fix_normal_b(b)
c---------------------------------------------------------------------
c      include 'incurv.h'

      real b(nx,ny,nz,3)

c normal x components
      do 10 j=2,ny-1
         do 10 k=2,nz-1
            b(2,j,k,1) = dx*(b(2,j+1,k,2) - b(2,j,k,2))/dy +
     x                    dx*(b(2,j,k+1,3) - b(2,j,k,3))/dz_grid(k) +
     x                    b(3,j,k,1)
c            b(1,j,k,1) = b(2,j,k,1)

            b(nx-1,j,k,1) = b(nx-2,j,k,1) -
     x               dx*(b(nx-2,j+1,k,2) - b(nx-2,j,k,2))/dy -
     x               dx*(b(nx-2,j,k+1,3) - b(nx-2,j,k,3))/dz_grid(k)
 10         continue

c normal y components
c      do 20 i=2,nx-1
c         do 20 k=2,nz-1
c            b(i,2,k,2) = dy*(b(i+1,2,k,1) - b(i,2,k,1))/dx + 
c     x                    dy*(b(i,2,k+1,3) - b(i,2,k,3))/dz_grid(k) + 
c     x                    b(i,3,k,2)

c            b(i,ny,k,2) = b(i,ny-1,k,2) -
c     x                  dy*(b(i+1,ny-1,k,1) - b(i,ny-1,k,1))/dx -
c     x                  dy*(b(i,ny-1,k+1,3) - b(i,ny-1,k,3))/dz_grid(k)
c 20         continue

c normal z components
c      do 30 i=2,nx-1
c         do 30 j=2,ny-1
c            b(i,j,2,3) = dz_grid(2)*(b(i+1,j,2,1) - b(i,j,2,1))/dx + 
c     x                   dz_grid(2)*(b(i,j+1,2,2) - b(i,j,2,2))/dy +
c     x                   b(i,j,3,3)

c            b(i,j,nz,3) = b(i,j,nz-1,3) -
c     x               dz_grid(nz)*(b(i+1,j,nz-1,1) - b(i,j,nz-1,1))/dx +
c     x               dz_grid(nz)*(b(i,j+1,nz-1,2) - b(i,j,nz-1,2))/dy

c 30         continue


      return
      end SUBROUTINE fix_normal_b
c---------------------------------------------------------------------

cc---------------------------------------------------------------------
c      SUBROUTINE smooth_boundary(b)
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real b(nx,ny,nz,3)

cc x surfaces
c      do 10 j=2,ny-1
c         do 10 k=2,nz-1
c            do 10 m=1,3
c               b(1,j,k,m) = b(2,j,k,m)       !tangential
c               b(nx,j,k,m) = b(nx-1,j,k,m)
c 10         continue

cc y surfaces
c       do 20 i=1,nx
c          do 20 k=1,nz
c             do 20 m=1,3
c                b(i,1,k,m) = b(i,2,k,m)      !tangential
c                b(i,ny,k,m) = b(i,ny-1,k,m)
c 20          continue

cc z surfaces
c       do 30 i=1,nx
c          do 30 j=1,ny
c             do 30 m=1,3
c                b(i,j,1,m) = b(i,j,2,m)      !tangential
c                b(i,j,nz,m) = b(i,j,nz-1,m)
c 30          continue

c      return
c      end
cc---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE fix_tangential_E(E)
c---------------------------------------------------------------------
c      include 'incurv.h'

      real E(nx,ny,nz,3)

c     i = 2 & i = nx
      do 10 j=2,ny     !periodic boundary conditions
         do 10 k=2,nz
cc            E(2,j,k,1) = E(nx-1,j,k,1)  !normal component
cc            E(2,j,k,2) = E(nx-1,j,k,2)
cc            E(2,j,k,3) = E(nx-1,j,k,3)
c            E(nx,j,k,1) = E(3,j,k,1)  !normal component
c            E(nx,j,k,2) = E(3,j,k,2)
c            E(nx,j,k,3) = E(3,j,k,3)
c            E(nx,j,k,1) =   !normal component
            E(nx,j,k,2) = -2.3
            E(nx,j,k,3) = 0.0

 10         continue

c      write(*,*) 'E bnd...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)

cc     j = 2 & j = ny
c      do 20 i=2,nx
c         do 20 k=2,nz
c            E(i,2,k,1) = E(i,3,k,1)
cc            E(i,2,k,2) = E(i,3,k,2)   !normal component
c            E(i,2,k,3) = E(i,3,k,3)
c            E(i,ny,k,1) = E(i,ny-1,k,1)
c            E(i,ny,k,2) = E(i,ny-1,k,2)  !normal component
c            E(i,ny,k,3) = E(i,ny-1,k,3)
c 20         continue

cc     k = 2 & k = nz
c      do 30 i=2,nx-1
c         do 30 j=2,ny-1
cc            E(i,j,2,1) = E(i,j,3,1)
cc            E(i,j,2,2) = E(i,j,3,2)
cc            E(i,j,nz,1) = E(i,j,nz-1,1)
cc            E(i,j,nz,2) = E(i,j,nz-1,2)
c            E(i,j,2,1) = 0.0
c            E(i,j,2,2) = 0.0
cc            E(i,j,2,3) = 0.0   !normal component
c            E(i,j,nz-1,1) = 0.0
c            E(i,j,nz-1,2) = 0.0
cc            E(i,j,nz,3) = 0.0  !normal component
c 30         continue

c      call periodic(E)

      return
      end SUBROUTINE fix_tangential_E
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE boundaries(b)
c---------------------------------------------------------------------
c      include 'incurv.h'

      real b(nx,ny,nz,3)

c      call periodic(aa)
c      call extrapolated(aa)
c      call edges_corners(aa)

c      call normal_B(b)
c      call tangent_B_zero(b)
cc      call edges_corners_2(b)
c      call edges_corners_1(b)

c      call copy_to_boundary(b)
c      call fix_normal_b(b)
c      call smooth_boundary(b)

      return
      end SUBROUTINE boundaries
c---------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE Neut_Center(cx,cy,cz)
c Locates the cartisian coords of the center of the neutral cloud
c at a given time t.
c----------------------------------------------------------------------
CVD$F NOVECTOR
c      include 'incurv.h'

c      t = m*dt + tstart          !0.2 reflects canister evacuation time
c      cx = qx(ri) + vsat*(t-0.2) !release point + cloud expansion
c      cx = qx(ri) + vsat*t       !release point 
c      cy = qy(rj) + dy/1e10      !second term to avoid division
c      cz = qz(rk)                !by zero.  That is to avoid

      x0 = dx_grid(nx/2+ri0)/2
      y0 = dy_grid(ny/2)/2
      z0 = dz_grid(nz/2)/2
      
      cx = qx(int(nx/2+ri0)) + x0
      cy = qy(ny/2) + y0
c      cz = qz(rk/2) + io_proc*qz(nz) + z0 !defines second proc from bottom
      cz = procnum*qz(nz-1)/2 + z0
c      cz = qz(1) + io_proc*qz(nz-1) + z0 !defines second proc from bottom
                                    !in global coordinates


                                 !centering the sat track on 
                                 !whole grid points, otherwise
                                 !danger of r = 0.
      return
      end SUBROUTINE Neut_Center
c----------------------------------------------------------------------




      end MODULE boundary















      


