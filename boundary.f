      MODULE boundary

      USE global
c      USE chem_rates
     
      contains

c---------------------------------------------------------------------
      SUBROUTINE periodic(b)
c---------------------------------------------------------------------
CVD$F VECTOR
c      include 'incurv.h'

      integer ierr
      real b(nx,ny,nz,3)
      real us(ny,nz,3) ! upstream condition

      integer stats(MPI_STATUS_SIZE,2), reqs(2)
      integer dest, source
      real out_buf_z(nx,ny,3)
      real in_buf_z(nx,ny,3)
      integer cnt_buf_z
      cnt_buf_z = nx*ny*3

c x direction periodic
      b(1,:,:,:) = b(nx-1,:,:,:)
      b(nx,:,:,:) = b(2,:,:,:)
c y direction periodic
      b(:,1,:,:) = b(:,ny-1,:,:)
      b(:,ny,:,:) = b(:,2,:,:)
            
c z direction periodic, but has domain decomp
            
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

      SUBROUTINE periodic_scalar(b)
c---------------------------------------------------------------------
c      include 'incurv.h'

      integer ierr
      real b(nx,ny,nz)

      integer stats(MPI_STATUS_SIZE,2), reqs(2)
      integer dest, source
      real out_buf_z(nx,ny)
      real in_buf_z(nx,ny)
      integer cnt_buf_z
      cnt_buf_z = nx*ny

c x surfaces
      b(1,:,:) = b(nx-1,:,:)
      b(nx,:,:) = b(2,:,:)
c y surfaces
      b(:,1,:) = b(:,ny-1,:)
      b(:,ny,:) = b(:,2,:)

c z surfaces
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

c----------------------------------------------------------------------
      SUBROUTINE Neut_Center(cx,cy,cz)
c Locates the cartisian coords of the center of the neutral cloud
c at a given time t.
c----------------------------------------------------------------------
CVD$F NOVECTOR
c      include 'incurv.h'
      real cx, cy, cz
      real x0, y0, z0

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
