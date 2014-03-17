c---------------------------------------------------------------------
      SUBROUTINE periodic(b)
c---------------------------------------------------------------------
CVD$F VECTOR
      include 'incurv.h'

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
            
      do 30 i=1,nx
         do 30 j=1,ny
            do 30 m=1,3
               b(i,j,1,m) = b(i,j,nz-1,m)
               b(i,j,nz,m) = b(i,j,2,m)
 30   continue
            
                        
      out_buf_z(:,:,:) = b(:,:,nz-1,:)         
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      dest = nbrs(n_up)
      source = nbrs(n_down)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,1,:) = in_buf_z
      
      
      out_buf_z(:,:,:) = b(:,:,2,:)         
      
      dest = nbrs(n_down)
      source = nbrs(n_up)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,nz,:) = in_buf_z


      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE obstacle_boundary(E)
c---------------------------------------------------------------------
      include 'incurv.h'

      real E(nx,ny,nz,3)
      real r,cx,cy,cz

      call Neut_Center(cx,cy,cz)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               r = sqrt((qx(i)-cx)**2 + (qy(j)-cy)**2 + (gz(k)-cz)**2)
               if (r .le. RIo) then 
                  do m = 1,3
                     E(i,j,k,m) = 0.0 !E(i,j,k,m)*exp(-1.)
                  enddo
               endif
 10   continue
            
      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE obstacle_boundary_nu(nu)
c---------------------------------------------------------------------
      include 'incurv.h'

      real nu(nx,ny,nz)
      real r,cx,cy,cz

      call Neut_Center(cx,cy,cz)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               r = sqrt((qx(i)-cx)**2 + (qy(j)-cy)**2 + (gz(k)-cz)**2)
               nu(i,j,k) = nu(i,j,k) + 10.0*nu_init*exp(-r**2/RIo**2)
 10   continue
            
      return
      end
c---------------------------------------------------------------------



c---------------------------------------------------------------------
      SUBROUTINE obstacle_boundary_B(b0,b1p2)
c---------------------------------------------------------------------
      include 'incurv.h'

      real b0(nx,ny,nz,3)
      real b1p2(nx,ny,nz,3)
      real r,cx,cy,cz

      call Neut_Center(cx,cy,cz)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               r = sqrt((qx(i)-cx)**2 + (qy(j)-cy)**2 + (gz(k)-cz)**2)
               if (r .lt. 0.5*RIo) then
                  do m = 1,3
                     b1p2(i,j,k,m) = 0.0 !b0(i,j,k,m)
                  enddo
               endif
c               if (r .ge. 0.5*RIo) then
c                  do m = 1,3
c                     b1p2(i,j,k,m) = b1p2(i,j,k,m) - 
c     x                    b1p2(i,j,k,m)*exp(-(r-RIo)**2/(4*dx)**2)
c                  enddo
c               endif

 10   continue
            
      return
      end
c---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE damp(aa)
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real aa(nx,ny,nz,3)

c      do 10 j=1,ny
c         do 10 k=1,nz
c            do 10 m=1,3

c               aa(2,j,k,m) = aa(2,j,k,m)*0.8 
c               aa(3,j,k,m) = aa(3,j,k,m)*0.85
c               aa(4,j,k,m) = aa(4,j,k,m)*0.9
c               aa(5,j,k,m) = aa(5,j,k,m)*0.95
c               aa(nx,j,k,m) = aa(nx,j,k,m)*0.8
c               aa(nx-1,j,k,m) = aa(nx-1,j,k,m)*0.85 
c               aa(nx-2,j,k,m) = aa(nx-2,j,k,m)*0.9 
c               aa(nx-3,j,k,m) = aa(nx-3,j,k,m)*0.95 

c 10            continue


c      do 20 i=2,nx-1
c         do 20 k=1,nz
c            do 20 m=1,3

c               aa(i,2,k,m) = aa(i,2,k,m)*0.8 
c               aa(i,3,k,m) = aa(i,3,k,m)*0.85
c               aa(i,4,k,m) = aa(i,4,k,m)*0.9
c               aa(i,5,k,m) = aa(i,5,k,m)*0.95
c               aa(i,ny,k,m) = aa(i,ny,k,m)*0.8
c               aa(i,ny-1,k,m) = aa(i,ny-1,k,m)*0.85 
c               aa(i,ny-2,k,m) = aa(i,ny-2,k,m)*0.9 
c               aa(i,ny-3,k,m) = aa(i,ny-3,k,m)*0.95 

c 20            continue


c      do 30 i=2,nx-1
c         do 30 j=2,ny-1
c            do 30 m=1,3

c               aa(i,j,2,m) = aa(i,j,2,m)*0.8 
c               aa(i,j,3,m) = aa(i,j,3,m)*0.85
c               aa(i,j,4,m) = aa(i,j,4,m)*0.9
c               aa(i,j,5,m) = aa(i,j,5,m)*0.95
c               aa(i,j,nz-1,m) = aa(i,j,nz-1,m)*0.8
c               aa(i,j,nz-2,m) = aa(i,j,nz-2,m)*0.85 
c               aa(i,j,nz-3,m) = aa(i,j,nz-3,m)*0.9 
c               aa(i,j,nz-4,m) = aa(i,j,nz-4,m)*0.95 

c 30            continue

c      return
c      end
cc---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE extrapolated(aa)
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real aa(nx,ny,nz,3)

c      do 40 j=2,ny-1 
c         do 40 k=2,nz-1
c            do 40 m=1,3
c               aa(1,j,k,m) = 3.0*aa(2,j,k,m) - 3.0*aa(3,j,k,m)  
c     x                          + aa(4,j,k,m)
c               aa(nx,j,k,m) = 3.0*aa(nx-1,j,k,m) - 3.0*aa(nx-2,j,k,m)  
c     x                           + aa(nx-3,j,k,m)
c 40            continue


c      do 50 i=2,nx-1
c         do 50 k=2,nz-1
c            do 50 m=1,3
c               aa(i,1,k,m) = -(3.0*aa(i,2,k,m) - 3.0*aa(i,3,k,m) 
c     x                          + aa(i,4,k,m))
c               aa(i,ny,k,m) = 3.0*aa(i,ny-1,k,m) - 3.0*aa(i,ny-2,k,m)
c     x                           + aa(i,ny-3,k,m)
c 50            continue    

c      do 60 i=2,nx-1
c         do 60 j=2,ny-1
c            do 60 m=1,3

cc               aa(i,j,1,m) = -(3.0*aa(i,j,2,m) - 3.0*aa(i,j,3,m)
cc     x                          + aa(i,j,4,m))
cc               aa(i,j,nz,m) = 3.0*aa(i,j,nz-1,m) - 3.0*aa(i,j,nz-2,m)
cc     x                           + aa(i,j,nz-3,m)
c 60            continue

c      return
c      end
cc---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE edges_corners_1(aa)
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real aa(nx,ny,nz,3)
c      real f2,f3

c      f2=(1.0/2.0)    !damp edges slightly...should be 1/2
c      f3=(1.0/3.0)    !damp corners slightly...should be 1/3

cc x edge 

c      do 10 i=2,nx-1
c         do 10 m=1,3
c            aa(i,1,1,m)   = f2*(aa(i,2,1,m)     + aa(i,1,2,m))
c            aa(i,1,nz,m)  = f2*(aa(i,2,nz,m)    + aa(i,1,nz-1,m))
c            aa(i,ny,1,m)  = f2*(aa(i,ny,2,m)    + aa(i,ny-1,1,m))
c            aa(i,ny,nz,m) = f2*(aa(i,ny,nz-1,m) + aa(i,ny-1,nz,m))
c 10      continue

cc y edge

c      do 20 j=2,ny-1
c         do 20 m=1,3
c            aa(1,j,1,m)   = f2*(aa(2,j,1,m)     + aa(1,j,2,m))
c            aa(1,j,nz,m)  = f2*(aa(2,j,nz,m)    + aa(1,j,nz-1,m))
c            aa(nx,j,1,m)  = f2*(aa(nx,j,2,m)    + aa(nx-1,j,1,m))
c            aa(nx,j,nz,m) = f2*(aa(nx,j,nz-1,m) + aa(nx-1,j,nz,m))
c 20      continue

cc z edge

c      do 30 k=2,nz-1
c         do 30 m=1,3
c            aa(1,1,k,m)   = f2*(aa(2,1,k,m)     + aa(1,2,k,m))
c            aa(1,ny,k,m)  = f2*(aa(2,ny,k,m)    + aa(1,ny-1,k,m))
c            aa(nx,1,k,m)  = f2*(aa(nx,2,k,m)    + aa(nx-1,1,k,m))
c            aa(nx,ny,k,m) = f2*(aa(nx,ny-1,k,m) + aa(nx-1,ny,k,m))
c 30      continue

cc corners

c      do 40 m=1,3
c         aa(1,1,1,m) = f3*(aa(1,1,2,m) + aa(1,2,1,m) + aa(2,1,1,m))
c         aa(1,1,nz,m) = f3*(aa(1,1,nz-1,m) + aa(1,2,nz,m) 
c     x                      + aa(2,1,nz,m))
c         aa(1,ny,1,m) = f3*(aa(1,ny,2,m) + aa(1,ny-1,1,m)
c     x                      + aa(2,ny,1,m))
c         aa(nx,1,1,m) = f3*(aa(nx,1,2,m) + aa(nx,2,1,m) 
c     x                      + aa(nx-1,1,1,m))
c         aa(1,ny,nz,m) = f3*(aa(1,ny,nz-1,m) + aa(1,ny-1,nz,m)
c     x                       + aa(2,ny,nz,m))
c         aa(nx,1,nz,m) = f3*(aa(nx,1,nz-1,m) + aa(nx,2,nz,m)
c     x                       + aa(nx-1,1,nz,m))
c         aa(nx,ny,1,m) = f3*(aa(nx,ny,2,m) + aa(nx,ny-1,1,m)
c     x                       + aa(nx-1,ny,1,m))
c         aa(nx,ny,nz,m) = f3*(aa(nx,ny,nz-1,m) + aa(nx,ny-1,nz,m)
c     x                        + aa(nx-1,ny,nz,m))

c 40      continue

c      return
c      end
cc---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE edges_corners_2(aa)
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real aa(nx,ny,nz,3)
c      real f3

c      f3=(1.0/3.0)

cc x edge 

c      do 10 i=3,nx-1
c         do 10 m=1,3
c            aa(i,2,2,m)   = 0.5*(aa(i,3,2,m)     + aa(i,2,3,m))
c            aa(i,2,nz,m)  = 0.5*(aa(i,3,nz,m)    + aa(i,2,nz-1,m))
c            aa(i,ny,2,m)  = 0.5*(aa(i,ny,3,m)    + aa(i,ny-1,2,m))
c            aa(i,ny,nz,m) = 0.5*(aa(i,ny,nz-1,m) + aa(i,ny-1,nz,m))
c 10      continue

cc y edge

c      do 20 j=3,ny-1
c         do 20 m=1,3
c            aa(2,j,2,m)   = 0.5*(aa(3,j,2,m)     + aa(2,j,3,m))
c            aa(2,j,nz,m)  = 0.5*(aa(3,j,nz,m)    + aa(2,j,nz-1,m))
c            aa(nx,j,2,m)  = 0.5*(aa(nx,j,3,m)    + aa(nx-1,j,2,m))
c            aa(nx,j,nz,m) = 0.5*(aa(nx,j,nz-1,m) + aa(nx-1,j,nz,m))
c 20      continue

cc z edge

c      do 30 k=3,nz-1
c         do 30 m=1,3
c            aa(2,2,k,m)   = 0.5*(aa(3,2,k,m)     + aa(2,3,k,m))
c            aa(2,ny,k,m)  = 0.5*(aa(3,ny,k,m)    + aa(2,ny-1,k,m))
c            aa(nx,2,k,m)  = 0.5*(aa(nx,3,k,m)    + aa(nx-1,2,k,m))
c            aa(nx,ny,k,m) = 0.5*(aa(nx,ny-1,k,m) + aa(nx-1,ny,k,m))
c 30      continue

cc corners

c      do 40 m=1,3
c         aa(2,2,2,m) = f3*(aa(2,2,3,m) + aa(2,3,2,m) + aa(3,2,2,m))
c         aa(2,2,nz,m) = f3*(aa(2,2,nz-1,m) + aa(2,3,nz,m) 
c     x                      + aa(3,2,nz,m))
c         aa(2,ny,2,m) = f3*(aa(2,ny,3,m) + aa(2,ny-1,2,m)
c     x                      + aa(3,ny,2,m))
c         aa(nx,2,2,m) = f3*(aa(nx,2,3,m) + aa(nx,3,2,m) 
c     x                      + aa(nx-1,2,2,m))
c         aa(2,ny,nz,m) = f3*(aa(2,ny,nz-1,m) + aa(2,ny-1,nz,m)
c     x                       + aa(3,ny,nz,m))
c         aa(nx,2,nz,m) = f3*(aa(nx,2,nz-1,m) + aa(nx,3,nz,m)
c     x                       + aa(nx-1,2,nz,m))
c         aa(nx,ny,2,m) = f3*(aa(nx,ny,3,m) + aa(nx,ny-1,2,m)
c     x                       + aa(nx-1,ny,2,m))
c         aa(nx,ny,nz,m) = f3*(aa(nx,ny,nz-1,m) + aa(nx,ny-1,nz,m)
c     x                        + aa(nx-1,ny,nz,m))

c 40      continue

c      return
c      end
cc---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE normal_B(b)  !divergence B
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real b(nx,ny,nz,3)

cc copy normal components to boundaries first
cc x surfaces
c      do 4 j=2,ny-1
c         do 4 k=2,nz-1

c            b(2,j,k,1) = b(3,j,k,1)       !normal
cc            b(1,j,k,1) = b(3,j,k,1)
c            b(nx,j,k,1) = b(nx-1,j,k,1)

c 4       continue

cc y surfaces
c       do 6 i=2,nx-1
c          do 6 k=2,nz-1

c             b(i,2,k,2) = b(i,3,k,2)      !normal
cc             b(i,1,k,2) = b(i,3,k,2)
c             b(i,ny,k,2) = b(i,ny-1,k,2)

c 6        continue

cc z surfaces
c       do 8 i=2,nx-1
c          do 8 j=2,ny-1

c             b(i,j,2,3) = b(i,j,3,3)      !normal
cc             b(i,j,1,3) = b(i,j,3,3)
c             b(i,j,nz,3) = b(i,j,nz-1,3)

c 8        continue

c      call edges_corners_2(b)

cc normal x components
c      do 10 j=2,ny-1
c         do 10 k=2,nz-1
c            b(2,j,k,1) = dx*(b(2,j+1,k,2) - b(2,j,k,2))/dy +
c     x                    dx*(b(2,j,k+1,3) - b(2,j,k,3))/dz_grid(k) +
c     x                    b(3,j,k,1)

c            b(nx,j,k,1) = b(nx-1,j,k,1) -
c     x                  dx*(b(nx-1,j+1,k,2) - b(nx-1,j,k,2))/dy -
c     x                  dx*(b(nx-1,j,k+1,3) - b(nx-1,j,k,3))/dz_grid(k)
c 10         continue

cc normal y components
c      do 20 i=2,nx-1
c         do 20 k=2,nz-1
c            b(i,2,k,2) = dy*(b(i+1,2,k,1) - b(i,2,k,1))/dx + 
c     x                    dy*(b(i,2,k+1,3) - b(i,2,k,3))/dz_grid(k) + 
c     x                    b(i,3,k,2)

c            b(i,ny,k,2) = b(i,ny-1,k,2) -
c     x                  dy*(b(i+1,ny-1,k,1) - b(i,ny-1,k,1))/dx -
c     x                  dy*(b(i,ny-1,k+1,3) - b(i,ny-1,k,3))/dz_grid(k)
c 20         continue

cc normal z components
c      do 30 i=2,nx-1
c         do 30 j=2,ny-1
c            b(i,j,2,3) = dz_grid(2)*(b(i+1,j,2,1) - b(i,j,2,1))/dx + 
c     x                   dz_grid(2)*(b(i,j+1,2,2) - b(i,j,2,2))/dy +
c     x                   b(i,j,3,3)

c            b(i,j,nz,3) = b(i,j,nz-1,3) -
c     x               dz_grid(nz)*(b(i+1,j,nz-1,1) - b(i,j,nz-1,1))/dx -
c     x               dz_grid(nz)*(b(i,j+1,nz-1,2) - b(i,j,nz-1,2))/dy

c 30         continue

c      return 
c      end
cc---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE tangent_B_zero(b) !normal derivative = 0
cc The normal derivative of the tangential components is used to 
cc determine the tangential boundary values.  ALSO, the normal
cc derivative of the normal components temporarily sets the boundary
cc values of the normal components.  These values are undated later
cc by the requirement that divB=0.  This helps with the values on
cc corners and edges.
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real b(nx,ny,nz,3)

cc x surfaces
c      do 10 j=2,ny-1
c         do 10 k=2,nz-1

c            b(1,j,k,1) = b(2,j,k,1)       !normal
cc            b(1,j,k,1) = b(3,j,k,1)
cc            b(nx,j,k,1) = b(nx-1,j,k,1)

c            b(1,j,k,2) = b(2,j,k,2)       !tangential
c            b(1,j,k,3) = b(2,j,k,3)
c            b(nx,j,k,2) = b(nx-1,j,k,2)
c            b(nx,j,k,3) = b(nx-1,j,k,3)


c 10         continue

cc y surfaces
c       do 20 i=2,nx-1
c          do 20 k=2,nz-1

c             b(i,1,k,2) = b(i,2,k,2)      !normal
cc             b(i,1,k,2) = b(i,3,k,2)
cc             b(i,ny,k,2) = b(i,ny-1,k,2)

c             b(i,1,k,1) = b(i,2,k,1)      !tangential
c             b(i,1,k,3) = b(i,2,k,3)
c             b(i,ny,k,1) = b(i,ny-1,k,1)
c             b(i,ny,k,3) = b(i,ny-1,k,3)

c 20          continue

cc z surfaces
c       do 30 i=2,nx-1
c          do 30 j=2,ny-1

c             b(i,j,1,3) = b(i,j,2,3)      !normal
cc             b(i,j,1,3) = b(i,j,3,3)
cc             b(i,j,nz,3) = b(i,j,nz-1,3)

c             b(i,j,1,1) = b(i,j,2,1)      !tangential
c             b(i,j,1,2) = b(i,j,2,2)
c             b(i,j,nz,1) = b(i,j,nz-1,1)
c             b(i,j,nz,2) = b(i,j,nz-1,2)
             
c 30          continue


c      return
c      end
cc---------------------------------------------------------------------


cc---------------------------------------------------------------------
c      SUBROUTINE copy_to_boundary(b)
cc---------------------------------------------------------------------
c      include 'incurv.h'

c      real b(nx,ny,nz,3)

cc x surfaces      !periodic
c      do 10 j=1,ny
c         do 10 k=1,nz
cc            do 10 m=1,3
cc               b(2,j,k,1) = b(nx-1,j,k,1)       !normal
c               b(1,j,k,1) = b(nx-1,j,k,1)       !normal
c               b(nx,j,k,1) = b(2,j,k,1)

c               b(1,j,k,2) = b(nx-1,j,k,2)       !tangential
cc               b(2,j,k,2) = b(nx-1,j,k,2)       !tangential
c               b(nx,j,k,2) = b(2,j,k,2)

c               b(1,j,k,3) = b(nx-1,j,k,3)       !tangential
cc               b(2,j,k,2) = b(nx-1,j,k,2)       !tangential
c               b(nx,j,k,3) = b(2,j,k,3)
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
      SUBROUTINE periodic_scalar(b)
c---------------------------------------------------------------------
      include 'incurv.h'

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
       do 30 i=1,nx
          do 30 j=1,ny
c             do 30 m=1,3
                b(i,j,1) = b(i,j,nz-1)      !tangential
                b(i,j,nz) = b(i,j,2)
 30          continue

      out_buf_z(:,:) = b(:,:,nz-1)         

      dest = nbrs(n_up)
      source = nbrs(n_down)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,1) = in_buf_z
      
      out_buf_z(:,:) = b(:,:,2)         

      dest = nbrs(n_down)
      source = nbrs(n_up)
      call MPI_ISEND(out_buf_z, cnt_buf_z , MPI_REAL, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_buf_z, cnt_buf_z, MPI_REAL, source, tag,
     x     cartcomm, reqs(2), ierr)

      call MPI_WAITALL(2, reqs, stats, ierr)
      B(:,:,nz) = in_buf_z
      




      return
      end
c---------------------------------------------------------------------



c---------------------------------------------------------------------
      SUBROUTINE fix_normal_b(b)
c---------------------------------------------------------------------
      include 'incurv.h'

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
      end
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
      include 'incurv.h'

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
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE boundaries(b)
c---------------------------------------------------------------------
      include 'incurv.h'

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
      end
c---------------------------------------------------------------------




















