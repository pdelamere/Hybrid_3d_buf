      MODULE grid_interp

      USE global
      USE dimensions
      USE boundary

      contains

c----------------------------------------------------------------------
      SUBROUTINE edge_to_center(bt,btc)
c----------------------------------------------------------------------
      real bt(nx,ny,nz,3),   !main cell edge
     x     btc(nx,ny,nz,3)   !main cell center

c      real zrat           !ratio for doing linear interpolation
c                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b1,b2
      real btmf(nx,ny,nz,3)

      call periodic(bt)

      call edge_to_face(bt,btmf)
      call face_to_center(btmf,btc)

      call periodic(btc)
      
      return
      end SUBROUTINE edge_to_center
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE edge_to_face(bt,btmf)
c----------------------------------------------------------------------
      real bt(nx,ny,nz,3),   !main cell edge
     x     btmf(nx,ny,nz,3)
      real btc(nx,ny,nz,3)  !main cell center

c      real zrat           !ratio for doing linear interpolation
c                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b1,b2

      call periodic(bt)

      do 10 k=2,nz
         do 10 j=2,ny
            do 10 i=2,nx

c               ip = i+1
c               jp = j+1
c               kp = k+1
               im = i-1
               jm = j-1
               km = k-1

c               b2 = 0.5*(bt(i,j,k,1)+bt(i,jm,k,1))
c               b1 = 0.5*(bt(i,j,km,1)+bt(i,jm,km,1))
               b2 = bt(i,jm,k,1) + yrat(j)*(bt(i,j,k,1)-bt(i,jm,k,1))
               b1 = bt(i,jm,km,1) + yrat(j)*(bt(i,j,km,1)-bt(i,jm,km,1))
 
               btc(i,j,k,1) = b1 + zrat(k)*(b2-b1)

c               b2 = 0.5*(bt(i,j,k,2)+bt(im,j,k,2))
c               b1 = 0.5*(bt(i,j,km,2)+bt(im,j,km,2))
 
               b2 = bt(im,j,k,2) + xrat(i)*(bt(i,j,k,2)-bt(im,j,k,2))
               b1 = bt(im,j,km,2) + xrat(i)*(bt(i,j,km,2)-bt(im,j,km,2))

               btc(i,j,k,2) = b1 + zrat(k)*(b2-b1)

               b2 = bt(i,jm,k,3) + yrat(j)*(bt(i,j,k,3)-bt(i,jm,k,3))
               b1 = bt(im,jm,k,3) + yrat(j)*(bt(im,j,k,3)-bt(im,jm,k,3))

               btc(i,j,k,3) = b1 + xrat(i)*(b2-b1)

c               btc(i,j,k,3) = 0.25*(bt(i,j,k,3)+bt(im,j,k,3)+ 
c     x              bt(im,jm,k,3) + bt(i,jm,k,3))

      
 10         continue

      call periodic(btc)

      do 20 k=2,nz-1
         do 20 j=2,ny-1
            do 20 i=2,nx-1

               btmf(i,j,k,1) = 0.5*(btc(i,j,k,1)+btc(i+1,j,k,1))
               btmf(i,j,k,2) = 0.5*(btc(i,j,k,2)+btc(i,j+1,k,2))
               btmf(i,j,k,3) = 0.5*(btc(i,j,k,3)+btc(i,j,k+1,3))
      
 20         continue
      
      return
      end SUBROUTINE edge_to_face
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE face_to_center(v,vc)
c----------------------------------------------------------------------
      !!include 'incurv.h'

      real v(nx,ny,nz,3)        !vector at contravarient position
      real vc(nx,ny,nz,3)       !vector at cell center
c      real zfrc(nz)             !0.5*dz_grid(k)/dz_cell(k)

c      do 5 k=1,nz
c         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
c 5       continue

      call periodic(v)

      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz

               im = i-1     
               jm = j-1     
               km = k-1

c               if (im .lt. 1) then im = nx-1
c               if (jm .lt. 1) then jm = ny-1
c               if (km .lt. 1) then km = nz-1

c               vc(i,j,k,1) = 0.5*(v(i,j,k,1) + v(im,j,k,1))
               vc(i,j,k,1) = xrat(i)*(v(i,j,k,1) - v(im,j,k,1)) + 
     x                                v(im,j,k,1)

c               vc(i,j,k,2) = 0.5*(v(i,j,k,2) + v(i,jm,k,2))
               vc(i,j,k,2) = yrat(j)*(v(i,j,k,2) - v(i,jm,k,2)) + 
     x                                v(i,jm,k,2)
               vc(i,j,k,3) = zrat(k)*(v(i,j,k,3) - v(i,j,km,3)) + 
     x                                v(i,j,km,3)
 10            continue

      call periodic(vc)

      return
      end SUBROUTINE face_to_center
c----------------------------------------------------------------------

      end MODULE grid_interp
