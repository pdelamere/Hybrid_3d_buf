      MODULE initial

      USE global
c      USE dimensions
c      USE inputs

      contains


c----------------------------------------------------------------------
      SUBROUTINE grd6_setup(b0,bt,b12,b1,b1p2,nu)
c----------------------------------------------------------------------
c      include 'incurv.h'

c      real eoverm
c      parameter(eoverm = q/mO)
      real mO_q

      real vol

      real b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     nu(nx,ny,nz)

      real b0r,a1,a2,omegar,kr,nfr,kdz
      real b0_1x, b0_2x, b0_1y, b0_2y
      real phi

      mO_q = mO/q

      phi = 0.0*PI/180.

      b0_1x = b0_top*eoverm*sin(phi)
      b0_2x = b0_bottom*eoverm*sin(phi)

      b0_1y = b0_top*eoverm*cos(phi)
      b0_2y = b0_bottom*eoverm*cos(phi)  

c      kr = 2.0/delz

c      a1 = kr**2*b0r/(alpha*nfr)
c      a2 = kr**2*b0r**2/(alpha*nfr)

c      omegar = 0.5*(a1 + sqrt(a1**2 + 4*a2))

      do i = 1,nx
         do j = 1,ny 
            do k = 1,nz
c               b0(i,j,k,1) = 0.5*(b0_1x+b0_2x) + 
c     x              0.5*(b0_1x-b0_2x)*tanh((qz(k)-qz(nz/2))/(Lo))
c               b0(i,j,k,2) = 0.5*(b0_1y+b0_2y) + 
c     x              0.5*(b0_1y-b0_2y)*tanh((qz(k)-qz(nz/2))/(Lo))
               b0(i,j,k,2) = 1.0*b0_top*eoverm
c               b0(i,j,k,1) = 0.0
               b0(i,j,k,3) = 0.0*b0_top*eoverm
            enddo
         enddo
      enddo
      
c      do 10 k=1,nz
c         kdz = 2.0/dz_grid(k)
cc         b0(k) = b0r
cc         b0(k) = 1.2*(0.5*(-omegar + (omegar/kdz)*
cc     x                 sqrt(kdz**2 + 4*alpha*nfr)))
c 10      continue

      do 20 i=1,nx
         do 20 j=1,ny
            do 20 k=1,nz
c               nf(i,j,k) = nfr  !/km^3
               nu(i,j,k) = nu_init
               do 20 m=1,3

                  bt(i,j,k,m) = b0(i,j,k,m)
c                  b12(i,j,k,m) = 0.0 !b0(i,j,k,m)
c                  b1(i,j,k,m) = 0.0 !b0(i,j,k,m)
c                  b1p2(i,j,k,m) = 0.0 !b0(i,j,k,m)
                  vol = dx*dy*dz_cell(k)*km_to_m**3
                  input_Eb = input_Eb + 
     x                      (vol/(2.0*mu0))*(mO_q*b0(i,j,k,m))**2 

 20            continue


c      nu(1,:,:) = nu_init*50.0
c      nu(2,:,:) = nu_init*50.0
c      nu(3,:,:) = nu_init*10.0
c      nu(4,:,:) = nu_init*5.0
c      nu(5,:,:) = nu_init*2.0

c      nu(nx,:,:) = nu_init*50.0
c      nu(nx-1,:,:) = nu_init*50.0
c      nu(nx-2,:,:) = nu_init*10.0
c      nu(nx-3,:,:) = nu_init*5.0
c      nu(nx-4,:,:) = nu_init*2.0


c set resistive boundary 

c      do 50 i=1,nx
c         do 50 k=1,nz
c            nu(i,4,k) = 5.0
c            nu(i,3,k) = 10.0
c            nu(i,2,k) = 20.0
c            nu(i,1,k) = 20.0
c            nu(i,ny,k) = 20.0
c            nu(i,ny-1,k) = 10.0
c            nu(i,ny-2,k) = 5.0
c 50         continue

      do 60 i=1,nx
       do 60 j=1,ny
          do 60 k=1,nz
             nu(i,j,k) = (q*b0_init/mproton)*
     x            exp(-(qx(nx)-qx(i))**2/(20.0*dx)**2) + nu_init
 60   continue

     
c      open(30,file='nf.dat',status='unknown',form='unformatted')
c      write(30) nz
c      write(30) nf
c      close(30)

      open(40,file='b0.dat',status='unknown',form='unformatted')
      write(40) nz
      write(40) b0
      close(40)

      return
      end SUBROUTINE grd6_setup
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE grd7()
c----------------------------------------------------------------------
c      include 'incurv.h'


      parameter(nrgrd = 0)

c      rk=nz/2
c      rj=ny/2
c      ri=nx-nx/2-20

      do 10 i=1,nx
         qx(i) = i*dx
 10            continue

      do 20 j=1,ny
         qy(j) = j*dy
 20            continue


c up from release
      do 32 k = rk,rk+nrgrd
         dz_grid(k) = delz
 32   continue
c CDIR@ NEXTSCALAR
      do 34 k = rk+nrgrd+1,nz
         dz_grid(k) = delz +
     x     0.0*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1)) 
c     x     2.0*sin((k-(rk+nrgrd+1))*0.5*pi/(nz-(rk+nrgrd+1)))**2 
c                                !dz_grid(k-1) + 0.01*delz 
 34   continue


c down from release
      do 36 k = rk-nrgrd,rk-1
         dz_grid(k) = delz
 36      continue
c CDIR@ NEXTSCALAR
      do 37 k = 1,rk-nrgrd-1
         ind = rk-nrgrd-k
         dz_grid(ind) = delz + 
     x     0.0*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
c     x     2.0*sin((rk-nrgrd-1-ind)*(-0.5*pi)/(rk-nrgrd-1))**2 
c                                !dz_grid(ind+1) + 0.01*delz
 37   continue

      qz(1) = 0.0
      do 39 k=2,nz
c         write(*,*) 'dz_grid...',k,dz_grid(k)
         qz(k) = qz(k-1)+dz_grid(k)
 39   continue

      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      do 40 k=2,nz-1
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) -
     x                ((qz(k) + qz(k-1))/2.0)
 40   continue
      
      if (cart_rank .eq. procnum-1) then
         gz(:) = qz(:)
      endif
      
c      do i = 1,procnum 
      if (cart_rank .lt. procnum-1) then
         do k = 1,nz
            i = procnum - (cart_rank+1)
            gz(k)=qz(nz)*i + ((k-1)*delz) - (i*delz)
         enddo
      endif


c      call assign('assign -F system -N ultrix f:' //'c.coord.dat')
      open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',
     x         form='unformatted')



c      open(40,file='coord.dat',status='unknown',form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) qx
      write(40) qy
      write(40) qz
      write(40) dz_grid
      write(40) dz_cell
      close(40)

      return
      end SUBROUTINE grd7
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE grd9()
c----------------------------------------------------------------------
c      include 'incurv.h'

      parameter(nrgrd = 40)

      rk=(nz/2)! - 35
      rj=ny/2
      ri=nx/2

      do 10 i=1,nx
         qx(i) = i*dx
 10            continue

      do 20 j=1,ny
         qy(j) = j*dy
 20            continue

c up from release
      do 32 k = rk,rk+nrgrd
         dz_grid(k) = delz
 32   continue
c CDIR@ NEXTSCALAR
      do 34 k = rk+nrgrd+1,nz
         dz_grid(k) = delz +
     x     0.0*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1)) 
c     x     2.0*sin((k-(rk+nrgrd+1))*0.5*pi/(nz-(rk+nrgrd+1)))**2 
c                                !dz_grid(k-1) + 0.01*delz 
 34   continue

c down from release
      do 36 k = rk-nrgrd,rk-1
         dz_grid(k) = delz
 36      continue
c CDIR@ NEXTSCALAR
      do 37 k = 1,rk-nrgrd-1
         ind = rk-nrgrd-k
         dz_grid(ind) = delz + 
     x     0.0*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
c     x     2.0*sin((rk-nrgrd-1-ind)*(-0.5*pi)/(rk-nrgrd-1))**2 
c                                !dz_grid(ind+1) + 0.01*delz
 37   continue

      qz(1) = delz
      do 39 k=2,nz
c         write(*,*) 'dz_grid...',k,dz_grid(k)
         qz(k) = qz(k-1)+dz_grid(k)
 39   continue

      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      do 40 k=2,nz-1
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) -
     x                ((qz(k) + qz(k-1))/2.0)
 40            continue


c      call assign('assign -F system -N ultrix f:' //'c.coord.dat')
      open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',
     x         form='unformatted')

c      open(40,file='coord.dat',status='unknown',form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) qx
      write(40) qy
      write(40) qz
      write(40) dz_grid
      write(40) dz_cell
      close(40)

      return
      end SUBROUTINE grd9
c----------------------------------------------------------------------



      end module initial
