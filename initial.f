      MODULE initial

      USE global
      USE inputs
      USE chem_rates
      USE boundary, only: Neut_Center
c      USE dimensions
c      USE inputs

      contains

c----------------------------------------------------------------      
      subroutine initparameters()
c----------------------------------------------------------------      
      integer*4 getpid, pid
      pid = getpid()
      write(*,*) 'PID....', pid

      mion = ion_amu*1.67e-27

      lambda_i = (3e8/
     x            sqrt((nf_init/1e9)*q*q/(8.85e-12*mion)))/1e3

      dx = lambda_i*dx_frac 
      dy = lambda_i*dx_frac   !units in km
      delz = lambda_i*dx_frac          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz


      print*, 'dx....',dx

      dx_buf = 3*dx


      dt = dt_frac*mion/(q*b0_init)     !main time step
      dtsub_init = dt/ntsub !subcycle time step 

      print*, 'dt...',dt,dtsub_init

      vtop = vsw
      vbottom = -vsw

      ! Initialize upstream solar wind velocities
      vsw_us(:,:,1) = -vsw
      vsw_us(:,:,2) = 0
      vsw_us(:,:,3) = 0

      Ni_tot_0 = Ni_max*Ni_tot_frac

      write(*,*) 'Ni_tot_0...',Ni_tot_0, Ni_max,Ni_tot_frac

      mO = mion    !mass of H (kg)

      eoverm = q/mO

      m_heavy = 1.0
      np_top = nf_init
      np_bottom = nf_init/m_heavy
      f_proton_top = 0.5       !fraction relative to top
      b0_top = 1.0*b0_init
      b0_bottom = b0_init
      vth_top = vth
      vth_bottom = vth
      vth_max = 3*vth
      m_top = mion
      m_bottom = mion
      Lo = 4.0*dx             !gradient scale length of boundary

      nu_init = nu_init_frac*q*b0_init/mion

      alpha = (mu0/1e3)*q*(q/mion) !mH...determines particle scaling

      moment = surf_field * Rpluto**3

      imf_theta = (pi/180)*imf_theta
      imf_phi = (pi/180)*imf_phi

      ! Store the upstream condition of B given imf_theta and imf_phi
      b0_us(:,:,1) = cos(imf_phi)*sin(imf_theta)*b0_top*eoverm
      b0_us(:,:,2) = sin(imf_phi)*sin(imf_theta)*b0_top*eoverm
      b0_us(:,:,3) = cos(imf_theta)*b0_top*eoverm

      max_ion_density = PCE_coef*sqrt(atmosphere(r_thin))
      write(*,*) 'max_ion_density...', max_ion_density

      ! Variable needs to be initialized so it can act as an accumulator.
      bndry_Eflux = 0

      end subroutine initparameters
c----------------------------------------------------------------      

c----------------------------------------------------------------------
      SUBROUTINE grd6_setup(b0,b1,bt)
c----------------------------------------------------------------------
      real mO_q

      real vol

      real b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b1(nx,ny,nz,3)

      real b0r,a1,a2,omegar,kr,nfr,kdz
      real b0_1x, b0_2x, b0_1y, b0_2y
      real phi

      real cx, cy, cz

      mO_q = mO/q

      phi = 0.0*PI/180.

      b0_1x = b0_top*eoverm*sin(phi)
      b0_2x = b0_bottom*eoverm*sin(phi)

      b0_1y = b0_top*eoverm*cos(phi)
      b0_2y = b0_bottom*eoverm*cos(phi)  

      call Neut_Center(cx,cy,cz)

      do i = 1,nx
         do j = 1,ny 
            do k = 1,nz
               x = qx(i)-cx
               y = qy(j)-cy
               z = gz(k)-cz ! global z

               r = sqrt(x**2 + y**2 + z**2)
               ! Dipole
               b0(i,j,k,1) = (3.*moment*x*y/r**5)*eoverm
               b0(i,j,k,2) = (moment*(3.*y**2 - r**2)/r**5)*eoverm
               b0(i,j,k,3) = (3.*moment*z*y/r**5)*eoverm
               ! plus IMF
               b0(i,j,k,1) = b0(i,j,k,1)
     x                 + cos(imf_phi)*sin(imf_theta)*b0_top*eoverm
               b0(i,j,k,2) = b0(i,j,k,2)
     x                 + sin(imf_phi)*sin(imf_theta)*b0_top*eoverm
               b0(i,j,k,3) = b0(i,j,k,3)
     x                 + cos(imf_theta)*b0_top*eoverm
            enddo
         enddo
      enddo

      ! Initially there is no perturbation
      b1 = 0.0

      do 20 i=1,nx
         do 20 j=1,ny
            do 20 k=1,nz
               do 20 m=1,3

                  bt(i,j,k,m) = b0(i,j,k,m)
                  vol = dx*dy*dz_cell(k)*km_to_m**3
                  input_E = input_E + 
     x                      (vol/(2.0*mu0))*(mO_q*b0(i,j,k,m))**2 

 20            continue




     
      open(40,file='b0.dat',status='unknown',form='unformatted')
      write(40) nz
      write(40) b0
      close(40)

      return
      end SUBROUTINE grd6_setup
c----------------------------------------------------------------------

      SUBROUTINE get_nu(nu)
      real nu(nx,ny,nz)
      do 60 i=1,nx
       do 60 j=1,ny
          do 60 k=1,nz
             nu(i,j,k) = (q*b0_init/mproton)*
     x            exp(-(qx(nx)-qx(i))**2/(10.0*dx)**2) + nu_init
 60   continue
      end SUBROUTINE


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
      SUBROUTINE grd8()
c----------------------------------------------------------------------
c      include 'incurv.h'

      integer nrgrd

      real zplus,zminus,xplus, xminus, yplus, yminus
      real zsf,xsf,ysf

      rk=(nz/2)! - 35
      rj=ny/2
      ri=nx/2


c==============stretch y direction=====================================
               
      ysf = 0.5
      rj = ny/2
      nrgrd = 10
c up from center
      do 42 j = rj,rj+nrgrd
         dy_grid(j) = dy
 42   continue
      do 44 j = rj+nrgrd+1,ny
         dy_grid(j) = dy +
     x     ysf*dy*(j-(rj+nrgrd+1))/(ny-(rj+nrgrd+1)) 
 44   continue

c down from center
      do 46 j = rj-nrgrd,rj-1
         dy_grid(j) = dy
 46      continue
      do 47 j = 1,rj-nrgrd-1
         ind = rj-nrgrd-j
         dy_grid(ind) = dy + 
     x     ysf*dy*(rj-nrgrd-1-ind)/(rj-nrgrd-1)
 47   continue

      qy(1) = dy
      do 48 j=2,ny
         qy(j) = qy(j-1)+dy_grid(j)
 48   continue

      do 49 j = 1,ny-1
         dy_grid(j) = qy(j+1)-qy(j)
 49   continue
      dy_grid(ny) = dy_grid(ny-1)
c======================================================================


c      do 20 j=1,ny
c         qy(j) = j*dy
c         dy_grid(j) = dy
c 20            continue

c      do 10 i=1,nx
c         qx(i) = i*dx
c         dx_grid(i) = dx
c 10            continue

c==============stretch x direction=====================================
               
      xsf = 1.0
      ri = nx/2 + ri0
      nrgrd = 8
c up from center
      do 12 i = ri,ri+nrgrd
         dx_grid(i) = dx
 12   continue
      do 14 i = ri+nrgrd+1,nx
         dx_grid(i) = dx +
     x     0.0*xsf*dx*(i-(ri+nrgrd+1))/(nx-(ri+nrgrd+1)) 
 14   continue

c down from center
      do 16 i = ri-nrgrd,ri-1
         dx_grid(i) = dx
 16      continue
      do 17 i = 1,ri-nrgrd-1
         ind = ri-nrgrd-i
         dx_grid(ind) = dx + 
     x     xsf*dx*(ri-nrgrd-1-ind)/(ri-nrgrd-1)
 17   continue

      qx(1) = dx
      do 18 i=2,nx
         qx(i) = qx(i-1)+dx_grid(i)
 18   continue

      do 19 i = 1,nx-1
         dx_grid(i) = qx(i+1)-qx(i)
 19   continue
      dx_grid(nx) = dx_grid(nx-1)
c======================================================================

c      print*,'dx...',dx_grid


c==============stretch z direction=====================================

      zsf = 0.0  !z stretch factor
      nrgrd = 0
c up from center
      do 32 k = rk,rk+nrgrd
         dz_grid(k) = delz
 32   continue
      do 34 k = rk+nrgrd+1,nz
         dz_grid(k) = delz +
     x     zsf*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1)) 
c     x     2.0*sin((k-(rk+nrgrd+1))*0.5*pi/(nz-(rk+nrgrd+1)))**2 
c                                !dz_grid(k-1) + 0.01*delz 
 34   continue

c down from center
      do 36 k = rk-nrgrd,rk-1
         dz_grid(k) = delz
 36      continue
      do 37 k = 1,rk-nrgrd-1
         ind = rk-nrgrd-k
         dz_grid(ind) = delz + 
     x     zsf*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
c     x     2.0*sin((rk-nrgrd-1-ind)*(-0.5*pi)/(rk-nrgrd-1))**2 
c                                !dz_grid(ind+1) + 0.01*delz
 37   continue

      qz(1) = 0.0
      do 38 k=2,nz
c         write(*,*) 'dz_grid...',k,dz_grid(k)
         qz(k) = qz(k-1)+dz_grid(k)
 38   continue

      do 39 k = 1,nz-1
         dz_grid(k) = qz(k+1)-qz(k)
 39   continue
      dz_grid(nz) = dz_grid(nz-1)
c======================================================================


      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      zrat(1) = 0.5
      zrat(nz) = 0.5
      do 40 k=2,nz-1
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) -
     x                ((qz(k) + qz(k-1))/2.0)
               
         zplus = (qz(k+1) + qz(k))/2.0
         zminus = (qz(k) + qz(k-1))/2.0
         zrat(k) = (qz(k) - zminus)/(zplus - zminus)
 40   continue


      dx_cell(1) = dx_grid(1)
      dx_cell(nx) = dx_grid(nx)
      xrat(1) = 0.5
      xrat(nx) = 0.5
      do 50 i=2,nx-1
         dx_cell(i) = ((qx(i+1) + qx(i))/2.0) -
     x                ((qx(i) + qx(i-1))/2.0)
         xplus = (qx(i+1) + qx(i))/2.0
         xminus = (qx(i) + qx(i-1))/2.0
         xrat(i) = (qx(i) - xminus)/(xplus - xminus)
 50   continue

      dy_cell(1) = dy_grid(1)
      dy_cell(ny) = dy_grid(ny)
      yrat(1) = 0.5
      yrat(ny) = 0.5
      do 60 j=2,ny-1
         dy_cell(j) = ((qy(j+1) + qy(j))/2.0) -
     x                ((qy(j) + qy(j-1))/2.0)
         yplus = (qy(j+1) + qy(j))/2.0
         yminus = (qy(j) + qy(j-1))/2.0
         yrat(j) = (qy(j) - yminus)/(yplus - yminus)
 60   continue

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
      end subroutine grd8
c----------------------------------------------------------------------

      SUBROUTINE grd_no_strech()
c----------------------------------------------------------------------

      integer nrgrd

      real zplus,zminus,xplus, xminus, yplus, yminus
      real zsf,xsf,ysf

      rj=ny/2
      ri=nx/2
      rk=nz/2


c==============stretch y direction=====================================
               
      ysf = 0.0
      rj = ny/2
      nrgrd = 0
c up from center
      do 42 j = rj,rj+nrgrd
         dy_grid(j) = dy
 42   continue
      do 44 j = rj+nrgrd+1,ny
         dy_grid(j) = dy +
     x     ysf*dy*(j-(rj+nrgrd+1))/(ny-(rj+nrgrd+1)) 
 44   continue

c down from center
      do 46 j = rj-nrgrd,rj-1
         dy_grid(j) = dy
 46      continue
      do 47 j = 1,rj-nrgrd-1
         ind = rj-nrgrd-j
         dy_grid(ind) = dy + 
     x     ysf*dy*(rj-nrgrd-1-ind)/(rj-nrgrd-1)
 47   continue

      qy(1) = dy
      do 48 j=2,ny
         qy(j) = qy(j-1)+dy_grid(j)
 48   continue

      do 49 j = 1,ny-1
         dy_grid(j) = qy(j+1)-qy(j)
 49   continue
      dy_grid(ny) = dy_grid(ny-1)

c==============stretch x direction=====================================
               
      xsf = 0.0
      ri = nx/2 + ri0
      nrgrd = 0
c up from center
      do 12 i = ri,ri+nrgrd
         dx_grid(i) = dx
 12   continue
      do 14 i = ri+nrgrd+1,nx
         dx_grid(i) = dx +
     x     0.0*xsf*dx*(i-(ri+nrgrd+1))/(nx-(ri+nrgrd+1)) 
 14   continue

c down from center
      do 16 i = ri-nrgrd,ri-1
         dx_grid(i) = dx
 16      continue
      do 17 i = 1,ri-nrgrd-1
         ind = ri-nrgrd-i
         dx_grid(ind) = dx + 
     x     xsf*dx*(ri-nrgrd-1-ind)/(ri-nrgrd-1)
 17   continue

      qx(1) = dx
      do 18 i=2,nx
         qx(i) = qx(i-1)+dx_grid(i)
 18   continue

      do 19 i = 1,nx-1
         dx_grid(i) = qx(i+1)-qx(i)
 19   continue
      dx_grid(nx) = dx_grid(nx-1)


c==============stretch z direction=====================================

      zsf = 0.0  !z stretch factor
      nrgrd = 0
c up from center
      do 32 k = rk,rk+nrgrd
         dz_grid(k) = delz
 32   continue
      do 34 k = rk+nrgrd+1,nz
         dz_grid(k) = delz +
     x     zsf*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1)) 
c     x     2.0*sin((k-(rk+nrgrd+1))*0.5*pi/(nz-(rk+nrgrd+1)))**2 
c                                !dz_grid(k-1) + 0.01*delz 
 34   continue

c down from center
      do 36 k = rk-nrgrd,rk-1
         dz_grid(k) = delz
 36      continue
      do 37 k = 1,rk-nrgrd-1
         ind = rk-nrgrd-k
         dz_grid(ind) = delz + 
     x     zsf*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
c     x     2.0*sin((rk-nrgrd-1-ind)*(-0.5*pi)/(rk-nrgrd-1))**2 
c                                !dz_grid(ind+1) + 0.01*delz
 37   continue

      qz(1) = 0.0
      do 38 k=2,nz
c         write(*,*) 'dz_grid...',k,dz_grid(k)
         qz(k) = qz(k-1)+dz_grid(k)
 38   continue

      do 39 k = 1,nz-1
         dz_grid(k) = qz(k+1)-qz(k)
 39   continue
      dz_grid(nz) = dz_grid(nz-1)
c======================================================================


      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      zrat(1) = 0.5
      zrat(nz) = 0.5
      do 40 k=2,nz-1
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) -
     x                ((qz(k) + qz(k-1))/2.0)
               
         zplus = (qz(k+1) + qz(k))/2.0
         zminus = (qz(k) + qz(k-1))/2.0
         zrat(k) = (qz(k) - zminus)/(zplus - zminus)
 40   continue


      dx_cell(1) = dx_grid(1)
      dx_cell(nx) = dx_grid(nx)
      xrat(1) = 0.5
      xrat(nx) = 0.5
      do 50 i=2,nx-1
         dx_cell(i) = ((qx(i+1) + qx(i))/2.0) -
     x                ((qx(i) + qx(i-1))/2.0)
         xplus = (qx(i+1) + qx(i))/2.0
         xminus = (qx(i) + qx(i-1))/2.0
         xrat(i) = (qx(i) - xminus)/(xplus - xminus)
 50   continue

      dy_cell(1) = dy_grid(1)
      dy_cell(ny) = dy_grid(ny)
      yrat(1) = 0.5
      yrat(ny) = 0.5
      do 60 j=2,ny-1
         dy_cell(j) = ((qy(j+1) + qy(j))/2.0) -
     x                ((qy(j) + qy(j-1))/2.0)
         yplus = (qy(j+1) + qy(j))/2.0
         yminus = (qy(j) + qy(j-1))/2.0
         yrat(j) = (qy(j) - yminus)/(yplus - yminus)
 60   continue

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
      end subroutine grd_no_strech



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
