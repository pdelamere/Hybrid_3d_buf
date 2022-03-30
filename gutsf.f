
      MODULE gutsf

      USE global
      USE boundary
      USE grid_interp
      USE iso_fortran_env, only: error_unit
      
      contains


c----------------------------------------------------------------------
      SUBROUTINE f_update_tlev(b1,b12,b1p2,bt,b0)
c----------------------------------------------------------------------
      
      real b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b0(nx,ny,nz,3)

      b12 = b1
      b1 = b1p2
      bt = b0 + b1
      return
      end SUBROUTINE f_update_tlev
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE crossf2(aa, btc, aus, bus, cc)
c The cross product is formed at the main cell center.  aa and btc must
c be given already extrapolated to the main cell center.
c----------------------------------------------------------------------
      !include 'incurv.h'

      real aa(nx,ny,nz,3)        !main cell contravarient vector 
      real btc(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)
      real aus(ny,nz,3)
      real bus(ny,nz,3)

      real cus(ny,nz,3)

      real ax,ay,az,bx,by,bz    !dummy vars

      integer i,j,k

      call boundaries(aa, aus)
      call boundaries(btc, bus)

      ct(:,:,:,1) = aa(:,:,:,2)*btc(:,:,:,3) - aa(:,:,:,3)*btc(:,:,:,2)
      ct(:,:,:,2) = aa(:,:,:,3)*btc(:,:,:,1) - aa(:,:,:,1)*btc(:,:,:,3)
      ct(:,:,:,3) = aa(:,:,:,1)*btc(:,:,:,2) - aa(:,:,:,2)*btc(:,:,:,1)


c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
      
      do 60 k=2,nz-1
         do 60 j=2,ny-1
            do 60 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

      cus = ct(nx,:,:,:)
      call boundaries(cc, cus)

      return
      end SUBROUTINE crossf2
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE curlB(b1,np,aj)
c Calculates curl B / n*alpha.  The resulting "current" is called aj
c which is used in several other places in the code.  This curl is 
c performed on the main cell where B is covarient.  The resulting
c current is main cell contravarient.  Note that dz_cell is used for
c the cell dimensions since dz_grid is not equal to dz_cell on non-
c uniform grid.
c----------------------------------------------------------------------
CVD$R VECTOR
      !include 'incurv.h'

      real b1(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)
      real us(ny,nz,3)

      real curl_B(nx,ny,nz,3)      !dummy for holding curl vector
      real ntot(nx,ny,nz,3)        !total density, np + nf
      integer i,j,k
      real minden
      minden = nf_init/30.0

c      call periodic_scalar(np)
c      call periodic_scalar(nf)
      us = 0.0
      call boundaries(b1,us)
cc     call fix_normal_b(b1)

      do 10 k=2,nz-1   
         do 10 j=2,ny-1
            us(j,k,1) = 0.5*(np(nx,j,k)+nf_init)
            us(j,k,2) = 0.5*(np(nx,j,k)+np(nx,jp,k))
            us(j,k,3) = 0.5*(np(nx,j,k)+np(nx,j,kp))
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1
               ntot(i,j,k,1) = 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(i,j,k,2) = 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(i,j,k,3) = 0.5*(np(i,j,k)+np(i,j,kp))
 10    continue
      call boundaries(ntot, us)

      do 20 k=2,nz-1   
         do 20 j=2,ny-1
            us(j,k,1) = (b1(nx,j,k,3) - 
     x              b1(nx,j-1,k,3))/dy_cell(j) +
     x              (b1(nx,j,k-1,2) - b1(nx,j,k,2))/dz_cell(k)
            us(j,k,2) = (b1(nx,j,k,1) - 
     x              b1(nx,j,k-1,1))/dz_cell(k) +
     x              (b1(nx-1,j,k,3) - b1(nx,j,k,3))/dx_cell(nx)
            us(j,k,3) = (b1(nx,j,k,2) - 
     x              b1(nx-1,j,k,2))/dx_cell(nx) + 
     x              (b1(nx,j-1,k,1) - b1(nx,j,k,1))/dy_cell(j)
            do 20 i=2,nx-1
               curl_B(i,j,k,1) = (b1(i,j,k,3) - 
     x              b1(i,j-1,k,3))/dy_cell(j) +
     x              (b1(i,j,k-1,2) - b1(i,j,k,2))/dz_cell(k)
               curl_B(i,j,k,2) = (b1(i,j,k,1) - 
     x              b1(i,j,k-1,1))/dz_cell(k) +
     x              (b1(i-1,j,k,3) - b1(i,j,k,3))/dx_cell(i)
               curl_B(i,j,k,3) = (b1(i,j,k,2) - 
     x              b1(i-1,j,k,2))/dx_cell(i) + 
     x              (b1(i,j-1,k,1) - b1(i,j,k,1))/dy_cell(j)
 20            continue
      call boundaries(curl_B, us)

      where(ntot .gt. minden)
      aj = curl_B/(ntot*alpha)
      elsewhere
      aj = curl_B/(minden*alpha)
      endwhere
      

      return
      end SUBROUTINE curlB
c----------------------------------------------------------------------




c----------------------------------------------------------------------
      SUBROUTINE curlE2(E,curl_E)
c E is dual cell covarient, and curl_E will be returned as main
c cell covarient...as all magnetic fields are.  All i,j,k exclude
c boundaries.  Boundaries are taken care of in main fluid code.
c----------------------------------------------------------------------
c      include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges
      integer i,j,k

c      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               lx = qx(i+1) - qx(i)
               ly = qy(j+1) - qy(j)
               lz = qz(k+1) - qz(k)

               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/ly) - (E(i,j,k,3)/ly)
     x                       + (E(i,j,k,2)/lz) - (E(i,j,k+1,2)/lz)
               curl_E(i,j,k,2) =  (E(i,j,k,3)/lx) - (E(i+1,j,k,3)/lx)
     x                       + (E(i,j,k+1,1)/lz) - (E(i,j,k,1)/lz)
               curl_E(i,j,k,3) =  (E(i,j,k,1)/ly) - (E(i,j+1,k,1)/ly)
     x                       + (E(i+1,j,k,2)/lx) - (E(i,j,k,2)/lx)

 10          continue

c      call periodic(curl_E)

      return
      end SUBROUTINE curlE2
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE curlE(E,curl_E)
c E is dual cell covarient, and curl_E will be returned as main
c cell covarient...as all magnetic fields are.  All i,j,k exclude
c boundaries.  Boundaries are taken care of in main fluid code.
c----------------------------------------------------------------------
      !include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges
      integer i,j,k

c      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               curl_E(i,j,k,1) =  (E(i,j+1,k,3)- 
     x              E(i,j,k,3))/dy_grid(j)
     x              + (E(i,j,k,2)- 
     x              E(i,j,k+1,2))/dz_grid(k)
               curl_E(i,j,k,2) =  (E(i,j,k,3) - 
     x              E(i+1,j,k,3))/dx_grid(i)
     x              + (E(i,j,k+1,1) - 
     x              E(i,j,k,1))/dz_grid(k)
               curl_E(i,j,k,3) =  (E(i,j,k,1) - 
     x              E(i,j+1,k,1))/dy_grid(j)
     x              + (E(i+1,j,k,2) - 
     x              E(i,j,k,2))/dx_grid(i)



 10          continue

c      call periodic(curl_E)

      return
      end SUBROUTINE curlE
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE get_E(E,b0,bt,aj,up,np,nu)
      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nu(nx,ny,nz)

      real btc(nx,ny,nz,3)
      real aus(ny,nz,3)
      real bus(ny,nz,3)
      integer i,j,k,m



      real aa(nx,ny,nz,3)
      real us(ny,nz,3)

      us = 0.0 ! only correct when b upstream is curl-free
      call face_to_center(aj,aa,us)

      a = aa - up

      aus = a(nx,:,:,:)
      bus = b0(nx,:,:,:)
      call edge_to_center(bt,btc, bus)
      call crossf2(a,btc,aus,bus,c)

      !E = c + spread(nu, 4, 3)*aj
      ! c, and aj both handle their own boundaries, but since nu gets
      ! averaged we have to handle the high end boundaries separately.
      ! Hence, the loops below
      do 30 k=1,nz-1
      do 30 j=1,ny-1
      do 30 i=1,nx-1
        E(i,j,k,1) = c(i,j,k,1)
     x                + 0.5*(nu(i,j,k)+nu(i+1,j,k))*aj(i,j,k,1)
        E(i,j,k,2) = c(i,j,k,2)
     x                + 0.5*(nu(i,j,k)+nu(i,j+1,k))*aj(i,j,k,2)
        E(i,j,k,3) = c(i,j,k,3)
     x                + 0.5*(nu(i,j,k)+nu(i,j,k+1))*aj(i,j,k,3)
 30   continue

      ! We need to handle the high end boundaries. We just take the
      ! nearest nu value rather than trying to do some kind of
      ! averaging. It should be good enough, probably.
      do 31 k=1,nz
      do 31 j=1,ny
      do 31 m=1,3
        E(nx,j,k,m) = c(nx,j,k,m) + nu(nx,j,k)*aj(nx,j,k,m)
 31   continue
      do 32 k=1,nz
      do 32 i=1,nx
      do 32 m=1,3
        E(i,ny,k,m) = c(i,ny,k,m) + nu(i,ny,k)*aj(i,ny,k,m)
 32   continue
      do 33 j=1,ny
      do 33 i=1,nx
      do 33 m=1,3
        E(i,j,nz,m) = c(i,j,nz,m) + nu(i,j,nz)*aj(i,j,nz,m)
 33   continue



      return
      end SUBROUTINE get_E
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu)
c Predictor step in magnetic field update.
c----------------------------------------------------------------------
CVD$R VECTOR
c      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)   !curl of E
      real bus(ny,nz,3)
      integer i,j,k,m

c      call cov_to_contra(bt,btmf) 
c      call edge_to_center(bt,btc)
c      call get_E(E,b0,bt,btmf,aj,up,np,nu)  !E at time level m 
      call get_E(E,b0,bt,aj,up,np,nu)  !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
c      call periodic(E)
c      call fix_tangential_E(E)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3

c                  b1p2(i,j,k,m)=lww1*(b12(i+1,j,k,m)+
c     x                 b12(i-1,j,k,m)+
c     x                 b12(i,j+1,k,m)+b12(i,j-1,k,m)+
c     x                 b12(i,j,k+1,m)+b12(i,j,k-1,m))+
c     x                 lww2*b12(i,j,k,m) -
c     x                 2.0*dtsub*curl_E(i,j,k,m)

                  b1p2(i,j,k,m) = b12(i,j,k,m) - 
     x                            2.0*dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      bus = 0.0
      call boundaries(b1p2,bus)
c      call obstacle_boundary_B(b0,b1p2)
c      call fix_normal_b(b1p2)

      return
      end SUBROUTINE predict_B
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)
c The main feature here is that E must be calculated at time level
c m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
c calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
c m + 1/2, so they are used as is. 
c----------------------------------------------------------------------

      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nu(nx,ny,nz)

      real b1p1(nx,ny,nz,3)   !b1 at time level m + 1/2
      real btp1(nx,ny,nz,3)   !bt at time level m + 1/2
      real btp1mf(nx,ny,nz,3) !btp1 at contravarient position
      real btc(nx,ny,nz,3) 
      real aa(nx,ny,nz,3) 
      real aus(ny,nz,3)
      real bus(ny,nz,3)
    

      real ntot(3)            !total density np + nf
      real fnp(3),fnf(3)      !fraction np and nf of n
      real npave(3)
      integer i,j,k,m

      do 5 k=1,nz
         do 5 j=1,ny
            do 5 i=1,nx
               btp1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1)) +
     x              b0(i,j,k,1)
               b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
               btp1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2)) + 
     x              b0(i,j,k,2)
               b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
               btp1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3)) + 
     x              b0(i,j,k,3)
               b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
 5             continue

      call curlB(b1p1,np,aj)
      aus = 0.0!not always correct
      call face_to_center(aj,aa,aus)
      a = aa - up

      aus = a(nx,:,:,:)!not always correct
      bus = b0(nx,:,:,:)

      call edge_to_center(btp1,btc,bus)

      call crossf2(a,btc,aus,bus,c)
       
      !E = c + spread(nu,4,3)*aj
      ! c, and aj both handle their own boundaries, but since nu gets
      ! averaged we have to handle the high end boundaries separately.
      ! Hence, the loops below
      do 40 k=1,nz-1
      do 40 j=1,ny-1
      do 40 i=1,nx-1
        E(i,j,k,1) = c(i,j,k,1)
     x                + 0.5*(nu(i,j,k)+nu(i+1,j,k))*aj(i,j,k,1)
        E(i,j,k,2) = c(i,j,k,2)
     x                + 0.5*(nu(i,j,k)+nu(i,j+1,k))*aj(i,j,k,2)
        E(i,j,k,3) = c(i,j,k,3)
     x                + 0.5*(nu(i,j,k)+nu(i,j,k+1))*aj(i,j,k,3)
 40   continue

      ! We need to handle the high end boundaries. We just take the
      ! nearest nu value rather than trying to do some kind of
      ! averaging. It should be good enough, probably.
      do 41 k=1,nz
      do 41 j=1,ny
      do 41 m=1,3
        E(nx,j,k,m) = c(nx,j,k,m) + nu(nx,j,k)*aj(nx,j,k,m)
 41   continue
      do 42 k=1,nz
      do 42 i=1,nx
      do 42 m=1,3
        E(i,ny,k,m) = c(i,ny,k,m) + nu(i,ny,k)*aj(i,ny,k,m)
 42   continue
      do 43 j=1,ny
      do 43 i=1,nx
      do 43 m=1,3
        E(i,j,nz,m) = c(i,j,nz,m) + nu(i,j,nz)*aj(i,j,nz,m)
 43   continue
      return
      end SUBROUTINE get_Ep1
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE correct_B(b0,b1,b1p2,E,aj,up,np,nu)
c Corrector step in magnetic field update.
c----------------------------------------------------------------------
CVD$R VECTOR
c      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3),
c     x     bdp(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)            !curl of E
      real bus(ny,nz,3)
      integer i,j,k,m

      call get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)  
                                                   !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
c      call periodic(E)
c      call fix_tangential_E(E)

c      write(*,*) 'E cb...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3

c                  b1p2(i,j,k,m)=lww1*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+
c     x                 b1(i,j+1,k,m)+b1(i,j-1,k,m)+
c     x                 b1(i,j,k+1,m)+b1(i,j,k-1,m))+
c     x                 lww2*b1(i,j,k,m) -
c     x                 dtsub*curl_E(i,j,k,m)

                  b1p2(i,j,k,m) = b1(i,j,k,m) - 
     x                            dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      bus = 0.0
      call boundaries(b1p2, bus)
c      call obstacle_boundary_B(b0,b1p2)
c      call fix_normal_b(b1p2)

      return
      end SUBROUTINE correct_B
c----------------------------------------------------------------------

      SUBROUTINE update_nu(nu, nu_background, aj, bt)
      real nu(nx,ny,nz)
      real nu_background(nx,ny,nz)
      real aj(nx,ny,nz,3)
      real bt(nx,ny,nz,3)
      real btc(nx,ny,nz,3)
      real ajc(nx,ny,nz,3)       !aj at cell center
      real aj_par
      real btc_mag
      real nu_target
      real us(ny,nz,3)
      integer i,j,k

      ! nu(aj) = kappa*(aj^2 - ajc^2)*S(aj-ajc) + nu_init
      ! i.e.
      ! nu is constant below ajc
      ! nu(aj <= ajc) == nu_init
      ! and quadradic above it with
      ! nu(aj_high) == nu_high

      ! Critical current where resistivity "turns on" i.e. becomes
      ! quadradic (should be 0 - 1400 for pluto?)
      ! resistivity is constant below ajc
      real aj_crit

      ! nu(aj_high) == nu_high
      real aj_high
      real nu_high

      ! kappa is a constant calculated to make all this work
      real kappa 
      
      aj_crit = 100
      aj_high = 300
      nu_high = 0.4 ! lower hybrid frequency

      kappa = (nu_high - nu_init)/(aj_high**2 - aj_crit**2)

      call edge_to_center(bt,btc,b0_us)
      us = 0.0 ! only if upstream B is curl free
      call face_to_center(aj,ajc,us)

      do 60 i=1,nx
        do 60 j=1,ny
          do 60 k=1,nz
          btc_mag = sqrt(  btc(i,j,k,1)**2
     x                   + btc(i,j,k,2)**2
     x                   + btc(i,j,k,3)**2
     x            )
          aj_par = (ajc(i,j,k,1)*btc(i,j,k,1) 
     x           + ajc(i,j,k,2)*btc(i,j,k,2) 
     x           + ajc(i,j,k,3)*btc(i,j,k,3))/btc_mag
          if (aj_par .ge. aj_crit) then
              nu_target = nu_background(i,j,k)
     x                    + kappa*(aj_par**2 - aj_crit**2)
          else
              nu_target = nu_background(i,j,k)
          endif
          nu(i,j,k) = nu(i,j,k) + (nu_target - nu(i,j,k))/10
 60   continue
      end SUBROUTINE update_nu
     
c----------------------------------------------------------------      
      SUBROUTINE check_time_step(b0,b1,bt,np,step,error_file)
c----------------------------------------------------------------      
      
      real b0(nx,ny,nz,3)
      real b1(nx,ny,nz,3)
      real bt(nx,ny,nz,3)
      real np(nx,ny,nz)
      real ak, btot, a1, a2, womega, phi, deltat
      integer error_file
      integer step
      integer i,j,k

      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               ak = 2./dx
               btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + 
     x              bt(i,j,k,3)**2)
               a1 = ak**2*Btot/(alpha*(np(i,j,k)))
               a2 = (ak*Btot)**2/(alpha*(np(i,j,k)))
               womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
               phi = womega/ak
               deltat = dx/phi
               if(deltat .lt. 2.0*dtsub) then 
                   ! rank plus 1 is proc number
                  !write(*,*) 'time stepping error...', my_rank+1
                  !write(error_file) step
                  !write(error_file) my_rank+1
                  !write(error_file) i,j,k
                  !write(error_file) qx(i),qy(j),qz(k),gz(k)
                  !write(error_file) np
                  !write(error_file) b0
                  !write(error_file) b1
                  !write(error_file) bt

                  dtsub = dtsub/2.0
                  ntf = 2*ntf
               endif

            enddo
         enddo
      enddo

      if (ntf .gt. 50) then 
          write(error_unit,*) 
     x         'Excesive subcycling detected. ntf=', ntf
      elseif (ntf .gt. 100) then
          write(error_unit,*) 
     x    'Excesive subcycling pinned at 100. desired ntf was ntf=', ntf
          ntf = 100
          dtsub = dt/ntf
      endif
      
      return
      end subroutine check_time_step
c----------------------------------------------------------------      

      


      end MODULE gutsf













