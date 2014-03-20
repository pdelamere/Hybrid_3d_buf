
      MODULE misc

      USE global
c      USE gutsf
c      USE boundary

      contains


c----------------------------------------------------------------------
      real FUNCTION pad_ranf()
c This is the random number generator that works on foo.
c----------------------------------------------------------------------
      include 'incurv.h'

c      integer function irand
c      integer iflag, irnum
c      external irand

c      irnum = irand(iflag)

c      ranf = irnum/32767.0
c      ranf = irnum/2147483647.

      call random_number(pad_ranf)

      return
      end FUNCTION pad_ranf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE random_initialize ( seed_input )
c----------------------------------------------------------------------
!**********************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!     Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed
!     However, if the input value is 0, the routine will come up with
!     its own "suggestion", based on the system clock.
!     
c      implicit none
      
      integer count
      integer  count_max
      integer  count_rate
      logical, parameter :: debug = .false.
      integer  i
      integer  seed
      integer  seed_input
      integer, allocatable :: seed_vector(:)
      integer seed_size
      real    t
      integer, parameter :: warm_up = 100
      
      seed = seed_input
!     
!     Initialize the random seed routine.
!     
      call random_seed ( )
!     
!     Determine the size of the random number seed vector.
!     
      call random_seed ( size = seed_size )
!     
!     Allocate a vector of the right size to be used as a random seed.
!     
      allocate ( seed_vector(seed_size) )
!     
!     If the user supplied a SEED value, use that.
!     
!     Otherwise, use the system clock value to make up a value that is
!     likely to change based on when this routine is called.
!     
      if ( seed /= 0 ) then
         
         if ( debug ) then
            write (*,*) ' '
            write (*,*) 'RANDOM_INITIALIZE'
            write (*,*) 'Initialize RANDOM_NUMBER, 
     x                                 user SEED = ', seed
         end if
         
      else
         
         call system_clock ( count, count_rate, count_max )
         
      seed = count
      
      if ( debug ) then
         write (*,*) ' '
         write (*,* ) 'RANDOM_INITIALIZE'
         write (*,* ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ',
     x                 seed
      end if
      
      end if
!     
!     Set the seed vector.  We don't know the significance of the
!     individual entries of the internal seed vector, so we'll just set
!     all entries to SEED.
!     
      seed_vector(1:seed_size) = seed
!     
!     Now call RANDOM_SEED, and tell it to use this seed vector.
!     
      call random_seed ( put = seed_vector(1:seed_size) )
!     
!     Free up the seed space.
!     
      deallocate ( seed_vector )
!     
!     Call the random number routine a bunch of times just to "warm it up".
!     
      do i = 1, warm_up
         call random_number ( harvest = t )
      enddo

      return
      end SUBROUTINE random_initialize
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      real FUNCTION ranf()
c This is the random number generator that works on foo.
c----------------------------------------------------------------------
      include 'incurv.h'

c      integer function irand
c      integer iflag, irnum
c      external irand

c      irnum = irand(iflag)

c      ranf = irnum/32767.0
c      ranf = irnum/2147483647.

      call random_number(ranf)

      return
      end FUNCTION ranf
c----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE Col_Freq(nu,aj)
cc----------------------------------------------------------------------
c      include 'incurv.h'
      
c      real nu(nx,ny,nz),    !collision frequency
c     x     aj(nx,ny,nz,3)    !particle density
cc     x     up(nx,ny,nz,3),  !particle bulk velocity
cc     x     nf(nx,ny,nz),    !fluid density
cc     x     uf(nx,ny,nz,3)   !fluid velocity

cc      parameter(sigma = 1)   !need a reasonable coulomb collsion xsec
cc      real ntot             !total ion density
cc      real ui(3)
c       real maxaj

cc      maxaj = 0.0
cc      do 5 i=1,nx
cc         do 5 j=1,ny
cc            do 5 k=1,nz
cc               do 5 m=1,3
cc                  if (abs(aj(i,j,k,m)) .gt. maxaj) then
cc                     maxaj = abs(aj(i,j,k,m))
cc                     endif
cc 5             continue

cc      write(*,*) 'maxaj....',maxaj
cc      if (maxaj .lt. 1.0) maxaj = 1.0
cc      write(*,*) 'maxaj....',maxaj

c      do 10 i=1,nx
c         do 10 j=1,ny
c            do 10 k=1,nz    !Kelley p. 464
c               nu(i,j,k) = 65.0 
cc   + 100.0*abs(aj(i,j,k,1))/maxaj +
cc     x                             100.0*abs(aj(i,j,k,2))/maxaj +
cc     x                             100.0*abs(aj(i,j,k,3))/maxaj 
cc               ntot = np(i,j,k) + nf(i,j,k)
cc               do 20 m=1,3
cc                  ui(m) = (np(i,j,k)*up(i,j,k,m)/ntot) + 
cc     x                      (nf(i,j,k)*uf(i,j,k,m)/ntot))
cc 20               continue
cc                  nu(i,j,k) = sqrt(ui(1)**2 + ui(2)**2 + 
cc     x                             ui(3)**2)*sigma
c 10            continue

cc      do 20 j=1,ny
cc         do 20 k=1,nz
cc            nu(2,j,k) = 250.0
cc            nu(3,j,k) = 200.0
cc            nu(4,j,k) = 150.0
cc            nu(nx,j,k) = 250.0
cc            nu(nx-1,j,k) = 200.0
cc            nu(nx-2,j,k) = 150.0
cc 20         continue


cc      do 30 i=1,nx
cc         do 30 k=1,nz
cc            nu(i,2,k) = 250.0
cc            nu(i,3,k) = 2000.0
cc            nu(i,4,k) = 150.0
cc            nu(i,ny,k) = 250.0
cc            nu(i,ny-1,k) = 200.0
cc            nu(i,ny-2,k) = 150.0
cc 30         continue


cc      do 40 i=1,nx
cc         do 40 j=1,ny
cc            nu(i,j,2) = 250.0
cc            nu(i,j,3) = 200.0
cc            nu(i,j,4) = 150.0
cc            nu(i,j,nz) = 250.0
cc            nu(i,j,nz-1) = 200.0
cc            nu(i,j,nz-2) = 150.0
cc 40         continue

c      return
c      end
cc----------------------------------------------------------------------


c----------------------------------------------------------------------
c      SUBROUTINE periodic_bc(i,j,k,ip,im,jp,jm,kp,km)
c----------------------------------------------------------------------
c      include 'incurv.h'

c      if (i .eq. 1) then
c         im = nx
c      else
c         im = i-1
c      endif

c      if (j .eq. 1) then
c         jm = ny
c      else
c         jm = j-1
c      endif

c      if (k .eq. 1) then
c         km = nz
c      else
c         km = k-1
c      endif
 
c      if (i .eq. nx) then
c         ip = 1
c      else
c         ip = i+1
c      endif
 
c      if (j .eq. ny) then
c         jp = 1
c      else
c         jp = j+1
c      endif

c      if (k .eq. nz) then
c         kp = 1
c      else
c         kp = k+1
c      endif

         
c      return
c      end
c----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE check_np(np,nf)
cc----------------------------------------------------------------------
c      include 'incurv.h'
      
c      real np(nx,ny,nz)
c      real nf(nx,ny,nz)

c      real maxnp
c      real nfrac

c      maxnp = 0.0
      
c      do 10 i=1,nx
c         do 10 j=1,ny
c            do 10 k=1,nz
cc               nfrac = np(i,j,k)/nf(i,j,k)
cc               if (nfrac .gt. 10.0) then
cc                  np(i,j,k) = nf(i,j,k)*10.0
cc                  endif

c               if (np(i,j,k) .gt. maxnp) then
c                  maxnp = np(i,j,k)
c                  endif
c 10            continue

c      write(*,*) 'maxnp, nf....',maxnp,nf(1,1,1)
c      write(*,*) 'np/nf....',maxnp/nf(1,1,1)
c      write(*,*)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE fix_b1z(b1)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real b1(nx,ny,nz,3)
c      real offsetx,offsety,offsetz

c      offsetx = (b1(3,3,3,1) + b1(nx-1,3,3,1) + b1(3,ny-1,3,1) + 
c     x         b1(nx-1,ny-1,3,1) + b1(3,3,nz-1,1) + b1(nx-1,3,nz-1,1) +
c     x         b1(3,ny-1,nz-1,1) + b1(nx-1,ny-1,nz-1,1))/8.0


c      offsety = (b1(3,3,3,2) + b1(nx-1,3,3,2) + b1(3,ny-1,3,2) + 
c     x         b1(nx-1,ny-1,3,2) + b1(3,3,nz-1,2) + b1(nx-1,3,nz-1,2) +
c     x         b1(3,ny-1,nz-1,2) + b1(nx-1,ny-1,nz-1,2))/8.0


c      offsetz = (b1(3,3,3,3) + b1(nx-1,3,3,3) + b1(3,ny-1,3,3) + 
c     x         b1(nx-1,ny-1,3,3) + b1(3,3,nz-1,3) + b1(nx-1,3,nz-1,3) +
c     x         b1(3,ny-1,nz-1,3) + b1(nx-1,ny-1,nz-1,3))/8.0

c      write(*,*) 'b1x offset....',offsetx
c      write(*,*) 'b1y offset....',offsety
c      write(*,*) 'b1z offset....',offsetz
c      write(*,*)

c      do 10 i=1,nx
c         do 10 j=1,ny
c            do 10 k=1,nz
c               b1(i,j,k,1) = b1(i,j,k,1) - (0.0 + offsetx)
c               b1(i,j,k,2) = b1(i,j,k,2) - (0.0 + offsety)
c               b1(i,j,k,3) = b1(i,j,k,3) - (0.0 + offsetz)
c 10            continue

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_chex_rate(np,nn,up,m,chex_rate)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real np(nx,ny,nz)
c      real nn(nx,ny,nz)
c      real up(nx,ny,nz,3)
c      real vn(3)

c      real cx,cy,cz          !neutral cloud center
c      real rx,ry,rz  
c      real t            !run time
c      real vr           !relative velocity between ions and neutrals
c      real sigma_chex
c      parameter (sigma_chex = 1.0e-24)   !km^2  check this
c      real vol

c      call Neut_Center(m,t,cx,cy,cz)
      
c      chex_rate = 0.0
c      do 10 i=1,nx
c         do 10 j=1,ny
c            do 10 k=1,nz
c               rx = qx(i) - cx
c               ry = qy(j) - cy
c               rz = qz(k) - cz
c               vn(1) = vsat + rx/t
c               vn(2) = ry/t
c               vn(3) = rz/t
c               vr = sqrt((up(i,j,k,1) - vn(1))**2 + 
c     x                   (up(i,j,k,2) - vn(2))**2 +
c     x                   (up(i,j,k,3) - vn(3))**2)
c               vol = dx*dy*dz_grid(k)
c               chex_rate = chex_rate + 
c     x                     vr*sigma_chex*np(i,j,k)*nn(i,j,k)*vol
c 10            continue 


c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE charge_exchange(np,xp,vp,vp1,m,chex_rate,
c     x                           input_p)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real np(nx,ny,nz)
c      real xp(Ni_max,3)
c      real vp(Ni_max,3)
c      real vp1(Ni_max,3)
c      real vn(3)
c      real input_p(3)

c      real cx,cy,cz          !neutral cloud center
c      real rx,ry,rz,r  
c      real t            !run time
c      real vr           !relative velocity between ions and neutrals
c      real sigma_chex
c      parameter (sigma_chex = 1.0e-24)   !km^2  check this
c      real vol
c      real dNcx
c      integer*4 nchex
c      real rnd
c      real nn           !neutral density
c      real nconst
c      real initial_E

c      nconst = vth*sqrt(pi)

c      call Neut_Center(m,t,cx,cy,cz)
      
c      initial_E = input_E
c      chex_rate = 0.0
c      nchex = 0      
c      do 10 l=1,Ni_tot

c         i=ijkp(l,1)
c         j=ijkp(l,2)
c         k=ijkp(l,3)

c         rx = qx(i) - cx
c         ry = qy(j) - cy
c         rz = qz(k) - cz
c         r = sqrt(rx**2 + ry**2 + rz**2)
c         vn(1) = vsat + rx/t
c         vn(2) = ry/t
c         vn(3) = rz/t
c         vr = sqrt((vp(l,1) - vn(1))**2 + 
c     x             (vp(l,2) - vn(2))**2 +
c     x             (vp(l,3) - vn(3))**2)

c         if (r .gt. 2.33*t) then    !2.33 km/s as distbn limit
                                    !otherwise float underflow
c            nn = 0.0
c         else
c            nn = (No/(4*pi*r*r*t*nconst)) *
c     x                exp(-(r-vo*t)**2 / (vth*t)**2)
c         endif

c            dNcx = dt*vr*sigma_chex*nn
cc            write(*,*) 'dNcx....',dNcx,dNcx/beta
c            rnd = pad_ranf()
c            if (rnd .lt. dNcx) then
c               nchex = nchex + 1
c               vol = dx*dy*dz_grid(k)
c               do 30 m=1,3             !remove neutral energy
c                  vp1(l,m) = vp(l,m)
c                  input_E = input_E -  
c     x                      0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
c                  input_p(m) = input_p(m) - mBa*vp(l,m) / beta
c 30            continue
c               vp(l,1) = vn(1)
c               vp(l,2) = vn(2)
c               vp(l,3) = vn(3)
cc               xp(l,1) = xp(l,1)
cc               xp(l,2) = xp(l,2)
cc               xp(l,3) = xp(l,3)
c               do 40 m=1,3             !add ion energy
c                  vp1(l,m) = vp(l,m)
c                  input_E = input_E +  
c     x                      0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
c                  input_p(m) = input_p(m) + mBa*vp(l,m) / beta
c 40               continue
c               endif
         
c 10      continue 


cc      write(*,*) 'nchex,chex_rate...',real(nchex)/beta,
cc     x            chex_rate/beta

c      input_chex = input_chex + (input_E - initial_E)
c      chex_rate = (real(nchex))/(dt*beta)
c      write(*,*) 'Normalized charge exchange energy gain...',
c     x            input_chex/(input_E - input_chex - input_bill),
c     x            input_E,input_chex,input_bill 
c      write(*,*) 'Charge exchange rate...',chex_rate
c

c      return
c      end SUBROUTINE charge_exchange
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_gradP(gradP,np,nf,itstep,etemp)
cc----------------------------------------------------------------------
c      include 'incurv.h'
      
c      real gradP(nx,ny,nz,3),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     etemp(nx,ny,nz)

c      real etemp0
c      parameter (etemp0 = 1.0e5)
cc      parameter (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
c      real a0,ntot

cc      write(*,*) 'tstep...',itstep,dt
cc      etemp = etemp0*exp(-itstep*dt/0.5)

c      do 5 i=1,nx
c         do 5 j=1,ny
c            do 5 k=1,nz
c               etemp(i,j,k) = etemp0*np(i,j,k)/(nf(i,j,k)+np(i,j,k))
c 5             continue

c      do 10 i=2,nx-1
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1
c               ntot = np(i,j,k) + nf(i,j,k)
c               a0 = kboltz*etemp(i,j,k)/(mO*ntot)
c               a1 = kboltz*np(i,j,k)/(mO*ntot)               
c               gradP(i,j,k,1) = a0*( 0.5*(np(i+1,j,k)+np(i,j,k))
c     x                          - 0.5*(np(i,j,k)+np(i-1,j,k)) )/dx 
c     x                  + a1*( 0.5*(etemp(i+1,j,k)+etemp(i,j,k))
c     x                    - 0.5*(etemp(i,j,k)+etemp(i-1,j,k)) )/dx
c               gradP(i,j,k,2) = a0*( 0.5*(np(i,j+1,k)+np(i,j,k))
c     x                          - 0.5*(np(i,j,k)+np(i,j-1,k)) )/dy
c     x                  + a1*( 0.5*(etemp(i,j+1,k)+etemp(i,j,k))
c     x                    - 0.5*(etemp(i,j,k)+etemp(i,j-1,k)) )/dy
c               gradP(i,j,k,3) = a0*( 0.5*(np(i,j,k+1)+np(i,j,k))
c     x                      - 0.5*(np(i,j,k)+np(i,j,k-1)) )/dz_grid(k)
c     x                  + a1*( 0.5*(etemp(i,j,k+1)+etemp(i,j,k))
c     x                - 0.5*(etemp(i,j,k)+etemp(i,j,k-1)) )/dz_grid(k)
c 10         continue

c      return
c      end
cc----------------------------------------------------------------------



cc----------------------------------------------------------------------
c      SUBROUTINE get_np_at_sat(np,m,satnp)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real np(nx,ny,nz)
c      real t
c      real cx,cy,cz
c      real Rsx,Rsy,Rsz                !offset coords of satellite
c      parameter (Rsx = -1.83)         !G1 release (km)
c      parameter (Rsy = -3.42)
c      parameter (Rsz = 1.14)
      

c      call Neut_Center(m,t,cx,cy,cz)

c      ii=nint((cx+Rsx)/dx)
c      jj=rj + nint(Rsy/dy)
c      kk=rk + nint(Rsz/dz_grid(rk)) 
      
c      write(*,*) 'sat coords........',ii,jj,kk
c      write(*,*) 'neutral center....',cx,cy,cz

c      satnp = np(ii,jj,kk)

c      return
c      end
cc----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE get_bndry_Eflux(b1,E)
c----------------------------------------------------------------------
      include 'incurv.h'
     
      real b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3)

      real vol
      real uf_flux
      real exb_flux

      real mO_q
      parameter(mO_q = mO/q)


c Energy flux through boundary faces

cc i = 2 face 

      do 20 j=2,ny
         do 20 k=2,nz
            m=1
            i=2
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*(-1)*dy*dz_cell(k)*
     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !+ sign since pos is flux into domain
 20         continue

cc i = nx face

      do 30 j=2,ny
         do 30 k=2,nz
            m=1
            i=nx
c            vol = uf(i,j,k,m)*dtsub*dy*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dy*dz_cell(k)*
     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
     x           km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 30         continue
cc**********************
cc j = 2 face

      do 40 i=2,nx
         do 40 k=2,nz
            m=2
            j=2
c            vol = uf(i,j,k,m)*dtsub*dx*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*(-1)*dx*dz_cell(k)*
     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !+ sign since neg is flux into domain
 40         continue

cc j = ny face

      do 50 i=2,nx
         do 50 k=2,nz
            m=2
            j=ny
c            vol = uf(i,j,k,m)*dtsub*dx*dz_cell(k)
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dz_cell(k)*
     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 50         continue
cc****************
c k = 2 face

      do 60 i=2,nx
         do 60 j=2,ny
            m=3
c            k=rk-20
            k=2
c            vol = uf(i,j,k,m)*dtsub*dx*dy
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*(-1)*dx*dy*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 60         continue

c k = nz face

      do 70 i=2,nx
         do 70 j=2,ny
c            m=3
c            k=rk+20
            k=nz-1
c            vol = uf(i,j,k,m)*dtsub*dx*dy
c            uf_flux = 0.5*nf(i,j,k)*vol*mO*(uf(i,j,k,m)*km_to_m)**2
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 70         continue

      return
      end SUBROUTINE get_bndry_Eflux
c----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_ui(uf,nf,up,np,ui)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real uf(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
c     x     up(nx,ny,nz,3),
c     x     np(nx,ny,nz),
c     x     ui(nx,ny,nz,3)     

c      real ntot(3),fnp(3),fnf(3)

c      do 10 i=1,nx
c         do 10 j=1,ny
c            do 10 k=1,nz

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
c     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
c     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
c     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

c               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
c               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
c               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)

c               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
c               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
c               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)

c               do 10 m=1,3
c                  ui(i,j,k,m) = fnp(m)*up(i,j,k,m)+fnf(m)*uf(i,j,k,m)
c 10               continue

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_vp_out(vprest,vp_out)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real vprest(3),
c     x     vp_out(3)

c      real theta,phi,mtheta,mphi
c      real v1(3),v2(3)
c      real vdx,vdy,vdz,vdx1,vdy1,vdz1,vm
c      real Ein,delta_E
c      real rnd
c      real anglei

c      phi = 0.0
     
c      if ((vprest(1) .gt. 0) .and. (vprest(2) .gt. 0)) then 
c         phi = atan(vprest(2)/vprest(1))
cc         write(*,*) 'First quadrant....'
c      endif
      
c      if ((vprest(1) .gt. 0) .and. (vprest(2) .le. 0)) then 
c         phi = atan(vprest(2)/vprest(1))
cc         write(*,*) 'Forth quadrant....'
c      endif

c      if ((vprest(1) .le. 0) .and. (vprest(2) .gt. 0)) then  
c         phi = pi - atan(abs(vprest(2))/abs(vprest(1)))
cc         write(*,*) 'Second quadrant....'
c      endif

c      if ((vprest(1) .le. 0) .and. (vprest(2) .le. 0)) then 
c         phi = pi + atan(abs(vprest(2))/abs(vprest(1)))
cc         write(*,*) 'Third quadrant....'
c      endif

c      Ein = vprest(1)**2 + vprest(2)**2 + vprest(3)**2

c      theta = acos(vprest(3)/sqrt(Ein))

cc      write(*,*) 'theta,phi...',theta*180.0/pi, phi*180.0/pi

c      mtheta = -theta
c      mphi = phi

cc      v1(1) = cos(mphi)*vprest(1) + sin(mphi)*vprest(2)
cc      v1(2) = -sin(mphi)*vprest(1) + cos(mphi)*vprest(2)
cc      v1(3) = vprest(3)

cc      v2(1) = cos(mtheta)*v1(1) + sin(mtheta)*v1(3)
cc      v2(2) = v1(2)
cc      v2(3) = -sin(mtheta)*v1(1) + cos(mtheta)*v1(3)

c      rnd = ranf()
c      anglei = asin(rnd)
c      delta_E = Ein*cos(anglei)*cos(anglei)

c      vm = sqrt(delta_E)

c      rnd = ranf()*2*pi
c      vdx = vm*sin(anglei)*cos(rnd)
c      vdy = vm*sin(anglei)*sin(rnd)
c      vdz = vm*cos(anglei)

cc      vdm = sqrt(vdx^2 + vdy^2 + vdz^2)

c      vdx1 = cos(theta)*vdx + sin(theta)*vdz
c      vdy1 = vdy
c      vdz1 = -sin(theta)*vdx + cos(theta)*vdz

c      vp_out(1) = cos(-phi)*vdx1 + sin(-phi)*vdy1
c      vp_out(2) = -sin(-phi)*vdx1 + cos(-phi)*vdy1
c      vp_out(3) = vdz1

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE billiard(vp,vp1,input_p,nn,timestep,bill_rate)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real vp(Ni_max,3),
c     x     vp1(Ni_max,3),
c     x     input_p(3),
c     x     nn(nx,ny,nz)


c      real cx,cy,cz          !neutral cloud center
c      real rx,ry,rz,r  
c      real t            !run time
c      real vr           !relative velocity between ions and neutrals
c      real sigma_Ba
c      real vol
c      real dNcol
c      integer*4 ncol
c      real rnd
c      real nneut           !neutral density
c      real nconst
c      real anglei
c      real vprest(3)
c      real vp_out(3)
c      real initial_E, final_E
c      real vn(3)

c      parameter (sigma_Ba = 3.0e-26)   !km^2  check this
      
c      initial_E = input_E
c      nconst = vth*sqrt(pi)

c      call Neut_Center(cx,cy,cz)
      
c      ncol = 0
c      bill_rate = 0.0      
c      do 10 l=1,Ni_tot

c         i=ijkp(l,1)
c         j=ijkp(l,2)
c         k=ijkp(l,3)

cc         rx = qx(i) - cx
cc         ry = qy(j) - cy
cc         rz = qz(k) - cz
cc         r = sqrt(rx**2 + ry**2 + rz**2)
cc         vn(1) = vsat + rx/t
cc         vn(2) = ry/t
cc         vn(3) = rz/t
c         vn(1) = 0.0
c         vn(2) = 0.0
c         vn(3) = 0.0
c         vr = sqrt((vp(l,1) - vn(1))**2 + 
c     x             (vp(l,2) - vn(2))**2 +
c     x             (vp(l,3) - vn(3))**2)

cc         write(*,*) 'r, 2.33*t....',r,2.33*t
cc         if (r .gt. 2.33*t) then    !2.33 km/s as distbn limit
cc                                    !otherwise float underflow
cc            nneut = 0.0
cc         else
cc            nneut = (5.0*No/(4*pi*r*r*t*nconst)) *
cc     x                exp(-(r-vo*t)**2 / (vth*t)**2)
cc         endif
c         nneut = nn(i,j,k)

c         dNcol = dt*vr*sigma_Ba*nneut
cc     write(*,*) 'dNcol....',dNcol,nneut
c         rnd = ranf()
c         if (rnd .lt. dNcol) then
c            ncol = ncol + 1
cc     vprest(1) = vn(1) - vp(l,1)
cc     vprest(2) = vn(2) - vp(l,2)
cc     vprest(3) = vn(3) - vp(l,3)
c            vprest(1) = -vp(l,1)
c            vprest(2) = -vp(l,2)
c            vprest(3) = -vp(l,3)
            
c            call get_vp_out(vprest,vp_out)
            
c            do 35 m=1,3         !remove ion energy
cc     vp1(l,m) = vp(l,m)
c               input_E = input_E -  
c     x              0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
c               input_p(m) = input_p(m) - mBa*vp(l,m) / beta
c 35         continue
            
c            write(*,*) 'vp in...',vp(l,1),vp(l,2),vp(l,3),
c     x           vp_out(1),vp_out(2),vp_out(3)
c            vp(l,1) = vp_out(1)+vp(l,1)
c            vp(l,2) = vp_out(2)+vp(l,2)
c            vp(l,3) = vp_out(3)+vp(l,3)
c            write(*,*) 'vp out...',vp(l,1),vp(l,2),vp(l,3)
            
c            vol = dx*dy*dz_grid(k)
c            do 40 m=1,3         !add ion energy
cc     vp1(l,m) = vp(l,m)
c               input_E = input_E +  
c     x              0.5*mBa*(vp(l,m)*km_to_m)**2 /beta
c               input_p(m) = input_p(m) + mBa*vp(l,m) / beta
c 40         continue
c         endif
         
c 10   continue 
      
c      final_E = input_E - initial_E
cc     input_bill = input_bill + (input_E - initial_E)
c      bill_rate = (real(ncol))/(dt*beta)
c      write(*,*) 'Normalized collisional energy gain...',
c     x     input_bill/(input_E - input_bill - input_chex) 
c      write(*,*) 'Collision rate...',bill_rate
      
c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_etar(np,aj)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real np(nx,ny,nz)
c      real aj(nx,ny,nz,3)

c      real ajperp

c      do 10 i = 1,nx 
c         do 10 j = 1,ny
c            do 10 k = 1,nz
c               do 10 m = 1,3
c                  etar(i,j,k,m) = 0.0
c 10            continue

c      do 20 i = 1,nx 
c         do 20 j = 1,ny
c            do 20 k = 1,nz
c               ajperp = sqrt(aj(i,j,k,1)**2 + aj(i,j,k,2)**2)
c               if ((np(i,j,k) .gt. 0.0) .and. (ajperp .gt. 1.0)) then
c                  etar(i,j,k,1) = eta_init
c                  etar(i,j,k,2) = eta_init
c                  etar(i,j,k,3) = 0
c               endif
c 20            continue

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE update_nu(nu,np,nf)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real nu(nx,ny,nz)
c      real np(nx,ny,nz)
c      real nf(nx,ny,nz)

c      real ntot
c      real sigma_coulb
c      parameter (sigma_coulb = 1.4e-19)   !units of km^2
c      real vthermal
c      parameter (vthermal = 210.0)  !km/s 1000 K
c      real numax,ntotmax

c      numax = 0.0
cc      ntotmax = 0.0
c      do 10 i=1,nx
c         do 10 j=1,ny
c            do 10 k=1,nz
c               ntot = np(i,j,k)+nf(i,j,k)
c               nu(i,j,k) = (melec/mO)*ntot*sigma_coulb*vthermal
c               if (nu(i,j,k) .gt. numax) then 
c                  numax = nu(i,j,k)
c                  endif
cc               if (ntotmax .gt. ntot) then
cc                  ntotmax = ntot
cc                  endif
c 10            continue

c      write(*,*) 'Nu max.........',numax
cc      write(*,*) 'np+nf max.....',ntotmax

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_nuin(nuin,nn,uf)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real nuin(nx,ny,nz)
c      real nn(nx,ny,nz)
c      real uf(nx,ny,nz,3)

c      real sigma_in
c      parameter (sigma_in = 3.0e-26)

c      do i = 1,nx
c         do j = 1,ny
c            do k = 1,nz
c               nuin(i,j,k) = sqrt(uf(i,j,k,1)**2 + uf(i,j,k,2)**2 + 
c     x                       uf(i,j,k,3)**2)*sigma_in*nn(i,j,k)
c            enddo
c         enddo
c      enddo

c      return
c      end
cc----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_beta()
c----------------------------------------------------------------------
      include 'incurv.h'

c      real Np_tot
      real src(200)
      real r_xyz(200)
      real tau
      real Nofr(200),Np_tot

      tau = tau_photo
      
      src(2) = Qo !1.5e27   !mol/s at Pluto

      do i = 1,200 
         r_xyz(i) = i*dx 
      enddo

c      Np_tot = 0.0
c      do i = 1,S_radius
cc         Nofr(i) =  (N_o - src(i)*tau)*exp(-r_xyz(i)/(tau*vrad)) + 
cc     x                src(i)*tau         
c         Nofr(i) = Qo*dx/vrad
c         Np_tot = Np_tot + Nofr(i)*dt/tau
cc         write(*,*) 'S...',src(i)
cc         write(*,*) Nofr(i),Np_tot
c      enddo
cc      write(*,*) sum(Nofr),Np_tot

c divide particles up between procnum processors      
c      beta = (Ni_tot_sys/((nx*dx*ny*dy*nz*delz))/nf_init
      beta = (Ni_tot_sys/((qx(nx)-qx(1))*(qy(ny-1)-qy(1))*
     x                    (qz(nz-1)-qz(1))))/nf_init


c      beta = (5e6/(nx*dx*ny*dy*nz*delz))/nf_init

c      write(*,*) 'beta...',beta,Ni_tot,nx*dx*ny*dy*nz*delz,nf_init

c      dNi = (beta*Np_tot)/procnum
c      dNi = (beta*Np_tot)   ! use this for all pickup on one processor
      dNi_sw = nf_init*vsw*((ny-2)*dy*(nz-2)*delz)*(dt/2)*beta
      print *,'dNi_sw....',dNi_sw

c      write(*,*) 'dNi...',dNi
c      if (dNi .lt. 1.0) then 
c         write(*,*) 'dNi lt 1.0....'
c         stop
c      endif

c      write(*,*) 'dNi...',dNi,Np_tot
c      dNi = Ni_max
c      beta = (dNi/Np_tot)

      write(*,*) 'beta, dNi....',beta,dNi

      return
      end SUBROUTINE get_beta
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_dipole(bdp)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real bdp(nx,ny,nz,3)
c      real Bo,phi,lam,br,blam,qdpx,qdpy,qdpz
c      real x,y,z,r,M
c      real eoverm
c      parameter(eoverm = q/mO)

c      Bo = 300e-9
c      M = Bo*(5.0)**3
      
c      x0 = dx/2
c      y0 = dy/2
c      z0 = dz_grid(nz/2)/2
      
c      qdpx = qx(ri) + x0
c      qdpy = qy(rj) + y0
c      qdpz = qz(rk) + z0
      
c      do 10 i = 1,nx 
c         do 10 j = 1,ny
c            do 10 k = 1,nz
c               x = qx(i) - qdpx
c               y = qy(j) - qdpy
c               z = qz(k) - qdpz
c               r = sqrt(x**2 + y**2 + z**2)
c               phi = atan(y/x)
c               if ((x .ge. 0) .and. (y .ge. 0)) phi = atan(y/x)
c               if ((x .le. 0) .and. (y .ge. 0)) then 
c                  phi = pi/2 + atan(abs(x/y))
c               endif
c               if ((x .le. 0) .and. (y .le. 0)) then 
c                  phi = pi + atan(abs(y/x))
c               endif
c               if ((x .ge. 0) .and. (y .le. 0)) then 
c                  phi = (3.0*pi/2.0) + atan(abs(x/y))
c               endif

cc               write(*,*) 'phi...',phi,qx(i),qdpx
               
c               lam = asin(z/r)
               
c               br = -M*sin(lam)/r**3
c               blam = M*0.5*cos(lam)/r**3
c               bdp(i,j,k,3) = br*sin(lam) + blam*cos(lam)
c               bx = br*cos(lam) + blam*sin(lam)
c               bdp(i,j,k,1) = bx*cos(phi)
c               bdp(i,j,k,2) = bx*sin(phi)
c               bdp(i,j,k,3) = bdp(i,j,k,3)
c 10         continue
c            bdp(:,:,:,1) = bdp(:,:,:,1)*eoverm
c            bdp(:,:,:,2) = bdp(:,:,:,2)*eoverm
c            bdp(:,:,:,3) = bdp(:,:,:,3)*eoverm

c            open(10,file='bdp.dat',status='unknown',
c     x           form='unformatted')
c            write(10) bdp
c            close(10)
c            return
c            end
cc----------------------------------------------------------------------


      end MODULE misc
