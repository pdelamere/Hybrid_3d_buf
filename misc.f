
      MODULE misc

      USE global
      USE boundary
      USE random_utils, notpi => pi

      contains


c----------------------------------------------------------------------
      real FUNCTION pad_ranf()
c This is the random number generator that works on foo.
c----------------------------------------------------------------------

      call random_number(pad_ranf)

      return
      end FUNCTION pad_ranf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE random_initialize ( )
c----------------------------------------------------------------------
      integer, allocatable :: seed(:)
      integer :: n,istat

      call random_seed(size = n)
      allocate(seed(n))

      open(unit=999,file="/dev/urandom", access="stream",
     x form="unformatted", action="read", status="old", iostat=istat)

      if(istat == 0) then
      read(999) seed
      close(999)
      else
      write(*,*) "Random_Initialize error"
      stop
      endif


      call random_seed(put=seed)

      end SUBROUTINE random_initialize
c----------------------------------------------------------------------

!----------------------------------------------------------------------
      SUBROUTINE debug_random_initialize ( )
!----------------------------------------------------------------------
      integer, allocatable :: seed(:)
      integer :: n,istat
      logical :: ex

      call random_seed(size = n)
      allocate(seed(n))

      inquire(file="./seeds/"//int_to_str(my_rank)//"_seed",
     x exist=ex)


      if(ex) then
      open(unit=999,file="./seeds/"//int_to_str(my_rank)//"_seed",
     x access="stream", form="unformatted", action="read",
     x status="old", iostat=istat)
      read(999) seed
      close(999)
      call random_seed(put=seed)

      else

      call random_initialize()
      call random_seed(get=seed)
      open(unit=998,file="./seeds/"//int_to_str(my_rank)//"_seed",
     x form="unformatted")
      write(998) seed
      close(998)

      endif


      end SUBROUTINE debug_random_initialize
!----------------------------------------------------------------------
      pure function int_to_str(num)
      integer,intent(in) :: num
      character(len=12) :: temp
      character(len=:), allocatable :: int_to_str

      if(num .eq. 0) then

      allocate(character(len=1) :: int_to_str)
      int_to_str = "0"

      else if(num .lt. 0) then

      allocate(character(len=int(log10(float(-num)))+2) :: int_to_str)
      write(temp,*) num
      temp = adjustl(temp)
      int_to_str = trim(temp)

      else

      allocate(character(len=int(log10(float(num)))+1) :: int_to_str)
      write(temp,*) num
      temp = adjustl(temp)
      int_to_str = trim(temp)

      endif
      end function int_to_str

c----------------------------------------------------------------------
      SUBROUTINE get_bndry_Eflux(b1,E)
c----------------------------------------------------------------------
     
      real b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3)

      real vol
      real uf_flux
      real exb_flux

      real mO_q
      integer i,j,k,m

      mO_q = mO/q



c Energy flux through boundary faces

cc i = 2 face (outflow boundary)

      do 20 j=2,ny
         do 20 k=2,nz
            m=1
            i=2
           exb_flux = (mO_q)**2*(1.0/mu0)*dt*dy_cell(j)*dz_cell(k)*
     x           (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux + exb_flux 
                               !+ sign since pos is flux into domain
 20         continue

cc i = nx face

      do 30 j=2,ny
         do 30 k=2,nz
            m=1
            i=nx
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dy_cell(j)*dz_cell(k)*
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
           exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx_cell(i)*dz_cell(k)*
     x            (-E(i,j,k,1)*b1(i,j,k,3)+E(i,j,k,3)*b1(i,j,k,1))*
     x            km_to_m**2
            bndry_Eflux = bndry_Eflux + exb_flux 
                               !+ sign since neg is flux into domain
 40         continue

cc j = ny face

      do 50 i=2,nx
         do 50 k=2,nz
            m=2
            j=ny
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx_cell(i)*dz_cell(k)*
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
           exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx_cell(i)*dy_cell(j)*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux + exb_flux 
                               !+ sign since pos is flux into domain
 60         continue

c k = nz face

      do 70 i=2,nx
         do 70 j=2,ny
            k=nz-1
            exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx_cell(i)*dy_cell(j)*
     x             (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))*
     x             km_to_m**3
            bndry_Eflux = bndry_Eflux - exb_flux 
                               !- sign since neg is flux into domain
 70         continue

      return
      end SUBROUTINE get_bndry_Eflux
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE get_beta()
c----------------------------------------------------------------------
      ! This isn't really used anymore, the whole beta value is in the
      ! array beta_p, so we set beta=1.0


      beta = 1.0

      return
      end SUBROUTINE get_beta
cc----------------------------------------------------------------------

      SUBROUTINE get_grad_P(np, grad_P)
      ! grad_P is filled with the pressure term that is needed for the
      ! generalized ohms law in the normalized units of the simulation.
      ! It is NOT just the gradient of P. It is poorly named, so keep
      ! that in mind any time you see grad_P in the hybrid code.
      real np(nx,ny,nz)
      real grad_P(nx,ny,nz,3)

      real np1, gdnp, a0

      integer :: i,j,k

      do i=2,nx-1
      do j=2,ny-1
      do k=2,nz-1
        np1  = 0.5*(np(i+1,j,k)+np(i,j,k))
        gdnp = (np(i+1,j,k)-np(i,j,k))/dx_grid(i)
        a0 = kboltz*e_temperature/(mion*np1)
        grad_P(i,j,k,1) = a0*gdnp

        np1  = 0.5*(np(i,j+1,k)+np(i,j,k))
        gdnp = (np(i,j+1,k)-np(i,j,k))/dx_grid(i)
        a0 = kboltz*e_temperature/(mion*np1)
        grad_P(i,j,k,2) = a0*gdnp

        np1  = 0.5*(np(i,j,k+1)+np(i,j,k))
        gdnp = (np(i,j,k+1)-np(i,j,k))/dx_grid(i)
        a0 = kboltz*e_temperature/(mion*np1)
        grad_P(i,j,k,3) = a0*gdnp
      enddo
      enddo
      enddo
      call periodic(grad_P)
      

      end SUBROUTINE get_grad_P

      end MODULE misc
