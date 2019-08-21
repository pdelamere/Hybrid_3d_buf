        Module random_utils
        use arms_interface
        use, intrinsic :: iso_c_binding, only: c_ptr, c_double
        use, intrinsic :: iso_fortran_env, only: error_unit
        implicit none
        real :: pi
        parameter( pi = 3.14159 )
        
        contains
        real function ranf()
        ! function version of random_number. Samples a random number
        ! from the range 0 <= x < 1
        call random_number(ranf)
        end function ranf

        subroutine sphere(r, x, y, z)
        ! Select a random point from the surface of a sphere with radius
        ! r. Returns its coordinates in x, y, z.
        ! This is the method described in Marsaglia (1972)
        ! in Ann Math Stat Vol. 42
        real, intent(in) :: r
        real, intent(out) :: x, y, z
        real :: v1, v2, S
        real :: rt
        do
          v1 = 2*ranf() - 1
          v2 = 2*ranf() - 1
          S = v1**2 + v2**2
          if (S .lt. 1.0) exit
        end do

        rt = sqrt(1-S)

        x = r*2.0*v1*rt
        y = r*2.0*v2*rt
        z = r*(1.0-2.0*S)
        end subroutine

        real function unnorm_VS_pdf(w, b)
        ! This is the unnormalized Vasyliunas as Siscoe 'shell'
        ! distribution function
        real, intent(in) :: w, b
        real :: x
        if ((w .lt. 0.0) .or. (w .gt. 1.0)) then
            unnorm_VS_pdf = 0.0
        else
            x = w**(-3.0/2.0)
            unnorm_VS_pdf = x*exp(-4.407*x)
        endif
        end function

        real function M_VS_dist(b)
        ! Calculate the maximum of the unnorm_VS_pdf. Used for the
        ! rejection sampling algorithm, VS_dist.
        real :: b
        if(b .lt. 1.0) then
            M_VS_dist = exp(-1.0)/b
        else
            M_VS_dist = exp(-b)
        endif
        end function

        real function VS_dist(b)
        ! Sample a value from the Vasyliunas and Siscoe 'shell'
        ! distribution. The input shape parameter is 
        ! b = \frac{\lambda \theta}{r \sin(\theta)}. Notably, this does
        ! not depend on the ionization rate \beta. This is because \beta
        ! does not affect the shape of the distribution, only the final
        ! normalization (i.e. the density).
        real, intent(in) :: b
        real :: y, u
        do
        y = ranf()
        u = ranf()
        if (u .lt. unnorm_VS_pdf(y,b)/M_VS_dist(b)) then
            VS_dist = y
            return
        endif
        enddo
        end function

        subroutine shell(u_inj, b, vx, vy, vz)
        ! A sampling algorithm for the V&S shell distribution. It uses
        ! sphere sampling to sample isotropically, and sets a radius by
        ! sampling from the V&S 'shell' distribution.
        ! A complication is that the sampling algorithm is adaptive and
        ! requires us to hang onto a pointer unique for each b. There's
        ! some extra logic to handle that. It's set up so that
        ! arms_inst(i) is the pointer that cooresponds with the arms
        ! instance for bs(i).
        real, intent(in)  :: u_inj, b
        real, intent(out) :: vx, vy, vz
        real :: samp(1) ! Only make one sample at a time

        integer, parameter :: max_instances = 5
        integer, save :: N_instances = 0
        real, save :: bs(max_instances)
        type(c_ptr), save :: arms_inst(max_instances)
        integer i,ind

        ! Check if b is in bs
        ind = -1
        do i=1, N_instances
        if (b .eq. bs(i)) then
            ind = i
            exit
        endif
        enddo
        ! Setup an arms instance if b isn't found in bs
        if (ind .eq. -1) then
            N_instances = N_instances + 1
            ind = N_instances
            if(N_instances .gt. max_instances) then
            write(error_unit,*) "More than the maximum arms instances"
            stop
            endif
            call arms_vs_setup(3, real(b, c_double), arms_inst(ind))
            bs(ind) = b
        endif

        ! Actual random sampling happens here
        call arms_sample_sp(samp, 1, arms_inst(ind))
        call sphere(u_inj*samp(1), vx, vy, vz)
        end subroutine
        end Module
