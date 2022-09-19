        Module random_utils
        implicit none
        real :: pi
        parameter( pi = 3.14159 )
        
        contains

        function circle_sample() result(r)
            real, dimension(3) :: r
            do
                r(1) = 2*ranf() - 1
                r(2) = 2*ranf() - 1
                r(3) = r(1)**2 + r(2)**2
                if((r(3) .lt. 1) .and. (r(3) .ne. 0)) exit
            end do
        end function circle_sample

        function mpm_factor(s)
            real :: mpm_factor
            real, intent(in) :: s
            mpm_factor = sqrt(-2*log(s)/s)
        end function mpm_factor

        subroutine randn_one(n)
            real, intent(out) :: n
            real :: x,y,s,f
            real, dimension(3) :: c
            c = circle_sample()
            x = c(1)
            y = c(2)
            s = c(3)
            f = mpm_factor(s)
            n = x*f
        end subroutine randn_one

        subroutine randn_pair(n1,n2)
            real, intent(out) :: n1, n2
            real :: x,y,s,f
            real, dimension(3) :: c
            c = circle_sample()
            x = c(1)
            y = c(2)
            s = c(3)
            f = mpm_factor(s)
            n1 = x*f
            n2 = y*f
        end subroutine randn_pair

        function randn(n) result(r)
            real, dimension(n) :: r
            integer, intent(in) :: n
            integer :: i
            do i=1, n-1, 2
                call randn_pair(r(i), r(i+1))
            enddo
            if(mod(n,2) .eq. 1) then
                call randn_one(r(n))
            endif
        end function randn

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
        real, intent(in)  :: u_inj, b
        real, intent(out) :: vx, vy, vz
        call sphere(u_inj*VS_dist(b), vx, vy, vz)
        end subroutine
        end Module
