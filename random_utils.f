        Module random_utils
        implicit none
        real :: pi
        parameter( pi = 3.14159 )
        
        contains
        real function ranf()
        ! function version of random_number. Samples a random number
        ! from the range 0 <= x < 1
        call random_number(ranf)
        end function ranf

        real function std_norm()
        ! Sample a random value from a standard normal distribution
        ! using the Box-Muller method
        real :: u1, u2
        u1 = 1.0-ranf()
        u2 = 1.0-ranf()
        std_norm = sqrt(-2*alog(u1))*cos(2*pi*u2)
        end function

        real function std_MT(a)
        ! Use the standard Marsaglia-Tsang algorithm to generate Gamma
        ! random variables with a>=1. Probably don't need to call this
        ! directly. Call MT(a) instead since it allows a>=0. We use a
        ! squeeze to imporve efficiency.
        real, intent(in) :: a
        real :: d, c
        real :: x, v
        real :: u

        d = a - 1.0/3.0
        c = 1/sqrt(9.0*d)
        do
            do
                x = std_norm()
                v = (1.0 + c*x)
                if (v .gt. 0.0) exit
            end do
            v = v*v*v
            u = ranf()
            if (u .lt. 1.0 - 0.0331*(x*x)*(x*x)) then
                std_MT = d*v
                return
            endif
            if (log(u) .lt. 0.5*x*x + d*(1.0 - v + log(v))) then
                std_MT = d*v
                return
            endif
        end do
        end function

        real function MT(a)
        ! The standard MT algorithm (std_MT) only alows a>=1, but a
        ! simple transformation can boost into the range 0<=a<=1.
        ! Combined, we have an algorithm for gamma distributed random
        ! variables with a>=0.
        real, intent(in) :: a

        if (a .ge. 1) then
            MT = std_MT(a)
            return
        else
            MT = std_MT(a+1)*ranf()**(1.0/a)
            return
        endif
        end function

        real function gen_gamma(a,b,c,d)
        ! The generalized gamma adds another shape parameter, c, and
        ! allows shifting and scaling via b, and d. It's implemented as
        ! a simple transformation of a gamma random variable.
        real, intent(in) :: a,b,c,d
        gen_gamma = d + b*MT(a)**(1.0/c)
        end function

        subroutine sphere(r, x, y, z)
        ! Select a random point from the surface of a sphere with radius
        ! r. Returns its coordinates in x, y, z.
        ! This is the method due to Marsaglia (1972) in Ann Math Stat Vol. 42
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
        z = r*1.0-2.0*S

        end subroutine

        subroutine shell(u_inj, r, lambda, theta, vx, vy, vz)
        ! A sampling algorithm for the V&S shell distribution. It uses
        ! sphere sampling to sample isotropically, and sets a radius by
        ! sampling from a truncated generalized gamma distribution. The
        ! trucation is accomplished by basic rejection sampling.
        real, intent(in) :: u_inj, r, lambda, theta
        real, intent(out) :: vx, vy, vz
        real :: b, w
        b = r*sin(theta)/(lambda*theta)
        do
            w = gen_gamma(1.0/3.0, b, -3.0/2.0, 0.0)
            if (w .le. 1.0) exit
        end do
        call sphere(u_inj*w, vx, vy, vz)
        end subroutine
        end Module

