    MODULE cdf_normal_mod
! ----------------------------------------------------------------------

!                              cdf_normal_mod
!                              *=*=*=*=*=*=*=

!  -  SUBROUTINE CDF_NORMAL(WHICH, CUM, CCUM, X, MEAN, SD, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_NORMAL(X, MEAN, SD,  STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_NORMAL(X, MEAN, SD, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_NORMAL(CUM, CCUM, MEAN, SD, STATUS,
!                                       CHECK_INPUT)

!                             The Distribution
!                             ================

!         The density of the normal distribution is proportional to:

!                                 (         2)
!                                 ( (X-MEAN) )
!                              exp(----------)
!                                 (       2  )
!                                 (   2 SD   )

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    arguments is to be calculated.
!  Input Range: [ 1:4 ]
!     1.  CUM and CCUM
!     2.  X
!     3.  MEAN
!     4.  SD
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the normal distribution.
!  Range: [ 10^-10:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the normal
!    distribution.
!  Range: [ 10^-10:1-10^-10 ]
!  - REAL (dpkind) :: X. The upper limit of integration of the normal
!    density. The lower limit is -oo.
!  Range: [ -10^100:10^100 ]
!  - REAL (dpkind), OPTIONAL:: MEAN. The mean of the normal distribution.
!    If omitted, the value 0 is used.
!  Range: [ -10^100:10^100 ]
!  - REAL (dpkind), OPTIONAL:: SD. The standard deviation of the normal
!    distribution. If omitted, the value 1 is used.
!  Range: [ 10^-10:10^100 ]
!  - INTEGER, INTENT(OUT), OPTIONAL :: STATUS. Return code. Possible values:
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 X outside range
!     -5 MEAN outside range
!     -6 SD outside range
!      3 CUM + CCUM is not nearly one
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and .TRUE.
!    input argument values are not checked for validity.
! ----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
      USE biomath_mathlib_mod
      USE cdf_aux_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_normal, cdf_normal, cum_normal, inv_normal
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_normal(x,mean,sd,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_normal
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: mean, sd
        REAL (dpkind) :: x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_normal(which=1,ccum=ccum_normal,x=x,mean=mean,sd=sd, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_normal

!*********************************************************************

      SUBROUTINE cdf_normal(which,cum,ccum,x,mean,sd,status,check_input)
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum, mean, sd
        REAL (dpkind) :: x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Arrays ..
        REAL (dpkind) :: params(6)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: local_ccum, local_cum, local_mean, local_sd, z
        LOGICAL :: has_status, local_check_input
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        has_status = PRESENT(status)

! status = 0 means NO error

        IF (has_status) THEN
          status = 0
        END IF

! Check presence of arguments and fix up

        CALL check_complements(cum,ccum,the_normal%name,'cum','ccum', &
          local_cum,local_ccum,set_values=which/=1,bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF


        IF (PRESENT(mean)) THEN
          local_mean = mean
        ELSE
          local_mean = zero
        END IF

        IF (PRESENT(sd)) THEN
          local_sd = sd
        ELSE
          local_sd = one
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = x
        params(4) = local_mean
        params(5) = local_sd

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

!!! Check for consistency of which and answer

        IF (local_check_input) THEN

          IF (which==3 .AND. .NOT. PRESENT(mean)) CALL which_miss(which, &
            the_normal%name,'mean')

          IF (which==4 .AND. .NOT. PRESENT(sd)) CALL which_miss(which, &
            the_normal%name,'sd')

          CALL validate_parameters(the_normal,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF

        END IF

! ++++++++++          ++++++++++          ++++++++++
! Compute the Answers
! ++++++++++          ++++++++++          ++++++++++

        SELECT CASE (which)
        CASE (1)
          z = (x-local_mean)/local_sd

          CALL local_cum_normal(z,local_cum,local_ccum)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)
          z = local_inv_normal(local_cum,local_ccum)

          x = local_sd*z + local_mean

        CASE (3)
          z = local_inv_normal(local_cum,local_ccum)

          mean = x - local_sd*z

        CASE (4)
          z = local_inv_normal(local_cum,local_ccum)

          sd = (x-local_mean)/z
        END SELECT

      END SUBROUTINE cdf_normal

!*********************************************************************

      FUNCTION cum_normal(x,mean,sd,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_normal
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: mean, sd
        REAL (dpkind) :: x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_normal(which=1,cum=cum_normal,x=x,mean=mean,sd=sd, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_normal

!*********************************************************************

      FUNCTION inv_normal(cum,ccum,mean,sd,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_normal
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum, mean, sd
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_normal(which=2,cum=cum,ccum=ccum,x=inv_normal,mean=mean, &
          sd=sd,status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_normal

!*********************************************************************

      SUBROUTINE local_cum_normal(arg,result,ccum)
!----------------------------------------------------------------------

!                              Function

!     Computes the cumulative  of    the  normal   distribution,   i.e.,
!     the integral from -infinity to x of
!          (1/sqrt(2*pi)) exp(-u*u/2) du
!     X --> Upper limit of integration.
!                                        X is DOUBLE PRECISION
!     RESULT <-- Cumulative normal distribution.
!                                        RESULT is DOUBLE PRECISION
!     CCUM <-- Compliment of Cumulative normal distribution.
!                                        CCUM is DOUBLE PRECISION
!     Renaming of function ANORM from:
!     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
!     Package of Special Function Routines and Test Drivers"
!     acm Transactions on Mathematical Software. 19, 22-32.
!     with slight modifications to return ccum and to deal with
!     machine constants.
!----------------------------------------------------------------------
! Original Comments:
!------------------------------------------------------------------
! This function evaluates the normal distribution function:
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!   The main computation evaluates near-minimax approximations
!   derived from those in "Rational Chebyshev approximations for
!   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!   This transportable program uses rational functions that
!   theoretically approximate the normal distribution function to
!   at least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!*******************************************************************
!*******************************************************************
! Explanation of machine-dependent constants.
!   MIN   = smallest machine representable number.
!   EPS   = argument below which anorm(x) may be represented by
!           0.5  and above which  x*x  will not underflow.
!           A conservative value is the largest machine number X
!           such that   1.0 + X = 1.0   to machine precision.
!*******************************************************************
!*******************************************************************
! Error returns
!  The program returns  ANORM = 0     for  ARG .LE. XLOW.
! Intrinsic functions required are:
!     ABS, AINT, EXP
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!  Latest modification: March 15, 1992
!------------------------------------------------------------------
!------------------------------------------------------------------
!  Mathematical constants
!  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
!  THRSH is the argument for which anorm = 0.75.
!------------------------------------------------------------------
!------------------------------------------------------------------
!  Coefficients for approximation in first interval
!------------------------------------------------------------------
!------------------------------------------------------------------
!  Coefficients for approximation in second interval
!------------------------------------------------------------------
!------------------------------------------------------------------
!  Coefficients for approximation in third interval
!------------------------------------------------------------------
! .. Scalar Arguments ..
        REAL (dpkind) :: arg, ccum, result
! ..
! .. Local Scalars ..
        REAL (dpkind) :: del, eps, min, temp, x, xden, xnum, xsq, y
        INTEGER :: i
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, AINT, EPSILON, EXP, TINY
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: a(5) = (/ &
          2.2352520354606839287E00_dpkind, &
          1.6102823106855587881E02_dpkind, &
          1.0676894854603709582E03_dpkind, &
          1.8154981253343561249E04_dpkind, &
          6.5682337918207449113E-2_dpkind/)
        REAL (dpkind), PARAMETER :: b(4) = (/ &
          4.7202581904688241870E01_dpkind, &
          9.7609855173777669322E02_dpkind, &
          1.0260932208618978205E04_dpkind, &
          4.5507789335026729956E04_dpkind/)
        REAL (dpkind), PARAMETER :: c(9) = (/ &
          3.9894151208813466764E-1_dpkind, &
          8.8831497943883759412E00_dpkind, &
          9.3506656132177855979E01_dpkind, &
          5.9727027639480026226E02_dpkind, &
          2.4945375852903726711E03_dpkind, &
          6.8481904505362823326E03_dpkind, &
          1.1602651437647350124E04_dpkind, &
          9.8427148383839780218E03_dpkind, &
          1.0765576773720192317E-8_dpkind/)
        REAL (dpkind), PARAMETER :: d(8) = (/ &
          2.2266688044328115691E01_dpkind, &
          2.3538790178262499861E02_dpkind, &
          1.5193775994075548050E03_dpkind, &
          6.4855582982667607550E03_dpkind, &
          1.8615571640885098091E04_dpkind, &
          3.4900952721145977266E04_dpkind, &
          3.8912003286093271411E04_dpkind, &
          1.9685429676859990727E04_dpkind/)
        REAL (dpkind), PARAMETER :: p(6) = (/ &
          2.1589853405795699E-1_dpkind, 1.274011611602473639E-1_dpkind, &
          2.2235277870649807E-2_dpkind, 1.421619193227893466E-3_dpkind, &
          2.9112874951168792E-5_dpkind, 2.307344176494017303E-2_dpkind/)
        REAL (dpkind), PARAMETER :: q(5) = (/ &
          1.28426009614491121E00_dpkind, 4.68238212480865118E-1_dpkind, &
          6.59881378689285515E-2_dpkind, 3.78239633202758244E-3_dpkind, &
          7.29751555083966205E-5_dpkind/)
        REAL (dpkind), PARAMETER :: root32 = 5.656854248E0_dpkind
        REAL (dpkind), PARAMETER :: sixten = 1.6_dpkind
        REAL (dpkind), PARAMETER :: sqrpi = &
          3.9894228040143267794E-1_dpkind
        REAL (dpkind), PARAMETER :: thrsh = 0.66291E0_dpkind
! ..
!------------------------------------------------------------------
!  Machine dependent constants
!------------------------------------------------------------------
!        eps = spmpar(1)*0.5E0_dpkind
        eps = EPSILON(one)*half
!        min = spmpar(2)
        min = TINY(one)
!------------------------------------------------------------------
        x = arg
        y = ABS(x)
        IF (y<=thrsh) THEN
!------------------------------------------------------------------
!  Evaluate  anorm  for  |X| <= 0.66291
!------------------------------------------------------------------
          xsq = zero
          IF (y>eps) xsq = x*x
          xnum = a(5)*xsq
          xden = xsq
          xnum = (xnum+a(1))*xsq
          xden = (xden+b(1))*xsq
          xnum = (xnum+a(2))*xsq
          xden = (xden+b(2))*xsq
          xnum = (xnum+a(3))*xsq
          xden = (xden+b(3))*xsq
          result = x*(xnum+a(4))/(xden+b(4))
          temp = result
          result = half + temp
          ccum = half - temp
!------------------------------------------------------------------
!  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
!------------------------------------------------------------------
        ELSE IF (y<=root32) THEN
          xnum = c(9)*y
          xden = y

          DO i = 1, 7
            xnum = (xnum+c(i))*y
            xden = (xden+d(i))*y
          END DO

          result = (xnum+c(8))/(xden+d(8))
          xsq = AINT(y*sixten)/sixten
          del = (y-xsq)*(y+xsq)
          result = EXP(-xsq*xsq*half)*EXP(-del*half)*result
          ccum = one - result

          IF (x>zero) THEN
            temp = result
            result = ccum
            ccum = temp
          END IF

!------------------------------------------------------------------
!  Evaluate  anorm  for |X| > sqrt(32)
!------------------------------------------------------------------
        ELSE
          result = zero
          xsq = one/(x*x)
          xnum = p(6)*xsq
          xden = xsq
          xnum = (xnum+p(1))*xsq
          xden = (xden+q(1))*xsq
          xnum = (xnum+p(2))*xsq
          xden = (xden+q(2))*xsq
          xnum = (xnum+p(3))*xsq
          xden = (xden+q(3))*xsq
          xnum = (xnum+p(4))*xsq
          xden = (xden+q(4))*xsq
          result = xsq*(xnum+p(5))/(xden+q(5))
          result = (sqrpi-result)/y
          xsq = AINT(x*sixten)/sixten
          del = (x-xsq)*(x+xsq)
          result = EXP(-xsq*xsq*half)*EXP(-del*half)*result
          ccum = one - result

          IF (x>zero) THEN
            temp = result
            result = ccum
            ccum = temp
          END IF

        END IF

        IF (result<min) result = zero
        IF (ccum<min) ccum = zero

      END SUBROUTINE local_cum_normal

!*********************************************************************

      FUNCTION local_inv_normal(cum,ccum)
!----------------------------------------------------------------------
! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.
! cum is the desired area to the left of local_inv_normal
! ccum is the desired area to the right
! A higher level of code checks ranges on cum and ccum and assures
! that they sum to one.
! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine
! Code modified  (but not in logic) by Barry W. Brown
! May,2001 for incorporation into dcdflib90.
! The code is again a function
!----------------------------------------------------------------------
! .. Parameters ..
        REAL (dpkind), PARAMETER :: a(8) = (/ &
          3.3871328727963666080E0_dpkind, 1.3314166789178437745E+2_dpkind &
          , 1.9715909503065514427E+3_dpkind, &
          1.3731693765509461125E+4_dpkind, &
          4.5921953931549871457E+4_dpkind, &
          6.7265770927008700853E+4_dpkind, &
          3.3430575583588128105E+4_dpkind, &
          2.5090809287301226727E+3_dpkind/)
        REAL (dpkind), PARAMETER :: b(8) = (/ 1.0_dpkind, &
          4.2313330701600911252E+1_dpkind, &
          6.8718700749205790830E+2_dpkind, &
          5.3941960214247511077E+3_dpkind, &
          2.1213794301586595867E+4_dpkind, &
          3.9307895800092710610E+4_dpkind, &
          2.8729085735721942674E+4_dpkind, &
          5.2264952788528545610E+3_dpkind/)
        REAL (dpkind), PARAMETER :: c(8) = (/ &
          1.42343711074968357734E0_dpkind, &
          4.63033784615654529590E0_dpkind, &
          5.76949722146069140550E0_dpkind, &
          3.64784832476320460504E0_dpkind, &
          1.27045825245236838258E0_dpkind, &
          2.41780725177450611770E-1_dpkind, &
          2.27238449892691845833E-2_dpkind, &
          7.74545014278341407640E-4_dpkind/)
        REAL (dpkind), PARAMETER :: d(8) = (/ 1.0_dpkind, &
          2.05319162663775882187E0_dpkind, &
          1.67638483018380384940E0_dpkind, &
          6.89767334985100004550E-1_dpkind, &
          1.48103976427480074590E-1_dpkind, &
          1.51986665636164571966E-2_dpkind, &
          5.47593808499534494600E-4_dpkind, &
          1.05075007164441684324E-9_dpkind/)
        REAL (dpkind), PARAMETER :: e(8) = (/ &
          6.65790464350110377720E0_dpkind, &
          5.46378491116411436990E0_dpkind, &
          1.78482653991729133580E0_dpkind, &
          2.96560571828504891230E-1_dpkind, &
          2.65321895265761230930E-2_dpkind, &
          1.24266094738807843860E-3_dpkind, &
          2.71155556874348757815E-5_dpkind, &
          2.01033439929228813265E-7_dpkind/)
        REAL (dpkind), PARAMETER :: f(8) = (/ 1.0_dpkind, &
          5.99832206555887937690E-1_dpkind, &
          1.36929880922735805310E-1_dpkind, &
          1.48753612908506148525E-2_dpkind, &
          7.86869131145613259100E-4_dpkind, &
          1.84631831751005468180E-5_dpkind, &
          1.42151175831644588870E-7_dpkind, &
          2.04426310338993978564E-15_dpkind/)
        REAL (dpkind), PARAMETER :: const1 = 0.180625E0_dpkind
        REAL (dpkind), PARAMETER :: const2 = 1.6E0_dpkind
        REAL (dpkind), PARAMETER :: split1 = 0.425E0_dpkind
        REAL (dpkind), PARAMETER :: split2 = 5.E0_dpkind
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: ccum, cum
! ..
! .. Local Scalars ..
        REAL (dpkind) :: p, q, r
! ..
! .. Function Return Value ..
        REAL (dpkind) :: local_inv_normal
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, LOG, MIN, SQRT
! ..
        p = MIN(cum,ccum)

        IF (p==half) THEN
          local_inv_normal = zero

          RETURN
        END IF

        q = p - half

        IF (ABS(q)<=split1) THEN
          r = const1 - q*q
          local_inv_normal = q*evaluate_polynomial(a,r)/ &
            evaluate_polynomial(b,r)
        ELSE
          r = SQRT(-LOG(p))

          IF (r<=split2) THEN
            r = r - const2
            local_inv_normal = evaluate_polynomial(c,r)/ &
              evaluate_polynomial(d,r)
          ELSE
            r = r - split2
            local_inv_normal = evaluate_polynomial(e,r)/ &
              evaluate_polynomial(f,r)
          END IF

        END IF

        local_inv_normal = ABS(local_inv_normal)

        IF (cum<half) local_inv_normal = -local_inv_normal

        RETURN

      END FUNCTION local_inv_normal

!*********************************************************************

    END MODULE cdf_normal_mod
