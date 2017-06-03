    MODULE cdf_t_mod
! ----------------------------------------------------------------------

!                                cdf_t_mod
!                                *=*=*=*=*

!  -  SUBROUTINE CDF_T( WHICH, CUM, CCUM, T, DF,   STATUS, CHECK_INPUT )
!  -  REAL (dpkind) FUNCTION CUM_T(T, DF, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_T(T, DF, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_T(CUM, CCUM, T, DF, STATUS, CHECK_INPUT)

!                             The Distribution
!                             ================

! The density is proportional to:

!                            [      2 ]
!                            [     T  ](DF+1)/2
!                            [ 1 + -- ]
!                            [     DF ]

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    arguments is to be calculated.
!  Input Range: [ 1:3 ]
!     1.  CUM and CCUM
!     2.  T
!     3.  DF
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the noncentral t
!    distribution.
!  Range: [ 10^-10:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the noncentral
!    t distribution.
!  Range: [ 10^-10:1-10^-10 ]
!  - REAL (dpkind) :: T. The upper limit of integration of the noncentral
!    t density. The lower limit is -oo.
!  Range: [ -10^100:10^100 ]
!  - REAL (dpkind) :: DF. The degrees of freedom.
!  Range: [ 10^-3:10^10 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!      0 Problem solved successfully
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 T outside range
!     -5 DF outside range
!      3 CUM + CCUM is not nearly one
!     10 can not solve cdf_beta in local_cum_t
!    -50 Answer (if any) is BELOW the LOWER search bound
!     50 Answer (if any) is ABOVE the UPPER search bound
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
      PUBLIC :: ccum_t, cdf_t, cum_t, inv_t
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_t(t,df,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_t
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: df, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_t(which=1,ccum=ccum_t,t=t,df=df,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION ccum_t

!*********************************************************************

      SUBROUTINE cdf_t(which,cum,ccum,t,df,status,check_input)
! .. Use Statements ..
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: df, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Structures ..
        TYPE (zf_locals) :: local
! ..
! .. Local Scalars ..
        REAL (dpkind) :: fx, local_ccum, local_cum, try_ccum, try_cum
        INTEGER :: zf_status
        LOGICAL :: has_status, local_check_input, match_cum
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! .. Local Arrays ..
        REAL (dpkind) :: params(6)
! ..
        has_status = PRESENT(status)

! status = 0 means NO error

        IF (has_status) THEN
          status = 0
        END IF

        CALL check_complements(cum,ccum,the_t%name,'cum','ccum', &
          local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = t
        params(4) = df

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_t,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF
        END IF

! ++++++++++          ++++++++++          ++++++++++
! Compute the Answers
! ++++++++++          ++++++++++          ++++++++++

        IF (which>1) match_cum = (local_cum<=half)

        SELECT CASE (which)
        CASE (1)

! Solve for cum and ccum

          CALL local_cum_t(t,df,local_cum,local_ccum,status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)

! Solve for t

! Get start value of t
          t = dt1(local_cum,local_ccum,df)

          zf_status = 0

          CALL cdf_set_zero_finder(the_t,3,local)

          DO
            CALL rc_step_zf(zf_status,t,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_t(t,df,try_cum,try_ccum,status)

            IF (has_status) THEN
              IF (status/=0) THEN
                RETURN
              END IF
            END IF

            IF (match_cum) THEN
              fx = try_cum - local_cum
            ELSE
              fx = try_ccum - local_ccum
            END IF
          END DO

        CASE (3)

! Solve for df (degrees of freedom)

          df = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_t,4,local)

          DO
            CALL rc_step_zf(zf_status,df,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_t(t,df,try_cum,try_ccum,status)

            IF (has_status) THEN
              IF (status/=0) THEN
                RETURN
              END IF
            END IF

            IF (match_cum) THEN
              fx = try_cum - local_cum
            ELSE
              fx = try_ccum - local_ccum
            END IF

          END DO

        END SELECT

        IF (has_status) THEN
! Set the status of the zero finder
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_t

!*********************************************************************

      FUNCTION cum_t(t,df,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_t
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: df, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_t(which=1,cum=cum_t,t=t,df=df,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION cum_t

!*********************************************************************

      FUNCTION dt1(p,q,df)
!**********************************************************************
!     Double precision Initalize Approximation to
!           INVerse of the cumulative T distribution
!                              Function
!     Returns  the  inverse   of  the T   distribution   function, i.e.,
!     the integral from 0 to INVT of the T density is P. This is an
!     initial approximation.
!                              Arguments
!     P --> The p-value whose inverse from the T distribution is
!           desired.
!     Q --> 1-P.
!     DF --> Degrees of freedom of the T distribution.
!                              Method
!     Formula (26.7.5) from Abramowitz and Stegun is used to compute
!     the asymptotic expansion for the inverse function.
!**********************************************************************
! .. Use Statements ..
        USE cdf_normal_mod
! ..
! .. Function Return Value ..
        REAL (dpkind) :: dt1
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: df, p, q
! ..
! .. Local Scalars ..
        REAL (dpkind) :: denpow, sum, term, x, xp, xx
        INTEGER :: i, ideg
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, RESHAPE
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: coef(5,4) = RESHAPE((/one,one,zero, &
          zero,zero,three,16.0_dpkind,five,zero,zero,-15.0_dpkind, &
          17.0_dpkind,19.0_dpkind,three,zero,-945.0_dpkind, &
          -1920.0_dpkind,1482.0_dpkind,776.0_dpkind,79.0E0_dpkind/), &
          shape=(/5,4/))
        REAL (dpkind), PARAMETER :: denom(4) = (/ four, 96.0_dpkind, &
          384.0_dpkind, 92160.0_dpkind/)
! ..
        x = ABS(inv_normal(cum=p,check_input=.FALSE.))

        xx = x*x
        sum = x
        denpow = one

! DMS ideg(i) = i+1

        DO i = 1, 4
          ideg = i + 1
          term = x*evaluate_polynomial(coef(1:5,i),xx)
          denpow = denpow*df
          sum = sum + term/(denpow*denom(i))
        END DO

        IF (p<half) THEN
          xp = -sum
        ELSE
          xp = sum
        END IF

        dt1 = xp

        RETURN

      END FUNCTION dt1

!*********************************************************************

      FUNCTION inv_t(cum,ccum,df,status,check_input)
! Solves the INVERSE Poisson problem:
! Given cum (or/and ccum), and df, computes t.
! .. Function Return Value ..
        REAL (dpkind) :: inv_t
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind), INTENT (IN) :: df
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_t(which=2,cum=cum,ccum=ccum,t=inv_t,df=df,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION inv_t

!*********************************************************************

      SUBROUTINE local_cum_t(t,df,cum,ccum,status)
!----------------------------------------------------------------------

!                              Function

!     Computes the integral from -infinity to T of the t-density.
!     t --> Upper limit of integration of the t-density.
!     df --> Degrees of freedom of the t-distribution
!     cum<-- Cumulative t-distribution.
!     ccum <-- Compliment of Cumulative t-distribution.
!     status <-- status of the computation. 
!                Set to 10 if cdf_beta has no answer.

!                         Method

!     Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of
!     Mathematical Functions  (1966) is used to reduce the computation
!     of the cumulative distribution function to that of an incomplete
!     beta.
!------------------------------------------------------------------
! .. Use Statements ..
        USE cdf_beta_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: df, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a, dfptt, oma, tt, xx, yy
        INTEGER :: beta_status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        tt = t*t
        dfptt = df + tt
        xx = df/dfptt
        yy = tt/dfptt

        CALL cdf_beta(1,a,oma,xx,yy,half*df,half,status=beta_status)

        IF (beta_status/=0) THEN
! cdf_beta has NO answer.

          IF (PRESENT(status)) THEN
            status = 10

            RETURN
          ELSE
            WRITE (*,*) 'Error in local_cum_t call to cdf_beta'
            WRITE (*,*) 'Status: ', beta_status

            STOP 'Error in local_cum_t call to cdf_beta'
          END IF
        END IF

        IF (t<=zero) THEN
          cum = half*a
          ccum = oma + cum
        ELSE
          ccum = half*a
          cum = oma + ccum
        END IF

        RETURN

      END SUBROUTINE local_cum_t

!*********************************************************************

    END MODULE cdf_t_mod
