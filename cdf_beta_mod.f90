    MODULE cdf_beta_mod
! ----------------------------------------------------------------------

!                               cdf_beta_mod
!                               *=*=*=*=*=*=

!  -  SUBROUTINE CDF_BETA(WHICH, CUM, CCUM, X, CX, A, B, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_BETA(X, A, B, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_BETA(X, A, B, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_BETA(CUM, CCUM, A, B, STATUS, CHECK_INPUT)

!                             The Distribution
!                             ================

!    The density of the beta distribution is defined on x in [0,1] and is
!    proportional to:

!                                 a      b
!                                x  (1-x)

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    four arguments is to be calculated.
!  Input Range: [ 1:4 ]
!     1.  CUM and CCUM
!     2.  X and CX
!     3.  A
!     4.  B
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the beta distribution.
!  Range: [ 0:1 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the beta
!    distribution.
!  Range: [ 0:1 ]
!  - REAL (dpkind), OPTIONAL :: X. The upper limit of integration of the
!    beta density. The lower limit is 0.
!  Range: [ 0:1 ]
!  - REAL (dpkind), OPTIONAL :: CX. One minus the upper limit of
!    integration of the beta density. The lower limit is 0.
!  Range: [ 0:1 ]
!  - REAL (dpkind) :: A. The first parameter of the beta density.
!  Range: [ 10^-10:10^10 ]
!  - REAL (dpkind) :: B. The second parameter of the beta density.
!  Range: [ 10^-10:10^10 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code.  Possible values:
!      0 problem solved successfully
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 X outside range
!     -5 CX outside range
!     -6 A outside range
!     -7 B outside range
!      3 CUM + CCUM is not nearly one
!      4 X + CX is not nearly one
!    -50 Answer (if any) is BELOW the LOWER search bound
!     50 Answer (if any) is ABOVE the UPPER search bound
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.

!  NOTE: CUM and CCUM and also X and CX must add to (nearly) one.
! ----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_beta, cdf_beta, cum_beta, inv_beta
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_beta(x,cx,a,b,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_beta
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b
        REAL (dpkind), OPTIONAL, INTENT (IN) :: cx, x
        INTEGER, INTENT (OUT) :: status
        LOGICAL, INTENT (IN) :: check_input
! ..
        CALL cdf_beta(which=1,ccum=ccum_beta,x=x,cx=cx,a=a,b=b, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_beta

!*********************************************************************

      SUBROUTINE cdf_beta(which,cum,ccum,x,cx,a,b,status,check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: a, b
        REAL (dpkind), OPTIONAL :: ccum, cum, cx, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Structures ..
! .. Local Arrays
        TYPE (zf_locals) :: local
! ..
! .. Local Scalars ..
        REAL (dpkind) :: fx, local_ccum, local_cum, local_cx, local_x, &
          try_ccum, try_cum
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

! status = 0 means no error

        IF (has_status) THEN
          status = 0
        END IF

        CALL check_complements(cum,ccum,the_beta%name,'cum','ccum', &
          local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        CALL check_complements(x,cx,the_beta%name,'x','cx',local_x, &
          local_cx,set_values=(which/=2),bad_status=4,status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = local_x
        params(4) = local_cx
        params(5) = a
        params(6) = b

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

! ++++++++++          ++++++++++          ++++++++++
! Range check arguments and see that they add to one
! ++++++++++          ++++++++++          ++++++++++

        IF (local_check_input) THEN
! Assure that x + cx nearly one

          IF (which/=2) THEN
            IF ( .NOT. add_to_one(local_x,local_cx,the_beta%name,'x','cx' &
              ,4,status)) RETURN
          END IF

          CALL validate_parameters(the_beta,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF

        END IF

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute the Answers
!++++++++++++++++++++++++++++++++++++++++++++++++++

        IF (which>1) match_cum = (local_cum<=half)

        SELECT CASE (which)

        CASE (1)
! Calculate cum and ccum

          CALL local_cum_beta(local_x,local_cx,a,b,local_cum,local_ccum)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

          RETURN

        CASE (2)
! Calculate x and cx

          local_x = half
          zf_status = 0

          CALL cdf_set_zero_finder(the_beta,3,local)

          IF (match_cum) THEN
            DO
              CALL rc_interval_zf(zf_status,local_x,fx,local)

              IF (zf_status/=1) EXIT

              local_cx = one - local_x

              CALL local_cum_beta(local_x,local_cx,a,b,try_cum,try_ccum)

              fx = try_cum - cum
            END DO

          ELSE
            DO
              CALL rc_interval_zf(zf_status,local_cx,fx,local)

              IF (zf_status/=1) EXIT

              local_x = one - local_cx

              CALL local_cum_beta(local_x,local_cx,a,b,try_cum,try_ccum)

              fx = try_ccum - ccum
            END DO
          END IF

          IF (PRESENT(x)) x = local_x
          IF (PRESENT(cx)) cx = local_cx

        CASE (3)
! Calculate a

          a = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_beta,5,local)

          DO
            CALL rc_step_zf(zf_status,a,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_beta(local_x,local_cx,a,b,try_cum,try_ccum)

            IF (match_cum) THEN
              fx = try_cum - cum
            ELSE
              fx = try_ccum - ccum
            END IF
          END DO

        CASE (4)
! Calculate b

          b = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_beta,6,local)

          DO
            CALL rc_step_zf(zf_status,b,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_beta(local_x,local_cx,a,b,try_cum,try_ccum)

            IF (match_cum) THEN
              fx = try_cum - cum
            ELSE
              fx = try_ccum - ccum
            END IF
          END DO

        END SELECT

        IF (has_status) THEN
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_beta

!*********************************************************************

      FUNCTION cum_beta(x,cx,a,b,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_beta
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b
        REAL (dpkind), OPTIONAL, INTENT (IN) :: cx, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_beta(which=1,cum=cum_beta,x=x,cx=cx,a=a,b=b, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_beta

!*********************************************************************

      FUNCTION inv_beta(cum,ccum,a,b,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_beta
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b
        REAL (dpkind), OPTIONAL, INTENT (IN) :: ccum, cum
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_beta(which=2,cum=cum,ccum=ccum,x=inv_beta,a=a,b=b, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_beta

!*********************************************************************

      SUBROUTINE local_cum_beta(x,y,a,b,cum,ccum)
!----------------------------------------------------------------------
!          Double precision cUMulative incomplete BETa distribution

!                              Function

!     Calculates the cdf to X of the incomplete beta distribution
!     with parameters a and b.  This is the integral from 0 to x
!     of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)

!                              Arguments

!     X --> Upper limit of integration.
!     Y --> 1 - X.
!     A --> First parameter of the beta distribution.
!     B --> Second parameter of the beta distribution.
!     CUM <-- Cumulative incomplete beta distribution.
!     CCUM <-- Compliment of Cumulative incomplete beta distribution.

!                              Method

!     Cumulative distribution function  (CUM)  is calculated directly by
!     code associated with the following reference.
!     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
!     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
!     Trans. Math.  Softw. 18 (1993), 360-373.
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of CUM.   The search relies  on  the
!     monotinicity of CUM with the other parameter.
!----------------------------------------------------------------------
! .. Use Statements ..
        USE biomath_mathlib_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, x, y
        REAL (dpkind), INTENT (OUT) :: ccum, cum
! ..
! .. Local Scalars ..
        INTEGER :: ierr
! ..
        IF (x<=zero) THEN
          cum = zero
          ccum = one
          RETURN
        END IF

        IF (y<=zero) THEN
          cum = one
          ccum = zero
          RETURN
        END IF

! Because of the bounds on a, b, x, y
! ierr can not be <> 0

        CALL bratio(a,b,x,y,cum,ccum,ierr)

        RETURN

      END SUBROUTINE local_cum_beta

!*********************************************************************

    END MODULE cdf_beta_mod
