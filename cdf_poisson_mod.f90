    MODULE cdf_poisson_mod
! ----------------------------------------------------------------------

!                             cdf_poisson_mod

!                             *=*=*=*=*=*=*=*
!  -  SUBROUTINE CDF_POISSON(WHICH, CUM, CCUM, S, LAMBDA, STATUS,
!                            CHECK_INPUT )
!  -  REAL (dpkind) FUNCTION CUM_POISSON(S, LAMBDA, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_POISSON(S, LAMBDA, STATUS, CHECK_INPUT)
!  -  FUNCTION INV_POISSON(CUM, CCUM, LAMBDA, STATUS, CHECK_INPUT)

!                             The Distribution
!                             ================

!    The density of the Poisson distribution (probability of observing S
!    events) is:
!                                 S
!                           LAMBDA
!                           ------- exp(-LAMBDA)
!                             S!

!   The Poisson distribution is extended to non-integer values of S using
!   the relation between the cumulative distribution function of the
!   Poisson distribution and the gamma distribution.

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    arguments is to be calculated.
!  Input Range: [ 1:3 ]
!     1.  CUM and CCUM
!     2.  S
!     3.  LAMBDA
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the Poisson distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the Poisson
!    distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: S. The upper limit of summation of the Poisson
!    density. The lower limit is 0.
!  Range: [ 0:10^100 ]
!  - REAL (dpkind) :: LAMBDA. The mean of the Poisson distribution.
!    Range: [ 10^-10:10^100 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Posible values:
!      0 Problem solved successfully
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 X outside range
!     -5 MEAN outside range
!     -6 SD outside range
!      3 CUM + CCUM is not nearly one
!     10 could not solve cdf_gamma in local_cum_poisson
!    -50 Answer (if any) is BELOW the LOWER bound on range
!     50 Answer (if any) is ABOVE the UPPER bound on range
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and .TRUE.
!    input argument values are not checked for validity.
! ----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_poisson, cdf_poisson, cum_poisson, inv_poisson
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_poisson(s,lambda,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_poisson
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: lambda, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_poisson(which=1,ccum=ccum_poisson,s=s,lambda=lambda, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_poisson

!*********************************************************************

      SUBROUTINE cdf_poisson(which,cum,ccum,s,lambda,status,check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE cdf_gamma_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: lambda, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Arrays ..
        REAL (dpkind) :: params(6)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: local_ccum, local_cum, s_p_1
        LOGICAL :: has_status, local_check_input
! ..
! .. Intrinsic Functions ..
        INTRINSIC MAX, PRESENT
! ..
        has_status = PRESENT(status)

! status = 0 means NO error

        IF (has_status) THEN
          status = 0
        END IF

! Check that at least one of cum, ccum is present in the calling list

        CALL check_complements(cum,ccum,the_poisson%name,'cum','ccum', &
          local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = s
        params(4) = lambda

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_poisson,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF
        END IF

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute the Answers
!++++++++++++++++++++++++++++++++++++++++++++++++++

        SELECT CASE (which)

        CASE (1)
! Calculate cum and ccum

          CALL cdf_gamma(1,local_ccum,local_cum,x=lambda,shape=s+one, &
            scale=one,status=status,check_input=.FALSE.)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)
! Calculate s

          CALL cdf_gamma(3,local_ccum,local_cum,x=lambda,shape=s_p_1, &
            scale=one,status=status,check_input=.FALSE.)

          s = MAX(s_p_1-one,zero)

        CASE (3)
! Calculate lambda

          CALL cdf_gamma(2,local_ccum,local_cum,x=lambda,shape=s+one, &
            scale=one,status=status,check_input=.FALSE.)

        END SELECT

        RETURN

      END SUBROUTINE cdf_poisson

!*********************************************************************

      FUNCTION cum_poisson(s,lambda,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_poisson
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: lambda, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_poisson(which=1,cum=cum_poisson,s=s,lambda=lambda, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_poisson

!*********************************************************************

      FUNCTION inv_poisson(cum,ccum,lambda,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_poisson
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: lambda
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_poisson(which=2,cum=cum,ccum=ccum,s=inv_poisson, &
          lambda=lambda,status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_poisson

!*********************************************************************

    END MODULE cdf_poisson_mod
