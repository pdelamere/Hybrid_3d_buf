    MODULE cdf_chisq_mod
! ______________________________________________________________________

!                              cdf_chisq_mod
!                              *=*=*=*=*=*=*

!  -  SUBROUTINE CDF_CHISQ( WHICH, CUM, CCUM, X, DF, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_CHISQ(X, DF, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_CHISQ(X, DF, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_CHISQ( CUM, CCUM, DF, STATUS, CHECK_INPUT)

!                             The Distribution
!                             ================

!   The chi-squared distribution is the distribution of the sum of squares
!   of DF independent unit (mean=0, sd=1) normal deviates.The density is
!              defined on x in [0, +Infinity) and is proportional to :

!                            (DF-2)/2           
!                           x         exp(-x/2)

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    three arguments is to be calculated.
!  Input Range: [ 1:3 ]
!     1.  CUM and CCUM 
!     2.  X  
!     3.  DF 
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the chi-squared 
!    distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the chi-squared
!    distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: X. The upper limit of integration of the chi-squared
!    density. The lower limit is 0.
!  Range: [ 0:10^100 ]
!  - REAL (dpkind) :: DF. The degrees of freedom of the chi-squared
!    distribution.
!  Range: [ 10^-3:10^10 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!      0 problem successfully solved
!     -1 WHICH outside input range 
!     -2 CUM outside range 
!     -3 CCUM outside range 
!     -4 X outside range 
!     -5 DF outside range 
!      3 CUM + CCUM is not nearly one 
!      4 X + CX is not nearly one 
!     10 Problem could NOT be solved because
!           cdf_gamma (called in local_cum_chisq) has NO solution.
!    -50 Answer (if any) is BELOW the LOWER search bound 
!     50 Answer (if any) is ABOVE the UPPER search bound 
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.

!       NOTE: CUM and CCUM MUST add to (nearly) one.
! ______________________________________________________________________
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_chisq, cdf_chisq, cum_chisq, inv_chisq
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_chisq(x,df,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_chisq
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: df, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_chisq(which=1,ccum=ccum_chisq,x=x,df=df,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION ccum_chisq

!*********************************************************************

      SUBROUTINE cdf_chisq(which,cum,ccum,x,df,status,check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind), INTENT (INOUT) :: df, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Structures ..
        TYPE (zf_locals) :: local
! ..
! .. Local Scalars ..
! .. Local Arrays
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

        CALL check_complements(cum,ccum,the_chi_square%name,'cum','ccum', &
          local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = x
        params(4) = df

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_chi_square,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF

        END IF

!==========
! Calculate answers
!==========

        IF (which>1) match_cum = (local_cum<=half)

        SELECT CASE (which)

        CASE (1) ! Solve for CDF
          CALL local_cum_chisq(x,df,local_cum,local_ccum,status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2) ! Solve for x
          x = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_chi_square,3,local)

          DO
            CALL rc_step_zf(zf_status,x,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_chisq(x,df,try_cum,try_ccum,status)

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

        CASE (3) ! Solve for df
          df = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_chi_square,4,local)

          DO
            CALL rc_step_zf(zf_status,df,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_chisq(x,df,try_cum,try_ccum,status)

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
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_chisq

!*********************************************************************

      FUNCTION cum_chisq(x,df,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_chisq
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: df, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_chisq(which=1,cum=cum_chisq,x=x,df=df,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION cum_chisq

!*********************************************************************

      FUNCTION inv_chisq(cum,ccum,df,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_chisq
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: df
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_chisq(2,cum,ccum,inv_chisq,df,status,check_input)

        RETURN

      END FUNCTION inv_chisq

!*********************************************************************

      SUBROUTINE local_cum_chisq(x,df,cum,ccum,status)
!----------------------------------------------------------------------
!                    Cumulative Chi Square distribution

!                              Function

!     Computes the integral from 0 to x for the Chi-square distribution.

!                              Arguments

!     x    --> upper limit of integration
!     df   --> degrees of freedom
!     cum  <-- Cumulative chi-square distribution.
!     ccum <-- Complement of Cumulative distribution.

!     STATUS <-- status of the computation. 
!                Set to 10 if cdf_gamma has no answer.

!                          Method

!     Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of
!     Mathematical Functions   (1966) is used   to reduce the chisqure
!     distribution to the incomplete gamma distribution.
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotinicity of P with the other parameter.
!----------------------------------------------------------------------
! .. Use Statements ..
        USE cdf_gamma_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: df, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        INTEGER :: gamma_status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        CALL cdf_gamma(1,cum,ccum,half*x,half*df,one,status=gamma_status)

        IF (gamma_status/=0) THEN
! cdf_gamma has NO answer.

          IF (PRESENT(status)) THEN
            status = 10
          ELSE
            WRITE (*,*) 'Error in local_cum_chisq call to cdf_gamma'
            WRITE (*,*) 'Status: ', gamma_status

            STOP 'Error in local_cum_chisq call to cdf_gamma'
          END IF
        END IF

        RETURN

      END SUBROUTINE local_cum_chisq

!*********************************************************************

    END MODULE cdf_chisq_mod
