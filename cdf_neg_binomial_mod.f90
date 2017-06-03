    MODULE cdf_neg_binomial_mod
! ----------------------------------------------------------------------

!                           cdf_neg_binomial_mod
!                           *=*=*=*=*=*=*=*=*=*=

!  -  SUBROUTINE CDF_NEG_BINOMIAL(WHICH, CUM, CCUM, F, S, PR, CPR,
!                                 STATUS, BCHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_NEG_BINOMIAL(F, S, N, PR, CPR, STATUS,
!                                             CHECK_INPUT)
!  -  REAL (dpkind)FUNCTION CCUM_NEG_BINOMIAL(F, S, N, PR, CPR, STATUS,
!                                             CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_NEG_BINOMIAL(CUM, CCUM, N, PR, CPR,
!                                             STATUS, CHECK_INPUT)

!                             The Distribution
!                             ================

!   The density of the negative binomial distribution provides the
!   probability of precisely F failures before the S'th success in
!   independent binomial trials, each with probability of success PR.
!   The density is:

!                          (  F+S-1  )  S       F
!                          (   S-1   )PR  (1-PR)

!   The cumulative distribution function is the probability of F or fewer
!   failures before the F'th success.The negative binomial is extended to
!   non-integer values of F via the relation between the cumulative
!   distribution function of the negative binomial and the incomplete beta
!   function.
!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    four arguments is to be calculated.
!  Input Range: [ 1:4 ]
!     - 1 CUM and CCUM
!     - 2 F
!     - 3 S
!     - 4 PR and CPR
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the negative-binomial
!    distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL : CCUM. One minus the CDF of the binomial
!    distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: F. The number of failures before the S'th success.
!  Range: [ 0: ]
!  - REAL (dpkind) :: S. The number of successes to occur.
!  Input Range: [ 0: ]
!  - REAL (dpkind) :: PR. The probability of success in each independent
!    trial.
!  Range: [ 0:1 ]
!  - REAL (dpkind) :: CPR. One minus the probability of success in each
!    independent trial; the probability of failure in each trial.
!  Range: [ 0:1 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!        0 problem successfully solved
!       -1 WHICH outside input range
!       -2 CUM outside range
!       -3 CCUM outside range
!       -4 F outside range
!       -5 S outside range
!       -6 PR outside range
!       -7 CPR outside range
!        3 CUM + CCUM is not nearly one
!        4 PR + CPR is not nearly one
!       10 cdf_beta (in local_cum_neg_binomial) has no answer.
!      -50 Answer (if any) is BELOW the LOWER search bound
!       50 Answer (if any) is ABOVE the UPPER search bound
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!        .TRUE. input argument values are not checked for validity.

!  NOTE: CUM and CCUM and also PR and CPR must add to (nearly) one.
! ----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_neg_binomial, cdf_neg_binomial, cum_neg_binomial, &
        inv_neg_binomial
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_neg_binomial(f,s,pr,cpr,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_neg_binomial
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: cpr, pr
        REAL (dpkind), INTENT (IN) :: f, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_neg_binomial(which=1,cum=ccum_neg_binomial,f=f,s=s, &
          pr=pr,cpr=cpr,status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_neg_binomial

!*********************************************************************

      SUBROUTINE cdf_neg_binomial(which,cum,ccum,f,s,pr,cpr,status, &
          check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cpr, cum, pr
        REAL (dpkind) :: f, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Structures ..
        TYPE (zf_locals) :: local
! ..
! .. Local Arrays ..
        REAL (dpkind) :: params(6)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: fx, local_ccum, local_cpr, local_cum, local_pr, &
          try_ccum, try_cpr, try_cum, try_pr
        INTEGER :: zf_status
        LOGICAL :: has_status, local_check_input, match_cum, vary_pr
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        has_status = PRESENT(status)

! status = 0 means NO error

        IF (has_status) THEN
          status = 0
        END IF

! Check that at least one of cum, ccum is present in the calling list

        CALL check_complements(cum,ccum,the_negative_binomial%name,'cum', &
          'ccum',local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        CALL check_complements(pr,cpr,the_negative_binomial%name,'pr', &
          'cpr',local_pr,local_cpr,set_values=(which/=4),bad_status=4, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

!==========
! Check that probabilities add to one
!==========

        IF (which/=4) THEN
          IF ( .NOT. add_to_one(local_pr,local_cpr,the_binomial%name,'pr' &
            ,'cpr',4,status)) RETURN
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = f
        params(4) = s
        params(5) = local_pr
        params(6) = local_cpr

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_negative_binomial,which,params, &
            status)

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
! Calculate cum and ccum (cdf)

          CALL local_cum_neg_binomial(f,s,local_pr,local_cpr,local_cum, &
            local_ccum,status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)
! Calculate f

          f = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_negative_binomial,3,local)

          DO
            CALL rc_step_zf(zf_status,f,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_neg_binomial(f,s,local_pr,local_cpr,try_cum, &
              try_ccum,status)

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
! Calculate s

          s = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_negative_binomial,4,local)

          DO
            CALL rc_step_zf(zf_status,s,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_neg_binomial(f,s,local_pr,local_cpr,try_cum, &
              try_ccum,status)

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

        CASE (4)

! Solve for probability

!!! A difficult case -- we will try to match the lesser of cum and ccum
!!! If the answer <= 1/2 we will vary pr else we will vary cpr

          zf_status = 0

! Decide whether to vary pr or cpr

          CALL local_cum_neg_binomial(f,s,half,half,try_cum,try_ccum, &
            status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF

! Opposite than for BINOMIAL distribution!!!

          vary_pr = (try_cum>local_cum)

          IF (match_cum .AND. vary_pr) THEN

! We are trying to match cum varying pr

            CALL cdf_set_zero_finder(the_dummy_binomial,1,local)

            DO
              CALL rc_interval_zf(zf_status,try_pr,fx,local)

              IF (zf_status/=1) EXIT

              try_cpr = one - try_pr

              CALL local_cum_neg_binomial(f,s,try_pr,try_cpr,try_cum, &
                try_ccum,status)

              IF (has_status) THEN
                IF (status/=0) THEN
                  RETURN
                END IF
              END IF

              fx = local_cum - try_cum
            END DO

          ELSE IF (match_cum .AND. .NOT. vary_pr) THEN
! Try to match cum varying cpr

            CALL cdf_set_zero_finder(the_dummy_binomial,2,local)

            DO

              CALL rc_interval_zf(zf_status,try_cpr,fx,local)

              IF (zf_status/=1) EXIT

              try_pr = one - try_cpr

              CALL local_cum_neg_binomial(f,s,try_pr,try_cpr,try_cum, &
                try_ccum,status)

              IF (has_status) THEN
                IF (status/=0) THEN
                  RETURN
                END IF
              END IF

              fx = local_cum - try_cum
            END DO

          ELSE IF ( .NOT. match_cum .AND. vary_pr) THEN
! Try to match ccum varying pr

            CALL cdf_set_zero_finder(the_dummy_binomial,1,local)

            DO
              CALL rc_interval_zf(zf_status,try_pr,fx,local)

              IF (zf_status/=1) EXIT

              try_cpr = one - try_pr

              CALL local_cum_neg_binomial(f,s,try_pr,try_cpr,try_cum, &
                try_ccum,status)

              IF (has_status) THEN
                IF (status/=0) THEN
                  RETURN
                END IF
              END IF

              fx = local_ccum - try_ccum

            END DO

          ELSE

! Try to match ccum varying cpr

            CALL cdf_set_zero_finder(the_dummy_binomial,2,local)

            DO
              CALL rc_interval_zf(zf_status,try_cpr,fx,local)

              IF (zf_status/=1) THEN
                IF (PRESENT(pr)) pr = try_pr
                IF (PRESENT(cpr)) cpr = try_cpr

                EXIT
              END IF

              try_pr = one - try_cpr

              CALL local_cum_neg_binomial(f,s,try_pr,try_cpr,try_cum, &
                try_ccum,status)

              IF (has_status) THEN
                IF (status/=0) THEN
                  RETURN
                END IF
              END IF

              fx = local_ccum - try_ccum
            END DO

          END IF

        END SELECT

        IF (has_status) THEN
! Set the status of the zero finder
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_neg_binomial

!*********************************************************************

      FUNCTION cum_neg_binomial(f,s,pr,cpr,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_neg_binomial
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: cpr, pr
        REAL (dpkind), INTENT (IN) :: f, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_neg_binomial(which=1,cum=cum_neg_binomial,f=f,s=s,pr=pr, &
          cpr=cpr,status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_neg_binomial

!*********************************************************************

      FUNCTION inv_neg_binomial(cum,ccum,s,pr,cpr,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_neg_binomial
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: ccum, cpr, cum, pr
        REAL (dpkind), INTENT (IN) :: s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_neg_binomial(which=2,cum=cum,ccum=ccum, &
          f=inv_neg_binomial,s=s,pr=pr,cpr=cpr,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION inv_neg_binomial

!*********************************************************************

      SUBROUTINE local_cum_neg_binomial(f,s,pr,cpr,cum,ccum,status)

!                              METHOD

!     Formula  26.5.26    of   Abramowitz  and    Stegun,  Handbook   of
!     Mathematical   Functions (1966) is   used  to reduce the  negative
!     binomial distribution to the cumulative beta distribution.
! .. Use Statements ..
        USE cdf_beta_mod
! ..
! .. Scalar Arguments ..
! .. Local Scalars
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: cpr, f, pr, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! .. Local Scalars ..
        INTEGER :: beta_status
! ..
        CALL cdf_beta(1,cum,ccum,x=pr,cx=cpr,a=s,b=f+one, &
          check_input=.FALSE.,status=beta_status)

        IF (beta_status/=0) THEN
! cdf_beta has NO answer.

          IF (PRESENT(status)) THEN
            status = 10
          ELSE
            WRITE (*,*) &
              'Error in local_cum_neg_binomial call to cdf_beta'
            WRITE (*,*) 'Status: ', beta_status

            STOP 'Error in local_cum_neg_binomial call to cdf_beta'
          END IF
        END IF

        RETURN

      END SUBROUTINE local_cum_neg_binomial

!*********************************************************************

    END MODULE cdf_neg_binomial_mod
