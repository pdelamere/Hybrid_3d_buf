    MODULE cdf_binomial_mod
! ______________________________________________________________________
!
!                             cdf_binomial_mod
!                             *=*=*=*=*=*=*=*=
!
!  -  SUBROUTINE CDF_BINOMIAL(WHICH, CUM, CCUM, S, N, PR, CPR, STATUS,
!                             CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_BINOMIAL(S, N, PR, CPR, STATUS, CHECK_INPUT)
!  -  REAL (dpkind)FUNCTION CCUM_BINOMIAL(S, N, PR, CPR, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_BINOMIAL(CUM, CCUM, N, PR, CPR, STATUS,
!                                         CHECK_INPUT)
!
!                             The Distribution
!                             ================
!
!   The density of the binomial distribution provides the probability of S
!   successes in N independent trials, each with probability of success
!   PR.  The density is proportional to:
! 
!                                S       N-S 
!                              PR  (1-PR)    
!
!     The binomial is extended to non-integer values via the connection
!     between the cumulative binomial and the incomplete beta.
!
!                                Arguments
!                                =========
!
!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    four arguments is to be calculated.
!  Input Range: [ 1:4 ]
!     - 1 CUM and CCUM 
!     - 2 S 
!     - 3 N 
!     - 4 PR and CPR 
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the binomial distribution,
!    i.e., the probability of 0 to S successes in N trials.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the binomial
!    distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: S. The upper limit of summation of the binomial
!    density. Note that S must be less than or equal to N.
!  Range: [ 0:10^10 ]
!  - REAL (dpkind) :: N. The number of independent trials generating the
!    binomial density. N must be greater than or equal to S.
!  Range: [ 0:10^10 ]
!  - REAL (dpkind), OPTIONAL :: PR. The probability of success in each
!    independent trial.
!  Range: [ 0:1 ]
!  - REAL (dpkind), OPTIONAL :: CPR. One minus the probability of success
!    in each independent trial; the probability of failure in each trial.
!  Range: [ 0:1 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code.
!        0 problem solved successfully 
!       -1 WHICH outside input range 
!       -2 CUM outside range 
!       -3 CCUM outside range 
!       -4 S outside range 
!       -5 N outside range 
!       -6 PR outside range 
!       -7 CPR outside range 
!        3 CUM + CCUM is not nearly one 
!        4 PR + CPR is not nearly one 
!        5 S not between 0 and N.  
!       10 cdf_beta (in local_cum_binomial) has no answer.
!      -50 Answer (if any) is BELOW the LOWER search bound 
!       50 Answer (if any) is ABOVE the UPPER search bound 
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.
!
!  NOTE: CUM and CCUM and also PR and CPR must add to (nearly) one.
! ______________________________________________________________________
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Implicit None Statement ..
      IMPLICIT NONE
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_binomial, cdf_binomial, cum_binomial, inv_binomial
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_binomial(s,n,pr,cpr,status,check_input)
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: cpr, pr
        REAL (dpkind) :: n, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Function Return Value ..
        REAL (dpkind) :: ccum_binomial
! ..
! .. Executable Statements ..

        CALL cdf_binomial(which=1,ccum=ccum_binomial,s=s,n=n,pr=pr, &
          cpr=cpr,status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_binomial

!*********************************************************************

      SUBROUTINE cdf_binomial(which,cum,ccum,s,n,pr,cpr,status, &
          check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE zero_finder
! ..
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cpr, cum, pr
        REAL (dpkind) :: n, s
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
        INTRINSIC MAX, PRESENT
! ..
! .. Executable Statements ..

        
        has_status = PRESENT(status)

        ! status = 0 means NO error

        IF (has_status) THEN
           status = 0
        END IF

        ! Check that at least one of cum, ccum is present in the calling list

        CALL check_complements(cum,ccum,the_binomial%name,'cum', 'ccum', &
             local_cum,local_ccum,set_values=(which/=1), bad_status=3, &
             status=status)

        IF (has_status) THEN
           IF (status /= 0) THEN
              RETURN
           END IF
        END IF

        ! Check that at least one of pr, cpr is present in the calling list

        CALL check_complements(pr,cpr,the_binomial%name,'pr','cpr', &
             local_pr,local_cpr,set_values=(which/=4), bad_status=4, &
             status=status)

        IF (has_status) THEN
           IF (status /= 0) THEN
              RETURN
           END IF
        END IF

        !==========
        ! Check that probabilities add to one
        !==========

        IF (which /= 4) THEN
           IF ( .NOT. add_to_one(local_pr,local_cpr,the_binomial%name, &
                'pr','cpr',4,status)) RETURN
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = s
        params(4) = n
        params(5) = local_pr
        params(6) = local_cpr

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN

          CALL validate_parameters(the_binomial,which,params,status)

          IF (has_status) THEN
             IF (status /= 0) THEN
                RETURN
             END IF
          END IF

          IF (which/=2 .AND. which/=3) THEN

            ! Assure that s is between 0 and n

            IF ( .NOT. in_range(value=s, lo=zero, hi=n, routine_name= &
              the_binomial%name, arg_name='s', bad_status=5,status=status)) &

              RETURN
          END IF

        END IF

        !==========
        ! Calculate answers
        !==========

        IF (which>1) match_cum = (local_cum <= half)

        SELECT CASE (which)

        CASE (1)
           ! Solve for the cumulative binomial (cdf)

           CALL local_cum_binomial(s,n,local_pr,local_cpr,local_cum, &
                local_ccum, status)

           IF (has_status) THEN
              IF (status /= 0) THEN
                 ! Problem can NOT be solved because cdf_beta has NO solution.
                 RETURN
              END IF
           END IF

           IF (PRESENT(cum)) cum = local_cum
           IF (PRESENT(ccum)) ccum = local_ccum

           RETURN

        CASE (2)
           ! Solve for s (number of successes)

          s = half * n
          zf_status = 0

! DMS
! Clearly s <= n, so the UPPER bound for s should be n
! We can not call cdf_set_zero() as in the other cases, because
! the_binomial structures contains PARAMETERS which can NOT be modified.
! The values are set like in  cdf_set_zero()

          CALL set_zero_finder( &
               low_limit=the_binomial%parameters(3)%low_bound, &
               hi_limit=n, abs_step=half, rel_step=half, &
               step_multiplier=five, abs_tol=1.0E-50_dpkind, &
               rel_tol=1.0E-8_dpkind, local=local)

          DO
!            CALL rc_interval_zf(zf_status, s, fx, local)
             CALL rc_step_zf(zf_status, s, fx, local)

            IF (zf_status /= 1) EXIT

            CALL local_cum_binomial(s, n, local_pr, local_cpr, try_cum, &
                                    try_ccum, status)

            IF (has_status) THEN
               IF (status /= 0) THEN
                  RETURN
               END IF
            END IF

            IF (match_cum) THEN
              fx = try_cum - local_cum
            ELSE
              fx = local_ccum - try_ccum
            END IF

          END DO

        CASE (3)
          ! Solve for n (sample size)

          n = two * MAX(s,one)
          zf_status = 0

          CALL cdf_set_zero_finder(the_binomial,4,local)

          DO
            CALL rc_step_zf(zf_status,n,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_binomial(s,n,local_pr,local_cpr,try_cum, &
                                    try_ccum, status)

            IF (has_status) THEN
               IF (status /= 0) THEN
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
           ! Solve for the binomial probability (pr) or
           ! its complement (cpr)

          zf_status = 0

!!! A difficult  case -- we  will try to  match the lesser of  cum and
!!! ccum If the answer <= 1/2 we will vary pr else we will vary cpr

! Decide whether to vary pr or cpr

          CALL local_cum_binomial(s,n,half,half,try_cum,try_ccum, status)

          IF (has_status) THEN
             IF (status /= 0) THEN
                RETURN
             END IF
          END IF

          vary_pr = (try_cum<=local_cum)

          IF (match_cum .AND. vary_pr) THEN

             ! We are trying to match cum varying pr

            CALL cdf_set_zero_finder(the_dummy_binomial,1,local)

            DO
              CALL rc_interval_zf(zf_status,try_pr,fx,local)

              try_cpr = one - try_pr

              IF (zf_status/=1) EXIT

              CALL local_cum_binomial(s,n,try_pr,try_cpr,try_cum, &
                                      try_ccum, status)

              IF (has_status) THEN
                 IF (status /= 0) THEN
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

              try_pr = one - try_cpr

              IF (zf_status/=1) EXIT

              CALL local_cum_binomial(s,n,try_pr,try_cpr,try_cum, &
                                      try_ccum, status)

              IF (has_status) THEN
                 IF (status /= 0) THEN
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

              try_cpr = one - try_pr

              IF (zf_status/=1) EXIT

              CALL local_cum_binomial(s,n,try_pr,try_cpr,try_cum, &
                                      try_ccum, status)

              IF (has_status) THEN
                 IF (status /= 0) THEN
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

              try_pr = one - try_cpr

              IF (zf_status/=1) EXIT


              CALL local_cum_binomial(s,n,try_pr,try_cpr,try_cum, &
                                      try_ccum, status)

              IF (has_status) THEN
                 IF (status /= 0) THEN
                    RETURN
                 END IF
              END IF

              fx = local_ccum - try_ccum
            END DO

          END IF

          IF (PRESENT(pr)) pr = try_pr

          IF (PRESENT(cpr)) cpr = try_cpr

        END SELECT

        IF (has_status) THEN
           ! Set the status of the zero finder
           CALL cdf_finalize_status( local, status )
        END IF

        RETURN

      END SUBROUTINE cdf_binomial

!*********************************************************************

      FUNCTION cum_binomial(s,n,pr,cpr,status,check_input)
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: cpr, pr
        REAL (dpkind) :: n, s
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Function Return Value ..
        REAL (dpkind) :: cum_binomial
! ..
! .. Executable Statements ..

        CALL cdf_binomial(which=1,cum=cum_binomial,s=s,n=n,pr=pr,cpr=cpr, &
                          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_binomial

!*********************************************************************

      FUNCTION inv_binomial(cum,ccum,n,pr,cpr,status,check_input)
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cpr, cum, pr
        REAL (dpkind) :: n
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Function Return Value ..
        REAL (dpkind) :: inv_binomial
! ..
! .. Executable Statements ..

        CALL cdf_binomial(2,cum,ccum,inv_binomial,n,pr,cpr,status, &
                          check_input)

        RETURN

      END FUNCTION inv_binomial

!*********************************************************************

      SUBROUTINE local_cum_binomial(s,xn,pr,ompr,cum,ccum, status)
!----------------------------------------------------------------------
!                    CUMulative BINomial distribution
!
!                              Function
!
!     Returns the probability   of 0  to  S  successes in  XN   binomial
!     trials, each of which has a probability of success, PBIN.
!                              Arguments
!     S --> The upper limit of cumulation of the binomial distribution.
!
!     XN --> The number of binomial trials.
!
!     PR --> The probability of success in each binomial trial.
!
!     OMPR --> 1 - PBIN
!
!     CUM <-- Cumulative binomial distribution.
!
!     CCUM <-- Complement of Cumulative binomial distribution.
!
!     STATUS <-- status of the computation. 
!                Set to 10 if cdf_beta has no answer.
!
!                              Method
!
!     Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
!     Mathematical   Functions (1966) is   used  to reduce the  binomial
!     distribution  to  the  cumulative    beta distribution.
!----------------------------------------------------------------------
! .. Use Statements ..
        USE cdf_beta_mod
! ..
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: ompr, pr, s, xn
        INTEGER, OPTIONAL, INTENT(OUT) :: status
! ..
! .. Local Scalars
        INTEGER :: beta_status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! .. Executable Statements ..


        IF (s < xn) THEN
          CALL cdf_beta(1,ccum,cum,pr,ompr,s+one,xn-s, status=beta_status)

          IF (beta_status /= 0) THEN
             ! cdf_beta has NO answer.

             IF (PRESENT(status)) THEN
                status = 10
             ELSE
                WRITE (*,*) 'Error in local_cum_binomial call to cdf_beta'
                WRITE (*,*) 'Status: ', beta_status

                STOP 'Error in local_cum_binomial call to cdf_beta'
             END IF
          END IF
        ELSE
          cum = one
          ccum = zero
        END IF

        RETURN

      END SUBROUTINE local_cum_binomial

!*********************************************************************

    END MODULE cdf_binomial_mod
