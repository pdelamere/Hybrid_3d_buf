    MODULE cdf_nc_chisq_mod
!----------------------------------------------------------------------

!                             cdf_nc_chisq_mod
!                             *=*=*=*=*=*=*=*=

!  -  SUBROUTINE CDF_NC_CHISQ(WHICH, CUM, CCUM, X, DF, PNONC, STATUS,
!                             CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_NC_CHISQ(X, DF, PNONC, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_NC_CHISQ(X, DF, PNONC, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_NC_CHISQ(CUM, CCUM, DF, PNONC, STATUS,
!                                         CHECK_INPUT)

!                             The Distribution
!                             ================

!    The noncentral chi-squared distribution is the sum of DF independent
!     normal distributions with unit standard deviations, but possibly
! non-zero means . Let the mean of the ith normal be delta_i. Then PNONC =
!                             Sigma_i delta_i.

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    three arguments is to be calculated.
!     1.  CUM and CCUM
!     2.  DF
!     3.  PNONC
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the noncentral chi-squared
!    distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the noncentral
!    chi-squared distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: X. The upper limit of integration of the noncentral
!    chi-squared distribution. (The lower limit is 0.) Range: [ 0:10^100 ]
!  - REAL (dpkind) :: DF. The degrees of freedom of the noncentral
!    chi-squared distribution.
!  Input Range: [ 10^-3:10^10 ]
!  - REAL (dpkind) :: PNONC. The noncentrality parameter.
!  Range: [ 0:10^4 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!      0 Problem solved successfully 
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 X outside range
!     -5 DF outside range
!     -6 PNONC outside range
!      3 CUM + CCUM is not nearly one
!     10 no solution for cdf_chisq called in local_cum_nc_chisq
!    -50 Answer (if any) is below the lower search bound
!     50 Answer (if any) is above the upper search bound
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.

!  NOTE: CUM and CCUM must add to (nearly) one.
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
      PUBLIC :: ccum_nc_chisq, cdf_nc_chisq, cum_nc_chisq, inv_nc_chisq
! ..
! .. Parameters ..
      REAL (dpkind), PARAMETER :: ten4 = 1.0E4_dpkind
      REAL (dpkind), PARAMETER :: tiny = 1.0E-100_dpkind
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_nc_chisq(x,df,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_nc_chisq
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: df, pnonc, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_chisq(which=1,ccum=ccum_nc_chisq,x=x,df=df, &
          pnonc=pnonc,status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_nc_chisq

!*********************************************************************

      SUBROUTINE cdf_nc_chisq(which,cum,ccum,x,df,pnonc,status, &
          check_input)
! .. Use Statements ..
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind), INTENT (INOUT) :: df, pnonc, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Structures ..
! .. Local Arrays
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

        CALL check_complements(cum,ccum,the_non_central_chi_square%name, &
          'cum','ccum',local_cum,local_ccum,set_values=(which/=1), &
          bad_status=3,status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = x
        params(4) = df
        params(5) = pnonc

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_non_central_chi_square,which, &
            params,status)

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

        CASE (1)
! Solve for CDF

          CALL local_cum_nc_chisq(x,df,pnonc,local_cum,local_ccum,status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)
! Solve for x

          x = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_chi_square,3,local)

          DO
            CALL rc_step_zf(zf_status,x,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_chisq(x,df,pnonc,try_cum,try_ccum,status)

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
! Solve for df

          df = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_chi_square,4,local)

          DO
            CALL rc_step_zf(zf_status,df,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_chisq(x,df,pnonc,try_cum,try_ccum,status)

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
! Solve for pnonc

          pnonc = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_chi_square,5,local)

          DO
            CALL rc_step_zf(zf_status,pnonc,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_chisq(x,df,pnonc,try_cum,try_ccum,status)

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

      END SUBROUTINE cdf_nc_chisq

!*********************************************************************

      FUNCTION cum_nc_chisq(x,df,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_nc_chisq
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: df, pnonc, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_chisq(which=1,cum=cum_nc_chisq,x=x,df=df,pnonc=pnonc, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_nc_chisq

!*********************************************************************

      FUNCTION inv_nc_chisq(cum,ccum,df,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_nc_chisq
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: df, pnonc
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_chisq(2,cum,ccum,inv_nc_chisq,df,pnonc,status, &
          check_input)

        RETURN

      END FUNCTION inv_nc_chisq

!*********************************************************************

      SUBROUTINE local_cum_nc_chisq(x,df,pnonc,cum,ccum,status)
!----------------------------------------------------------------------
!                    Cumulative Non-Central Chi Square distribution

!                              Function

!     Computes the integral from 0 to x for the Non-Central
!     Chi-square distribution.

!                              Arguments

!     x    --> upper limit of integration
!     df   --> degrees of freedom
!     pnonc --> non-centrality parameter of the distribution
!     cum  <-- Cumulative chi-square distribution.
!     ccum <-- Compliment of Cumulative distribution.
!     status <-- status of the computation. 
!                Set to 10 if cdf_chisq has no answer.

!                              Method

!     Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of
!     Mathematical    Functions,  US   NBS   (1966)    to calculate  the
!     non-central chi-square.
!----------------------------------------------------------------------
! .. Use Statements ..
        USE cdf_chisq_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind) :: df, pnonc, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: adj, centaj, centwt, chid2, df0, dfd2, lcntaj, &
          lcntwt, lfact, pcent, pterm, sum, sumadj, term, wt, xnonc
        INTEGER :: chisq_status, i, icent
! ..
! .. Intrinsic Functions ..
        INTRINSIC EXP, INT, LOG, MIN, PRESENT, REAL
! ..
        IF (x<=zero) THEN
          cum = zero
          ccum = one

          RETURN
        END IF

        IF (pnonc<=1.0E-10_dpkind) THEN

!     When non-centrality parameter is (essentially) zero,
!     use cumulative chi-square distribution

          CALL cdf_chisq(which=1,cum=cum,ccum=ccum,x=x,df=df, &
            status=chisq_status)

          IF (chisq_status/=0) THEN
! cdf_chisq has NO answer.

            IF (PRESENT(status)) THEN
              status = 10
            ELSE
              WRITE (*,*) 'Error in local_cum_nc_chisq call to cdf_chisq'
              WRITE (*,*) 'Status: ', chisq_status

              STOP 'Error in local_cum_nc_chisq call to cdf_chisq'
            END IF
          END IF

          RETURN

        END IF

        xnonc = half*pnonc
!***********************************************************************

!     The following code calcualtes the weight, chi-square, and
!     adjustment term for the central term in the infinite series.
!     The central term is the one in which the poisson weight is
!     greatest.  The adjustment term is the amount that must
!     be subtracted from the chi-square to move up two degrees
!     of freedom.

!***********************************************************************
        icent = INT(xnonc)

        IF (icent==0) icent = 1

        chid2 = half*x

!     Calculate central weight term

        lfact = alngam(REAL(icent+1,kind=dpkind))
        lcntwt = (-xnonc) + icent*LOG(xnonc) - lfact
        centwt = EXP(lcntwt)

!     Calculate central chi-square
! DMS dg = df + 2*icent, could be > upper bound for the Chi-squared df

        df0 = dg(icent)

        df0 = MIN(df0,the_chi_square%parameters(4)%high_bound)

        CALL cdf_chisq(which=1,cum=pcent,ccum=ccum,x=x,df=df0, &
          status=chisq_status)

        IF (chisq_status/=0) THEN
! cdf_chisq has NO answer.

          IF (PRESENT(status)) THEN
            status = 10
            RETURN
          ELSE
            WRITE (*,*) 'Error in local_cum_nc_chisq call to cdf_chisq'
            WRITE (*,*) 'Status: ', chisq_status

            STOP 'Error in local_cum_nc_chisq call to cdf_chisq'
          END IF
        END IF

! Calculate central adjustment term

        dfd2 = half*df0
        lfact = alngam(one+dfd2)
        lcntaj = dfd2*LOG(chid2) - chid2 - lfact
        centaj = EXP(lcntaj)
        sum = centwt*pcent

!***********************************************************************

!     Sum backwards from the central term towards zero.
!     Quit whenever either
!     (1) the zero term is reached, or
!     (2) the term gets small relative to the sum, or

!***********************************************************************

        sumadj = zero
        adj = centaj
        wt = centwt
        i = icent

        GO TO 20

10      CONTINUE
        IF (qsmall(term) .OR. i==0) GO TO 30
20      CONTINUE

        dfd2 = half*dg(i)

!     Adjust chi-square for two fewer degrees of freedom.
!     The adjusted value ends up in PTERM.

        adj = adj*dfd2/chid2
        sumadj = sumadj + adj
        pterm = pcent + sumadj

!     Adjust Poisson weight for J decreased by one

        wt = wt*(i/xnonc)
        term = wt*pterm
        sum = sum + term
        i = i - 1
        GO TO 10

30      CONTINUE
        sumadj = centaj
!***********************************************************************

!     Now sum forward from the central term towards infinity.
!     Quit when either
!     (1) the term gets small relative to the sum, or

!***********************************************************************
        adj = centaj
        wt = centwt
        i = icent

        GO TO 50

40      CONTINUE

        IF (qsmall(term)) GO TO 60

!     Update weights for next higher J

50      CONTINUE
        wt = wt*(xnonc/(i+1))

!     Calculate PTERM and add term to sum

        pterm = pcent - sumadj
        term = wt*pterm
        sum = sum + term

!     Update adjustment term for DF for next iteration

        i = i + 1
        dfd2 = half*dg(i)
        adj = adj*chid2/dfd2
        sumadj = sumadj + adj
        GO TO 40

60      CONTINUE
        cum = sum
        ccum = half + (half-cum)

        RETURN

      CONTAINS

!.....................................................................

        FUNCTION dg(i)
! .. Function Return Value ..
          REAL (dpkind) :: dg
! ..
! .. Scalar Arguments ..
          INTEGER, INTENT (IN) :: i
! ..
          dg = df + two*REAL(i,kind=dpkind)

        END FUNCTION dg

!.....................................................................

        FUNCTION qsmall(xx)
! .. Function Return Value ..
          LOGICAL :: qsmall
! ..
! .. Parameters ..
          REAL (dpkind), PARAMETER :: eps = 1.0E-5_dpkind
! ..
! .. Scalar Arguments ..
          REAL (dpkind), INTENT (IN) :: xx
! ..
          qsmall = sum < 1.0E-20_dpkind .OR. xx < eps*sum

        END FUNCTION qsmall

!.....................................................................

      END SUBROUTINE local_cum_nc_chisq

!*********************************************************************

    END MODULE cdf_nc_chisq_mod
