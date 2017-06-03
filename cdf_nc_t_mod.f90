    MODULE cdf_nc_t_mod
! ----------------------------------------------------------------------

!                               cdf_nc_t_mod
!                               *=*=*=*=*=*=

!  -  SUBROUTINE CDF_NC_T(WHICH, CUM, CCUM, T, DF, PNONC, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_NC_T(T,  DF, PNONC, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_NC_T(T, DF, PNONC, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV__NC_T(CUM, CCUM, DF, PNONC,STATUS,CHECK_INPUT)

!                             The Distribution
!                             ================

!   The noncentral T is the distribution of the ratio of two independent
!   random variables. The numerator random variable is distributed as a
!   normal distribution with mean PNONC and variance 1. The denominator
!   random variable is distributed as a the square root of a (central)
!   chi-squared with DF degrees of freedom divided by DF.

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    arguments is to be calculated.
!  Input Range: [ 1:4 ]
!     1.  CUM and CCUM
!     2.  T
!     3.  DF
!     4.  PNONC
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
!  - REAL (dpkind) :: PNONC. The noncentrality parameter.
!  Range: [ 0:10^4 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!      0 Problem solved successfully
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 T outside range
!     -5 DF outside range
!     -6 PNONC outside range
!      3 CUM + CCUM is not nearly one
!     10 can not solve cdf_t or cdf_beta in local_cum_nc_t
!    -50 Answer (if any) is below the lower search bound
!     50 Answer (if any) is above the upper search bound
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
      PUBLIC :: ccum_nc_t, cdf_nc_t, cum_nc_t, inv_nc_t
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_nc_t(t,df,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_nc_t
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: df, pnonc, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_t(which=1,ccum=ccum_nc_t,t=t,df=df,pnonc=pnonc, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_nc_t

!*********************************************************************

      SUBROUTINE cdf_nc_t(which,cum,ccum,t,df,pnonc,status,check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: df, pnonc, t
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

! Check that at least one of cum, ccum is present in the calling list

        CALL check_complements(cum,ccum,the_non_central_t%name,'cum', &
          'ccum',local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
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
        params(5) = pnonc

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_non_central_t,which,params,status)

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
! Solve for cumulative (cdf) value

          CALL local_cum_nc_t(t,df,pnonc,local_cum,local_ccum,status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)

! Solve for t

! Start value for t:
          t = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_t,3,local)

          DO
            CALL rc_step_zf(zf_status,t,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_t(t,df,pnonc,try_cum,try_ccum,status)

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
! Solve for F

! Start value for f:
          df = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_t,4,local)

          DO
            CALL rc_interval_zf(zf_status,df,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_t(t,df,pnonc,try_cum,try_ccum,status)

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
! Start for the non-centrality parameter

! Start value for pnonc:
          pnonc = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_t,5,local)

          DO
            CALL rc_interval_zf(zf_status,pnonc,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_t(t,df,pnonc,try_cum,try_ccum,status)

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

      END SUBROUTINE cdf_nc_t

!*********************************************************************

      FUNCTION cum_nc_t(t,df,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_nc_t
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: df, pnonc, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_t(which=1,cum=cum_nc_t,t=t,df=df,pnonc=pnonc, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_nc_t

!*********************************************************************

      FUNCTION inv_nc_t(cum,ccum,df,pnonc,status,check_input)
! Solves the INVERSE Poisson problem:
! Given cum (or/and ccum), and df, computes t.
! .. Function Return Value ..
        REAL (dpkind) :: inv_nc_t
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind), INTENT (IN) :: df, pnonc
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_t(which=2,cum=cum,ccum=ccum,t=inv_nc_t,df=df, &
          pnonc=pnonc,status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_nc_t

!*********************************************************************

      SUBROUTINE local_cum_nc_t(t,df,pnonc,cum,ccum,status)
!----------------------------------------------------------------------

!                              Function

!     Upper tail    of  the  cumulative  noncentral t   using
!     formulae from page 532  of Johnson, Kotz,  Balakrishnan, Coninuous
!     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
!     This implementation starts the calculation at i = lambda,
!     which is near the largest Di.  It then sums forward and backward.
!     t     --> Upper limit of integration of the t-density.
!     df    --> Degrees of freedom of the t-distribution
!     pnonc --> the non-centrality parameter
!     cum   <-- Cumulative t-distribution.
!     ccum  <-- Compliment of Cumulative t-distribution.
!     status <-- status of the computation. 
!                Set to 10 if cdf_t or cdf_beta has no answer.

!                         Method

!     Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of
!     Mathematical Functions  (1966) is used to reduce the computation
!     of the cumulative distribution function to that of an incomplete
!     beta.
!------------------------------------------------------------------
! .. Use Statements ..
        USE biomath_mathlib_mod
        USE cdf_beta_mod
        USE cdf_normal_mod
        USE cdf_t_mod
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: conv = 1.0E-7_dpkind
        REAL (dpkind), PARAMETER :: onep5 = 1.5E0_dpkind
        REAL (dpkind), PARAMETER :: tiny = 1.0E-10_dpkind
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: df, pnonc, t
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: alghdf, b, bb, bbcent, bcent, cent, d, dcent, &
          dpnonc, dum1, dum2, e, ecent, halfdf, lambda, lnomx, lnx, omx, &
          pnonc2, s, scent, ss, sscent, t2, term, tt, twoi, x, xi, xlnd, &
          xlne
        INTEGER :: local_status
        LOGICAL :: has_status, qrevs
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, INT, LOG, MAX, MIN, PRESENT
! ..
        has_status = PRESENT(status)

! Case pnonc essentially zero

        IF (ABS(pnonc)<=tiny) THEN
          CALL cdf_t(1,cum,ccum,t,df,status=local_status)

          IF (local_status/=0) THEN
! cdf_t has NO answer.

            IF (has_status) THEN
              status = 10
            ELSE
              WRITE (*,*) 'Error in local_cum_nc_t call to cdf_t'
              WRITE (*,*) 'Status: ', local_status

              STOP 'Error in local_cum_nc_t call to cdf_t'

            END IF
          END IF

          RETURN
        END IF

        qrevs = (t<zero)
        IF (qrevs) THEN
          tt = -t
          dpnonc = -pnonc
        ELSE
          tt = t
          dpnonc = pnonc
        END IF

        pnonc2 = dpnonc*dpnonc
        t2 = tt*tt

        IF (ABS(tt)<=tiny) THEN
! NORMAL has NO errors (other than bounds -- but we checked)

          CALL cdf_normal(which=1,cum=cum,ccum=ccum,x=-pnonc)

          RETURN
        END IF

        lambda = half*pnonc2
        x = df/(df+t2)
        omx = one - x

        lnx = LOG(x)
        lnomx = LOG(omx)

        halfdf = half*df
        alghdf = gamln(halfdf)

! ***** Case i = lambda *****

        cent = INT(lambda)

        cent = MAX(one,cent)

! Compute d=T(2i) in log space and offset by exp(-lambda)

        xlnd = cent*LOG(lambda) - gamln(cent+one) - lambda

        dcent = EXP(xlnd)

! Compute e=t(2i+1) in log space offset by exp(-lambda)

        xlne = (cent+half)*LOG(lambda) - gamln(cent+onep5) - lambda
        ecent = EXP(xlne)

        IF (dpnonc<zero) ecent = -ecent

! Compute bcent=B(2*cent)

        CALL cdf_beta(which=1,cum=bcent,ccum=dum1,x=x,a=halfdf, &
          b=cent+half,check_input=.FALSE.,status=local_status)

        IF (local_status/=0) THEN
! cdf_beta has NO answer.

          IF (PRESENT(status)) THEN
            status = 10

            RETURN
          ELSE
            WRITE (*,*) 'Error in local_cum_nc_t call to cdf_beta'
            WRITE (*,*) 'Status: ', local_status

            STOP 'Error in local_cum_nc_t call to cdf_beta'

          END IF
        END IF

! Compute bbcent=B(2*cent+1)

        CALL cdf_beta(which=1,cum=bbcent,ccum=dum2,x=x,a=halfdf, &
          b=cent+one,check_input=.FALSE.,status=local_status)

        IF (local_status/=0) THEN
! cdf_beta has NO answer.

          IF (PRESENT(status)) THEN
            status = 10

            RETURN
          ELSE
            WRITE (*,*) 'Error in local_cum_nc_t call to cdf_beta'
            WRITE (*,*) 'Status: ', local_status

            STOP 'Error in local_cum_nc_t call to cdf_beta'

          END IF
        END IF

! Case bcent and bbcent are essentially zero
! Thus t is effectively infinite

        IF (bcent+bbcent<tiny) THEN
          IF (qrevs) THEN
            cum = zero
            ccum = one
          ELSE
            cum = one
            ccum = zero
          END IF

          RETURN

        END IF

! Case bcent and bbcent are essentially one
! Thus t is effectively zero

        IF (dum1+dum2<tiny) THEN
          CALL cdf_normal(which=1,cum=cum,ccum=ccum,x=-pnonc)

          RETURN
        END IF

! First term in ccum is D*B + E*BB

        ccum = dcent*bcent + ecent*bbcent

! Compute: s(cent) = B(2*(cent+1)) - B(2*cent))

        scent = gamln(halfdf+cent+half) - gamln(cent+onep5) - alghdf + &
          halfdf*lnx + (cent+half)*lnomx
        scent = EXP(scent)

! Compute: ss(cent) = B(2*cent+3) - B(2*cent+1)

        sscent = gamln(halfdf+cent+one) - gamln(cent+two) - alghdf + &
          halfdf*lnx + (cent+one)*lnomx
        sscent = EXP(sscent)

! ***** Sum FORWARD *****

        xi = cent + one
        twoi = two*xi

        d = dcent
        e = ecent
        b = bcent
        bb = bbcent
        s = scent
        ss = sscent

        DO
          b = b + s
          bb = bb + ss

          d = (lambda/xi)*d
          e = (lambda/(xi+half))*e

          term = d*b + e*bb
          ccum = ccum + term

          s = s*omx*(df+twoi-one)/(twoi+one)
          ss = ss*omx*(df+twoi)/(twoi+two)

          xi = xi + one
          twoi = two*xi

          IF (ABS(term)<=conv*ccum) THEN
            EXIT
          END IF

        END DO

! ***** Sum BACKWARD *****

        xi = cent
        twoi = two*xi

        d = dcent
        e = ecent
        b = bcent
        bb = bbcent

        s = scent*(one+twoi)/((df+twoi-one)*omx)
        ss = sscent*(two+twoi)/((df+twoi)*omx)

        DO
          b = b - s
          bb = bb - ss

          d = d*(xi/lambda)
          e = e*((xi+half)/lambda)

          term = d*b + e*bb
          ccum = ccum + term

          xi = xi - one

          IF (xi<half) THEN
            EXIT
          END IF

          twoi = two*xi

          s = s*(one+twoi)/((df+twoi-one)*omx)
          ss = ss*(two+twoi)/((df+twoi)*omx)

          IF (ABS(term)<=conv*ccum) THEN
            EXIT
          END IF

        END DO

        IF (qrevs) THEN
          cum = half*ccum
          ccum = one - cum
        ELSE
          ccum = half*ccum
          cum = one - ccum
        END IF

! Due to roundoff error the answer may not lie between zero and one
! Force it to do so

        cum = MAX(MIN(cum,one),zero)
        ccum = MAX(MIN(ccum,one),zero)

      END SUBROUTINE local_cum_nc_t

!*********************************************************************

    END MODULE cdf_nc_t_mod
