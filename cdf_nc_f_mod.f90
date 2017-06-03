    MODULE cdf_nc_f_mod
! ----------------------------------------------------------------------

!                               cdf_nc_f_mod
!                               *=*=*=*=*=*=

!  -  SUBROUTINE CDF_F(WHICH, CUM, CCUM, F, DFN, DFD, PNONC, STATUS,
!                      CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_F(F,  DFN, DFD, PNONC, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_F(F, DFN, DFD, PNONC, STATUS,  CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_F(CUM, CCUM, DFN, DFD, PNONC, STATUS,
!                                  CHECK_INPUT)

!                             The Distribution
!                             ================

!    The noncentral F is the distribution of the ratio of two independent
!    random variables. The numerator random variable is distributed as a
!    noncentral chi-squared with DFN degrees of freedom and noncentrality
!    parameter PNONC divided by DFN. The denominator random variable is
!    distributed as a (central) chi-squared with DFD degrees of freedom
!    divided by DFD.

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    five arguments is to be calculated.
!  Input Range: [ 1:2 ]
!     1.  CUM and CCUM
!     2.  F
!  NOTE: DFN and DFD will not be computed because CUM is not monotone in
!    either argument.
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the noncentral f
!    distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the noncentral
!    f distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: F. The upper limit of integration of the noncentral
!    f density. The lower limit is 0.
!  Range: [ 0:10^100 ]
!  - REAL (dpkind) :: DFN. The numerator degrees of freedom.
!  Range: [ 10^-3:10^10 ]
!  - REAL (dpkind) :: DFD. The denominator degrees of freedom.
!  Range: [ 10^-3:10^10 ]
!  - REAL (dpkind) :: PNONC. The noncentrality parameter.
!  Range: [ 0:10^4 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!      0 Problem solved successfully
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 F outside range
!     -5 DFN outside range
!     -6 DFD outside range
!      3 CUM + CCUM is not nearly one
!     10 can not solve cdf_f in local_cum_nc_f
!    -50 Answer (if any) is BELOW the LOWER search bound
!     50 Answer (if any) is ABOVE the UPPER search bound
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.

!  NOTE: CUM and CCUM must add to (nearly) one.
!----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
      USE biomath_mathlib_mod
      USE cdf_aux_mod
      USE zero_finder
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_nc_f, cdf_nc_f, cum_nc_f, inv_nc_f
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_nc_f(f,dfn,dfd,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_nc_f
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: dfd, dfn, f, pnonc
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_f(which=1,ccum=ccum_nc_f,f=f,dfn=dfn,dfd=dfd, &
          pnonc=pnonc,status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_nc_f

!*********************************************************************

      SUBROUTINE cdf_nc_f(which,cum,ccum,f,dfn,dfd,pnonc,status, &
          check_input)
! .. Use Statements ..
        USE cdf_gamma_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: dfd, dfn, f, pnonc
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
        REAL (dpkind) :: fx, local_ccum, local_cum, try_ccum, try_cum
        INTEGER :: local_status, zf_status
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

        CALL check_complements(cum,ccum,the_non_central_f%name,'cum', &
          'ccum',local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
          status=status)

        IF (has_status) THEN
          IF (status/=0) THEN
            RETURN
          END IF
        END IF

        params(1) = local_cum
        params(2) = local_ccum
        params(3) = f
        params(4) = dfn
        params(5) = dfd
        params(6) = pnonc

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_non_central_f,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF

        END IF

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute the Answers
!++++++++++++++++++++++++++++++++++++++++++++++++++

        local_status = 0

        SELECT CASE (which)

        CASE (1)
! Calculate cum and ccum

          CALL local_cum_nc_f(f,dfn,dfd,pnonc,local_cum,local_ccum, &
            status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)
! Calculate f

! Start value for f
          f = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_f,3,local)

          DO
            CALL rc_interval_zf(zf_status,f,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_nc_f(f,dfn,dfd,pnonc,try_cum,try_ccum,status)

            IF (has_status) THEN
              IF (status/=0) THEN
                RETURN
              END IF
            END IF

            fx = try_cum - local_cum
          END DO

! Can NOT solve for Degrees of Freedom (which=4 and 5)

        CASE (3)
! Solve for PNONC

! Start value for pnonc:
          pnonc = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_non_central_f,6,local)

          DO
            CALL rc_interval_zf(zf_status,pnonc,fx,local)


            IF (zf_status/=1) EXIT

            CALL local_cum_nc_f(f,dfn,dfd,pnonc,try_cum,try_ccum,status)

            IF (has_status) THEN
              IF (status/=0) THEN
                RETURN
              END IF
            END IF

            fx = try_cum - local_cum
          END DO

        END SELECT

        IF (has_status) THEN
! Set the status of the zero finder
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_nc_f

!*********************************************************************

      FUNCTION cum_nc_f(f,dfn,dfd,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_nc_f
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: dfd, dfn, f, pnonc
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_f(which=1,cum=cum_nc_f,f=f,dfn=dfn,dfd=dfd, &
          pnonc=pnonc,status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_nc_f

!*********************************************************************

      FUNCTION inv_nc_f(cum,ccum,dfn,dfd,pnonc,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_nc_f
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: dfd, dfn, pnonc
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_nc_f(which=2,cum=cum,ccum=ccum,f=inv_nc_f,dfn=dfn, &
          dfd=dfd,pnonc=pnonc,status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_nc_f

!*********************************************************************

      SUBROUTINE local_cum_nc_f(f,dfn,dfd,pnonc,cum,ccum,status)
!**********************************************************************
!               F -NON- -C-ENTRAL F DISTRIBUTION

!                              Function

!     Computes Non-Central F Distribution with DFN and DFD
!     Degrees of Freedom and Non-Centrality parameter PNONC

!                              Arguments

!     F --> upper limit of integration of Non-Central F
!     DFN --> degrees of freedom of numerator
!     DFD --> degrees of freedom of denominator
!     PNONC --> non-centrality parameter.
!     CUM <-- cumulative Non-Central F distribution
!     CCUM <-- compliment of cummulative (1-CUM)

!                              Method

!     Uses formula 26.6.20 of reference for infinite series.
!     Series is calculated backward and forward from j = lambda/2
!     (this is the term with the largest Poisson weight) until
!     the convergence criterion is met.
!     for speed, the incomplete beta functions are evaluated
!     by formula 26.5.16.

!               Reference

!     Handbood of Mathematical Functions
!     by Milton Abramowitz and Irene A. Stegun
!     NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55
!     March 1965
!     P 947, equations 26.6.17, 26.6.18

!                              Note

!     The sum continues until a succeeding term is less than EPS
!     times the sum (or the sum is less than 1.0e-20).  EPS is
!     set to 1.0E-4 in a data statement which can be changed.
!**********************************************************************
! .. Use Statements ..
        USE cdf_f_mod
        USE cdf_beta_mod
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: eps = 1.0E-4_dpkind
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: dfd, dfn, f, pnonc
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: adn, aup, b, betdn, betup, centwt, dnterm, dsum, &
          prod, sum, upterm, xmult, xnonc, xx, yy
        INTEGER :: f_status, i, icent
! ..
! .. Intrinsic Functions ..
        INTRINSIC EXP, LOG, PRESENT, REAL
! ..
        IF (f<=zero) THEN
          cum = zero
          ccum = one

          RETURN
        END IF

        IF (pnonc<1.0E-10_dpkind) THEN

!     Handle case in which the non-centrality parameter is
!     (essentially) zero.

          CALL cdf_f(which=1,cum=cum,ccum=ccum,f=f,dfn=dfn,dfd=dfd, &
            status=f_status)

          IF (f_status/=0) THEN
! cdf_f has NO answer.

            IF (PRESENT(status)) THEN
              status = 10
            ELSE
              WRITE (*,*) 'Error in local_cum_nc_f call to cdf_f'
              WRITE (*,*) 'Status: ', f_status

              STOP 'Error in local_cum_nc_f call to cdf_f'
            END IF
          END IF

          RETURN
        END IF

        xnonc = half*pnonc

! Calculate the central term of the poisson weighting factor.

        icent = xnonc
        IF (icent==0) icent = 1

! Compute central weight term

        centwt = EXP((-xnonc)+icent*LOG(xnonc)-alngam(REAL(icent+1, &
          kind=dpkind)))

!     Compute central incomplete beta term
!     Assure that minimum of arg to beta and 1 - arg is computed accurately.

        prod = dfn*f
        dsum = dfd + prod
        yy = dfd/dsum

        IF (yy>half) THEN
          xx = prod/dsum
          yy = one - xx
        ELSE
          xx = one - yy
        END IF

        betdn = cum_beta(x=xx,a=dfn*half+REAL(icent,kind=dpkind), &
          b=dfd*half,check_input=.FALSE.)

        adn = half*dfn + REAL(icent,kind=dpkind)
        aup = adn
        b = half*dfd
        betup = betdn
        sum = centwt*betdn

!     Now sum terms backward from icent until convergence or all done

        xmult = centwt
        i = icent
        dnterm = EXP(alngam(adn+b)-alngam(adn+one)-alngam(b)+adn*LOG(xx)+ &
          b*LOG(yy))

10      CONTINUE
        IF (qsmall(xmult*betdn) .OR. i<=0) GO TO 20
        xmult = xmult*(i/xnonc)
        i = i - 1
        adn = adn - 1
        dnterm = (adn+1)/((adn+b)*xx)*dnterm
        betdn = betdn + dnterm
        sum = sum + xmult*betdn
        GO TO 10

20      CONTINUE
        i = icent + 1

!     Now sum forwards until convergence

        xmult = centwt
        IF (aup-one+b==0) THEN
          upterm = EXP((-alngam(aup))-alngam(b)+(aup-one)*LOG(xx)+b*LOG( &
            yy))
        ELSE
          upterm = EXP(alngam(aup-one+b)-alngam(aup)-alngam(b)+ &
            (aup-one)*LOG(xx)+b*LOG(yy))
        END IF

        GO TO 40

30      CONTINUE
        IF (qsmall(xmult*betup)) GO TO 50

40      CONTINUE
        xmult = xmult*(xnonc/i)
        i = i + 1
        aup = aup + 1
        upterm = (aup+b-2.0E0_dpkind)*xx/(aup-1)*upterm
        betup = betup - upterm
        sum = sum + xmult*betup
        GO TO 30

50      CONTINUE
        cum = sum

        ccum = half + (half-cum)
        RETURN

      CONTAINS

!.....................................................................

        FUNCTION qsmall(x)
! .. Function Return Value ..
          LOGICAL :: qsmall
! ..
! .. Scalar Arguments ..
          REAL (dpkind) :: x
! ..
          qsmall = sum < 1.0E-20_dpkind .OR. x < eps*sum

          RETURN

        END FUNCTION qsmall

!.....................................................................

      END SUBROUTINE local_cum_nc_f

!*********************************************************************

    END MODULE cdf_nc_f_mod
