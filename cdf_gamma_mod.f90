    MODULE cdf_gamma_mod
! ----------------------------------------------------------------------

!                              cdf_gamma_mod
!                              *=*=*=*=*=*=*

!  -  SUBROUTINE CDF_GAMMA(WHICH, CUM, CCUM, X, SHAPE, SCALE, STATUS,
!                          CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_GAMMA(X, SHAPE, SCALE, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_GAMMA(X, SHAPE, SCALE, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_GAMMA(CUM, CCUM, SHAPE, SCALE, STATUS,
!                                      CHECK_INPUT)

!                             The Distribution
!                             ================

!         The density of the GAMMA distribution is proportional to:

!                               SHAPE-1
!                      (x/SCALE)        exp(-x/SCALE)

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    four arguments is to be calculated.
!  Input Range: [ 1:4 ]
!     1.  CUM and CCUM
!     2.  X
!     3.  SHAPE
!     4.  SCALE
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the gamma distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the gamma
!                               distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: X. The upper limit of integration of the gamma
!                        density. The lower limit is 0.
!  Range: [ 0:10^100 ]
!  - REAL (dpkind) :: SHAPE. The shape parameter of the distribution.
!  Range: [ 10^-10:10^100 ]
!  - REAL (dpkind) :: SCALE. The scale parameter of the distribution.
!  Range: [ 10^-10:10^100 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. Possible values:
!      0 Problem solved successfully
!     -1 WHICH outside input range
!     -2 CUM outside range
!     -3 CCUM outside range
!     -4 X outside range
!     -5 SHAPE outside range
!     -6 SCALE outside range
!      3 CUM + CCUM is not nearly one
!     10 can not solve gamma_inverse in cdf_gamma
!    -50 Answer (if any) is BELOW the LOWER search bound
!     50 Answer (if any) is ABOVE the UPPER search bound
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.

!  NOTE: CUM and CCUM must add to (nearly) one.
! ----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
      USE biomath_mathlib_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_gamma, cdf_gamma, cum_gamma, inv_gamma
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_gamma(x,shape,scale,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_gamma
! ..
! .. Scalar Arguments ..
! DMS OPTIONAL added
        REAL (dpkind), INTENT (IN) :: scale, shape, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_gamma(which=1,ccum=ccum_gamma,x=x,shape=shape, &
          scale=scale,status=status,check_input=check_input)

        RETURN

      END FUNCTION ccum_gamma

!*********************************************************************

      SUBROUTINE cdf_gamma(which,cum,ccum,x,shape,scale,status, &
          check_input)
! .. Use Statements ..
        USE zero_finder
        USE cdf_aux_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: scale, shape, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER, INTENT (IN) :: which
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
! .. Local Scalars ..
        REAL (dpkind) :: fx, local_ccum, local_cum, try_ccum, try_cum, xx
        INTEGER :: ierr, zf_status
        LOGICAL :: has_status, local_check_input, match_cum
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! .. Local Structures ..
        TYPE (zf_locals) :: local
! ..
! .. Local Arrays ..
        REAL (dpkind) :: params(6)
! ..
        has_status = PRESENT(status)

! status = 0 means NO error

        IF (has_status) THEN
          status = 0
        END IF

        CALL check_complements(cum,ccum,the_gamma%name,'cum','ccum', &
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
        params(4) = shape
        params(5) = scale

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_gamma,which,params,status)

          IF (has_status) THEN
            IF (status/=0) THEN
              RETURN
            END IF
          END IF
        END IF

        IF (which>1) match_cum = (local_cum<=half)

        SELECT CASE (which)

        CASE (1)
! Calculate cum and ccum

          CALL local_cum_gamma(x,shape,scale,local_cum,local_ccum)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

          RETURN

        CASE (2)
! Calculate x

          CALL gamma_inverse(shape,xx,-one,local_cum,local_ccum,ierr)

          IF (ierr<0) THEN
            IF (has_status) THEN
              status = 10

              RETURN
            ELSE
              WRITE (*,*) 'Error in cdf_gamma call to gamma_inverse'
              WRITE (*,*) 'Status: ', ierr

              STOP 'Error in cdf_gamma call to gamma_inverse'
            END IF

          END IF

          x = xx/scale

          RETURN

        CASE (3)
! Calculate SHAPE

          shape = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_gamma,4,local)

          DO

            CALL rc_step_zf(zf_status,shape,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_gamma(x,shape,scale,try_cum,try_ccum)

            IF (match_cum) THEN
              fx = try_cum - local_cum
            ELSE
              fx = try_ccum - local_ccum
            END IF

          END DO

        CASE (4)
! Calculate SCALE

          CALL gamma_inverse(shape,xx,-one,local_cum,local_ccum,ierr)

          IF (ierr/=0) THEN
            IF (has_status) THEN
              status = 10

              RETURN
            ELSE
              WRITE (*,*) 'Error in cdf_gamma call to gamma_inverse'
              WRITE (*,*) 'Status: ', ierr

              STOP 'Error in cdf_gamma call to gamma_inverse'
            END IF

          END IF

          scale = xx/x

          RETURN

        END SELECT

! Gets here ONLY for which=3

        IF (has_status) THEN
! Set the status of the zero finder
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_gamma

!*********************************************************************

      FUNCTION cum_gamma(x,shape,scale,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_gamma
! ..
! .. Scalar Arguments ..
! DMS OPTIONAL added
        REAL (dpkind), INTENT (IN) :: scale, shape, x
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_gamma(which=1,cum=cum_gamma,x=x,shape=shape,scale=scale, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION cum_gamma

!*********************************************************************

      SUBROUTINE gamma_inverse(a,x,x0,p,q,ierr)
! ----------------------------------------------------------------------
!            INVERSE INCOMPLETE GAMMA RATIO FUNCTION
!     GIVEN POSITIVE A, AND NON-NEGATIVE P AND Q WHERE P + Q = 1.
!     THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
!     ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
!     TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
!     PARTICULAR COMPUTER ARITHMETIC BEING USED.
!                      ------------
!     X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
!     AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
!     NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
!     A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
!     IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.
!     X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER
!     DOES NOT WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET
!     X0 .LE. 0.
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!     WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
!     VALUES ...
!       IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
!                    NOT USED.
!       IERR.GT.0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
!                    WERE PERFORMED.
!       IERR = -2    (INPUT ERROR) A .LE. 0
!       IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
!                    IS TOO LARGE.
!       IERR = -4    (INPUT ERROR) P + Q .NE. 1
!       IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
!                    RECENT VALUE OBTAINED FOR X IS GIVEN.
!                    THIS CANNOT OCCUR IF X0 .LE. 0.
!       IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
!                    THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
!       IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
!                    ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
!                    ITERATION CANNOT BE PERFORMED IN THIS
!                    CASE. IF X0 .LE. 0, THIS CAN OCCUR ONLY
!                    WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS
!                    POSITIVE THEN THIS CAN OCCUR WHEN A IS
!                    EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
!                    LARGE (SAY A .GE. 1.E20).
! ----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!     -------------------
!     -------------------
!     LN10 = LN(10)
!     C = EULER CONSTANT
!     -------------------
!     -------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, p, q, x0
        REAL (dpkind), INTENT (OUT) :: x
        INTEGER, INTENT (OUT) :: ierr
! ..
! .. Local Scalars ..
        REAL (dpkind) :: am1, amax, ap1, ap2, ap3, apn, b, c1, c2, c3, &
          c4, c5, d, e, e2, eps, g, h, pn, qg, qn, r, rta, s, s2, sum, t, &
          u, w, xmax, xmin, xn, y, z
        INTEGER :: iop
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EPSILON, EXP, HUGE, LOG, MAX, SQRT, TINY
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: ac(4) = (/ 3.31125922108741E0_dpkind, &
          11.6616720288968E0_dpkind, 4.28342155967104E0_dpkind, &
          .213623493715853E0_dpkind/)
        REAL (dpkind), PARAMETER :: amin(2) = (/ 500.0E0_dpkind, &
          100.0E0_dpkind/)
        REAL (dpkind), PARAMETER :: bc(5) = (/ one, &
          6.61053765625462E0_dpkind, 6.40691597760039E0_dpkind, &
          1.27364489782223E0_dpkind, .036117081018842E0_dpkind/)
        REAL (dpkind), PARAMETER :: bmin(2) = (/ 1.E-28_dpkind, &
          1.E-13_dpkind/)
        REAL (dpkind), PARAMETER :: dmin(2) = (/ 1.E-06_dpkind, &
          1.E-04_dpkind/)
        REAL (dpkind), PARAMETER :: emin(2) = (/ 2.E-03_dpkind, &
          6.E-03_dpkind/)
        REAL (dpkind), PARAMETER :: eps0(2) = (/ 1.E-10_dpkind, &
          1.E-08_dpkind/)
        REAL (dpkind), PARAMETER :: c = .577215664901533E0_dpkind
        REAL (dpkind), PARAMETER :: ln10 = 2.302585E0_dpkind
        REAL (dpkind), PARAMETER :: tol = 1.E-5_dpkind
! ..
!     -------------------
!     ****** E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
!            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
!            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
!            LARGEST POSITIVE NUMBER.
        e = EPSILON(one)
        xmin = TINY(one)
        xmax = HUGE(one)

        x = zero
        IF (a<=zero) GO TO 240
        t = p + q - one
        IF (ABS(t)>e) GO TO 260

        ierr = 0
        IF (p==zero) RETURN
        IF (q==zero) GO TO 220
        IF (a==one) GO TO 230

        e2 = two*e
        amax = 0.4E-10_dpkind/(e*e)
        iop = 1

        IF (e>1.E-10_dpkind) iop = 2

        eps = eps0(iop)
        xn = x0

        IF (x0>zero) GO TO 110

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A .LT. 1

        IF (a>one) GO TO 50

        g = gamma(a+one)
        qg = q*g
        IF (qg==zero) GO TO 300

        b = qg/a

        IF (qg>0.6E0_dpkind*a) GO TO 40
        IF (a>=0.30E0_dpkind .OR. b<0.35E0_dpkind) GO TO 10

        t = EXP(-(b+c))
        u = t*EXP(t)
        xn = t*EXP(u)
        GO TO 110

10      CONTINUE
        IF (b>=0.45E0_dpkind) GO TO 40
        IF (b==zero) GO TO 300
        y = -LOG(b)
        s = half + (half-a)
        z = LOG(y)
        t = y - s*z
        IF (b<0.15E0_dpkind) GO TO 20
        xn = y - s*LOG(t) - LOG(one+s/(t+one))
        GO TO 170

20      CONTINUE
        IF (b<=0.01E0_dpkind) GO TO 30
        u = ((t+two*(three-a))*t+(two-a)*(three-a))/((t+(five-a))*t+two)
        xn = y - s*LOG(t) - LOG(u)
        GO TO 170

30      CONTINUE
        c1 = -s*z
        c2 = -s*(one+c1)
        c3 = s*((half*c1+(two-a))*c1+(2.5E0_dpkind-1.5E0_dpkind*a))
        c4 = -s*(((c1/three+(2.5E0_dpkind-1.5E0_dpkind*a))*c1+((a-six)*a+ &
          seven))*c1+((11.0E0_dpkind*a-46)*a+47.0E0_dpkind)/six)
        c5 = -s*(((((-c1/four)+(11.0E0_dpkind*a-17.0E0_dpkind)/six)*c1+ &
          (((-three*a)+13.0E0_dpkind)*a-13.0E0_dpkind))*c1+half*(((two*a- &
          25.0E0_dpkind)*a+72.0E0_dpkind)*a-61.0E0_dpkind))*c1+((( &
          25.0E0_dpkind*a-195.0E0_dpkind)*a+477.0E0_dpkind)*a- &
          379.0E0_dpkind)/twelve)

        xn = ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y

        IF (a>one) GO TO 170
        IF (b>bmin(iop)) GO TO 170
        x = xn
        RETURN

40      CONTINUE
        IF (b*q<=1.E-8_dpkind) THEN
          xn = EXP(-(q/a+c))
        ELSE

          IF (p>0.9E0_dpkind) THEN
            xn = EXP((alnrel(-q)+gamln1(a))/a)
          ELSE

            xn = EXP(LOG(p*g)/a)
          END IF
        END IF

        IF (xn==zero) GO TO 250

        t = half + (half-xn/(a+one))
        xn = xn/t
        GO TO 110

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A .GT. 1

50      CONTINUE
        IF (q>half) THEN
          w = LOG(p)
        ELSE
          w = LOG(q)
        END IF

        t = SQRT(-two*w)

        s = t - evaluate_polynomial(ac,t)/evaluate_polynomial(bc,t)

        IF (q>half) s = -s

        rta = SQRT(a)
        s2 = s*s
        xn = a + s*rta + (s2-one)/three + s*(s2-seven)/(36.0E0_dpkind*rta &
          ) - ((three*s2+seven)*s2-16.0E0_dpkind)/(810.0E0_dpkind*a) + &
          s*((nine*s2+256.0E0_dpkind)*s2-433.0E0_dpkind)/ &
          (38880.0E0_dpkind*a*rta)
        xn = MAX(xn,zero)

        IF (a<amin(iop)) GO TO 60

        x = xn
        d = half + (half-x/a)

        IF (ABS(d)<=dmin(iop)) RETURN

60      CONTINUE
        IF (p<=half) GO TO 80
        IF (xn<three*a) GO TO 170

        y = -(w+gamln(a))
        d = MAX(two,a*(a-one))

        IF (y<ln10*d) GO TO 70

        s = one - a
        z = LOG(y)
        GO TO 30

70      CONTINUE
        t = a - one
        xn = y + t*LOG(xn) - alnrel(-t/(xn+one))
        xn = y + t*LOG(xn) - alnrel(-t/(xn+one))
        GO TO 170

80      CONTINUE
        ap1 = a + one

        IF (xn>0.70E0_dpkind*ap1) GO TO 120

        w = w + gamln(ap1)

        IF (xn>0.15E0_dpkind*ap1) GO TO 90

        ap2 = a + two
        ap3 = a + three
        x = EXP((w+x)/a)
        x = EXP((w+x-LOG(one+(x/ap1)*(one+x/ap2)))/a)
        x = EXP((w+x-LOG(one+(x/ap1)*(one+x/ap2)))/a)
        x = EXP((w+x-LOG(one+(x/ap1)*(one+(x/ap2)*(one+x/ap3))))/a)
        xn = x

        IF (xn>1.E-2_dpkind*ap1) GO TO 90
        IF (xn<=emin(iop)*ap1) RETURN
        GO TO 120

90      CONTINUE
        apn = ap1
        t = xn/apn
        sum = one + t
100     CONTINUE
        apn = apn + one
        t = t*(xn/apn)
        sum = sum + t
        IF (t>1.E-4_dpkind) GO TO 100
        t = w - LOG(sum)
        xn = EXP((xn+t)/a)
        xn = xn*(one-(a*LOG(xn)-xn-t)/(a-xn))
        GO TO 120

!                 SCHRODER ITERATION USING P

110     CONTINUE
        IF (p>half) GO TO 170

120     CONTINUE
        IF (p<=1.E10_dpkind*xmin) GO TO 290
        am1 = (a-half) - half

130     CONTINUE
        IF (a<=amax) GO TO 140
        d = half + (half-xn/a)
        IF (ABS(d)<=e2) GO TO 290

140     CONTINUE
        IF (ierr>=20) GO TO 270
        ierr = ierr + 1

        CALL gratio(a,xn,pn,qn,0)

        IF (pn==zero .OR. qn==zero) GO TO 290

        r = rcomp(a,xn)

        IF (r==zero) GO TO 290

        t = (pn-p)/r
        w = half*(am1-xn)

        IF (ABS(t)<=tenth .AND. ABS(w*t)<=tenth) GO TO 150

        x = xn*(one-t)

        IF (x<=zero) GO TO 280

        d = ABS(t)
        GO TO 160

150     CONTINUE
        h = t*(one+w*t)
        x = xn*(one-h)
        IF (x<=zero) GO TO 280
        IF (ABS(w)>=one .AND. ABS(w)*t*t<=eps) RETURN
        d = ABS(h)

160     CONTINUE
        xn = x
        IF (d>tol) GO TO 130
        IF (d<=eps) RETURN
        IF (ABS(p-pn)<=tol*p) RETURN
        GO TO 130

!                 SCHRODER ITERATION USING Q

170     CONTINUE
        IF (q<=1.E10_dpkind*xmin) GO TO 290
        am1 = (a-half) - half

180     CONTINUE
        IF (a<=amax) GO TO 190
        d = half + (half-xn/a)
        IF (ABS(d)<=e2) GO TO 290

190     CONTINUE
        IF (ierr>=20) GO TO 270
        ierr = ierr + 1

        CALL gratio(a,xn,pn,qn,0)

        IF (pn==zero .OR. qn==zero) GO TO 290

        r = rcomp(a,xn)
        IF (r==zero) GO TO 290
        t = (q-qn)/r
        w = half*(am1-xn)
        IF (ABS(t)<=tenth .AND. ABS(w*t)<=tenth) GO TO 200
        x = xn*(one-t)
        IF (x<=zero) GO TO 280
        d = ABS(t)
        GO TO 210

200     CONTINUE
        h = t*(one+w*t)
        x = xn*(one-h)
        IF (x<=zero) GO TO 280
        IF (ABS(w)>=one .AND. ABS(w)*t*t<=eps) RETURN
        d = ABS(h)

210     CONTINUE
        xn = x
        IF (d>tol) GO TO 180
        IF (d<=eps) RETURN
        IF (ABS(q-qn)<=tol*q) RETURN
        GO TO 180

!                       SPECIAL CASES

220     CONTINUE
        x = xmax
        RETURN

230     CONTINUE
        IF (q>=0.9E0_dpkind) THEN
          x = -alnrel(-p)
          RETURN

        END IF

        x = -LOG(q)
        RETURN

!                       ERROR RETURN

240     CONTINUE
        ierr = -2
        RETURN

250     CONTINUE
        ierr = -3
        RETURN

260     CONTINUE
        ierr = -4
        RETURN

270     CONTINUE
        ierr = -6
        RETURN

280     CONTINUE
        ierr = -7
        RETURN

290     CONTINUE
        x = xn
        ierr = -8
        RETURN

300     CONTINUE
        x = xmax
        ierr = -8

      END SUBROUTINE gamma_inverse

!*********************************************************************

      FUNCTION inv_gamma(cum,ccum,shape,scale,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: inv_gamma
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: scale, shape
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_gamma(2,cum=cum,ccum=ccum,x=inv_gamma,shape=shape, &
          scale=scale,status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_gamma

!*********************************************************************

      SUBROUTINE local_cum_gamma(x,shape,scale,cum,ccum)
!----------------------------------------------------------------------

!                              Function

!     Computes   the  cumulative        of    the     incomplete   gamma
!     distribution, i.e., the integral from 0 to X*SCALE of
!          (1/GAM(SHAPE))*EXP(-T)*T**(SHAPE-1) DT
!     where GAM(SHAPE) is the complete gamma function of SHAPE, i.e.,
!          GAM(SHAPE) = integral from 0 to infinity of

!                    EXP(-T)*T**(SHAPE-1) DT

!                              Arguments

!     X --> The upper limit of integration of the incomplete gamma.
!     SHAPE --> The shape parameter of the incomplete gamma.
!     SCALE --> The scale parameter of the incomplete gamma
!     CUM <-- Cumulative incomplete gamma distribution.
!     CCUM <-- Compliment of Cumulative incomplete gamma distribution.

!                              Method

!    The  cumulative  gamma and  its  inverse  are  from the  following
!    reference:
!    DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
!    gamma function  ratios  and their  inverse.   ACM  Trans.  Math.
!    Softw. 12 (1986), 377-393.
!----------------------------------------------------------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: scale, shape, x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: arg
! ..
        IF (x<=zero) THEN
          cum = zero
          ccum = one

          RETURN
        END IF

        arg = x*scale

! Suppose NOTHING bad happens in gratio

        CALL gratio(shape,arg,cum,ccum,0)

        RETURN

      END SUBROUTINE local_cum_gamma

!*********************************************************************

    END MODULE cdf_gamma_mod
