    MODULE biomath_mathlib_mod
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PUBLIC
! ..
! FUNCTION log_beta(a0,b0) -- An exact renamed copy of betaln
! FUNCTION log_gamma(a) -- An exact renamed copy of gamln
! FUNCTION log_bicoef( k, n )
! FUNCTION algdiv(a,b)
! FUNCTION alnrel(a)
! FUNCTION bcorr(a0,b0)
! FUNCTION betaln(a0,b0)
! FUNCTION brcomp(a,b,x,y)
! FUNCTION erf(x)
! FUNCTION erfc1(ind,x)
! FUNCTION exparg(l)
! FUNCTION gam1(a)
! FUNCTION gamln(a)
! FUNCTION gamln1(a)
! FUNCTION gsumln(a,b)
! FUNCTION rexp(x)
! FUNCTION rlog1(x)
    CONTAINS

!*********************************************************************

      FUNCTION algdiv(a,b)
!-----------------------------------------------------------------------
!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8
!                         --------
!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: algdiv
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b
! ..
! .. Local Scalars ..
        REAL (dpkind) :: c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: c0 = .833333333333333E-01_dpkind
        REAL (dpkind), PARAMETER :: c1 = -.277777777760991E-02_dpkind
        REAL (dpkind), PARAMETER :: c2 = .793650666825390E-03_dpkind
        REAL (dpkind), PARAMETER :: c3 = -.595202931351870E-03_dpkind
        REAL (dpkind), PARAMETER :: c4 = .837308034031215E-03_dpkind
        REAL (dpkind), PARAMETER :: c5 = -.165322962780713E-02_dpkind
! ..
!------------------------
        IF (a>b) THEN
          h = b/a
          c = one/(one+h)
          x = h/(one+h)
          d = a + (b-half)
        ELSE
          h = a/b
          c = h/(one+h)
          x = one/(one+h)
          d = b + (a-half)
        END IF

! SET SN = (1 - X**N)/(1 - X)

        x2 = x*x
        s3 = one + (x+x2)
        s5 = one + (x+x2*s3)
        s7 = one + (x+x2*s5)
        s9 = one + (x+x2*s7)
        s11 = one + (x+x2*s9)

! SET W = DEL(B) - DEL(A + B)

        t = (one/b)**2
        w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t + c0
        w = w*(c/b)

! COMBINE THE RESULTS

        u = d*alnrel(a/b)
        v = a*(LOG(b)-one)

        IF (u>v) THEN
          algdiv = (w-v) - u
        ELSE
          algdiv = (w-u) - v
        END IF

      END FUNCTION algdiv

!*********************************************************************

      FUNCTION alngam(x)
!**********************************************************************
!     DOUBLE PRECISION FUNCTION ALNGAM(X)
!                 double precision LN of the GAMma function
!                              Function
!     Returns the natural logarithm of GAMMA(X).
!                              Arguments
!     X --> value at which scaled log gamma is to be returned
!                    X is DOUBLE PRECISION
!                              Method
!     If X .le. 6.0, then use recursion to get X below 3
!     then apply rational approximation number 5236 of
!     Hart et al, Computer Approximations, John Wiley and Sons, NY, 1968.
!     If X .gt. 6.0, then use recursion to get X to at least 12 and
!     then use formula 5423 of the same source.
!**********************************************************************
! .. Parameters ..
        REAL (dpkind), PARAMETER :: coef(5) = (/ &
          0.83333333333333023564E-1_dpkind, - &
          0.27777777768818808E-2_dpkind, 0.79365006754279E-3_dpkind, &
          -0.594997310889E-3_dpkind, 0.8065880899E-3_dpkind/)
        REAL (dpkind), PARAMETER :: scoefd(4) = (/ &
          0.62003838007126989331E2_dpkind, 0.9822521104713994894E1_dpkind &
          , -0.8906016659497461257E1_dpkind, &
          0.1000000000000000000E1_dpkind/)
        REAL (dpkind), PARAMETER :: scoefn(9) = (/ &
          0.62003838007127258804E2_dpkind, &
          0.36036772530024836321E2_dpkind, &
          0.20782472531792126786E2_dpkind, 0.6338067999387272343E1_dpkind &
          , 0.215994312846059073E1_dpkind, 0.3980671310203570498E0_dpkind &
          , 0.1093115956710439502E0_dpkind, 0.92381945590275995E-2_dpkind &
          , 0.29737866448101651E-2_dpkind/)
        REAL (dpkind), PARAMETER :: hln2pi = &
          0.91893853320467274178E0_dpkind
        REAL (dpkind), PARAMETER :: two = 2.0_dpkind
! ..
! .. Function Return Value ..
        REAL (dpkind) :: alngam
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: offset, prod, xx
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC INT, LOG, REAL
! ..
        IF (x<=six) THEN
          prod = one
          xx = x
          IF (x>three) THEN
10          CONTINUE
            IF ( .NOT. xx>three) GO TO 20
            xx = xx - one
            prod = prod*xx
            GO TO 10

20          CONTINUE
          END IF

          IF (x<two) THEN
30          CONTINUE
            IF ( .NOT. xx<two) GO TO 40
            prod = prod/xx
            xx = xx + one
            GO TO 30

40          CONTINUE
          END IF

          alngam = evaluate_polynomial(scoefn,xx-two)/ &
            evaluate_polynomial(scoefd,xx-two)

!     COMPUTE RATIONAL APPROXIMATION TO GAMMA(X)

          alngam = alngam*prod
          alngam = LOG(alngam)
        ELSE

          offset = hln2pi

!     IF NECESSARY MAKE X AT LEAST 12 AND CARRY CORRECTION IN OFFSET

          n = INT(12.0E0_dpkind-x)

          IF (n>0) THEN
            prod = one
            DO i = 1, n
              prod = prod*(x+REAL(i-1,kind=dpkind))
            END DO

            offset = offset - LOG(prod)
            xx = x + REAL(n,kind=dpkind)
          ELSE

            xx = x

!     COMPUTE POWER SERIES

          END IF

          alngam = evaluate_polynomial(coef,one/xx**2)/xx
          alngam = alngam + offset + (xx-half)*LOG(xx) - xx
        END IF

      END FUNCTION alngam

!*********************************************************************

      FUNCTION alnrel(a)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: alnrel
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a
! ..
! .. Local Scalars ..
        REAL (dpkind) :: t, t2, w
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: p(4) = (/ one, &
          -.129418923021993E+01_dpkind, .405303492862024E+00_dpkind, &
          -.178874546012214E-01_dpkind/)
        REAL (dpkind), PARAMETER :: q(4) = (/ one, &
          -.162752256355323E+01_dpkind, .747811014037616E+00_dpkind, &
          -.845104217945565E-01_dpkind/)
! ..
!--------------------------
        IF (ABS(a)<=0.375_dpkind) THEN
          t = a/(a+two)
          t2 = t*t

          w = evaluate_polynomial(p,t2)/evaluate_polynomial(q,t2)

          alnrel = two*t*w

        ELSE IF (a<-one) THEN
          alnrel = -one
        ELSE
          alnrel = LOG(a+one)
        END IF

      END FUNCTION alnrel

!*********************************************************************

      FUNCTION apser(a,b,x,eps)
!-----------------------------------------------------------------------
!     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
!     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
!     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: apser
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: aj, bx, c, j, s, t, tol
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: g = .577215664901533_dpkind
! ..
        bx = b*x
        t = x - bx

        IF (b*eps<=2.E-2_dpkind) THEN
          c = LOG(x) + psi(b) + g + t
        ELSE
          c = LOG(bx) + g + t
        END IF

        tol = five*eps*ABS(c)
        j = one
        s = zero

        DO
          j = j + one
          t = t*(x-bx/j)
          aj = t/j
          s = s + aj

          IF (ABS(aj)<=tol) EXIT
        END DO

        apser = -a*(c+s)

        RETURN

      END FUNCTION apser

!*********************************************************************

      FUNCTION basym(a,b,lambda,eps)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
!     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
!     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
!     A AND B ARE GREATER THAN OR EQUAL TO 15.
!-----------------------------------------------------------------------
!------------------------
!     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
!            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
!            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
!------------------------
!     E0 = 2/SQRT(PI)
!     E1 = 2**(-3/2)
!------------------------
! .. Function Return Value ..
        REAL (dpkind) :: basym
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, lambda
! ..
! .. Local Scalars ..
        REAL (dpkind) :: bsum, dsum, f, h, h2, hn, j0, j1, r, r0, r1, s, &
          sum, t, t0, t1, u, w, w0, z, z0, z2, zn, znm1
        INTEGER :: i, im1, imj, j, m, mm1, mmj, n, np1
! ..
! .. Local Arrays ..
        REAL (dpkind) :: a0(21), b0(21), c(21), d(21)
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, SQRT
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: e0 = 1.12837916709551_dpkind
        REAL (dpkind), PARAMETER :: e1 = .353553390593274_dpkind
        INTEGER, PARAMETER :: num = 20
! ..
        basym = zero

        IF (a<b) THEN
          h = a/b
          r0 = one/(one+h)
          r1 = (b-a)/b
          w0 = one/SQRT(a*(one+h))
        ELSE
          h = b/a
          r0 = one/(one+h)
          r1 = (b-a)/a
          w0 = one/SQRT(b*(one+h))
        END IF

        f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
        t = EXP(-f)

        IF (t==zero) RETURN

        z0 = SQRT(f)
        z = half*(z0/e1)
        z2 = f + f

        a0(1) = (two/three)*r1
        c(1) = -half*a0(1)
        d(1) = -c(1)
        j0 = (half/e0)*erfc1(1,z0)
        j1 = e1
        sum = j0 + d(1)*w0*j1

        s = one
        h2 = h*h
        hn = one
        w = w0
        znm1 = z
        zn = z2
        DO n = 2, num, 2
          hn = h2*hn
          a0(n) = two*r0*(one+h*hn)/(n+two)
          np1 = n + 1
          s = s + hn
          a0(np1) = two*r1*s/(n+three)

          DO i = n, np1
            r = -half*(i+one)
            b0(1) = r*a0(1)
            DO m = 2, i
              bsum = zero
              mm1 = m - 1
              DO j = 1, mm1
                mmj = m - j
                bsum = bsum + (j*r-mmj)*a0(j)*b0(mmj)
              END DO
              b0(m) = r*a0(m) + bsum/m
            END DO

            c(i) = b0(i)/(i+one)

            dsum = zero
            im1 = i - 1
            DO j = 1, im1
              imj = i - j
              dsum = dsum + d(imj)*c(j)
            END DO
            d(i) = -(dsum+c(i))
          END DO

          j0 = e1*znm1 + (n-one)*j0
          j1 = e1*zn + n*j1
          znm1 = z2*znm1
          zn = z2*zn
          w = w0*w
          t0 = d(n)*w*j0
          w = w0*w
          t1 = d(np1)*w*j1
          sum = sum + (t0+t1)
          IF (ABS(t0)+ABS(t1)<=eps*sum) EXIT
        END DO

        u = EXP(-bcorr(a,b))

        basym = e0*t*u*sum

        RETURN

      END FUNCTION basym

!*********************************************************************

      FUNCTION bcorr(a0,b0)
!-----------------------------------------------------------------------
!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: bcorr
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: a0, b0
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a, b, c, h, s11, s3, s5, s7, s9, t, w, x, x2
! ..
! .. Intrinsic Functions ..
        INTRINSIC MAX, MIN
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: c0 = .833333333333333E-01_dpkind
        REAL (dpkind), PARAMETER :: c1 = -.277777777760991E-02_dpkind
        REAL (dpkind), PARAMETER :: c2 = .793650666825390E-03_dpkind
        REAL (dpkind), PARAMETER :: c3 = -.595202931351870E-03_dpkind
        REAL (dpkind), PARAMETER :: c4 = .837308034031215E-03_dpkind
        REAL (dpkind), PARAMETER :: c5 = -.165322962780713E-02_dpkind
! ..
!------------------------
        a = MIN(a0,b0)
        b = MAX(a0,b0)

        h = a/b
        c = h/(one+h)
        x = one/(one+h)
        x2 = x*x

!                SET SN = (1 - X**N)/(1 - X)

        s3 = one + (x+x2)
        s5 = one + (x+x2*s3)
        s7 = one + (x+x2*s5)
        s9 = one + (x+x2*s7)
        s11 = one + (x+x2*s9)

!                SET W = DEL(B) - DEL(A + B)

        t = (one/b)**2
        w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t + c0
        w = w*(c/b)

!                   COMPUTE  DEL(A) + W

        t = (one/a)**2

        bcorr = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a + w

      END FUNCTION bcorr

!*********************************************************************

      FUNCTION betaln(a0,b0)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
! .. Function Return Value ..
        REAL (dpkind) :: betaln
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: a0, b0
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a, b, c, h, u, v, w, z
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG, MAX, MIN
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: e = .918938533204673_dpkind
! ..
!--------------------------
        a = MIN(a0,b0)
        b = MAX(a0,b0)

        IF (a>=eight) GO TO 50
        IF (a>=one) GO TO 10

!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .LT. 1
!-----------------------------------------------------------------------
        IF (b<eight) THEN
          betaln = gamln(a) + (gamln(b)-gamln(a+b))
          RETURN

        END IF

        betaln = gamln(a) + algdiv(a,b)

        RETURN

!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 .LE. A .LT. 8
!-----------------------------------------------------------------------
10      CONTINUE
        IF (a>two) GO TO 20
        IF (b<=two) THEN
          betaln = gamln(a) + gamln(b) - gsumln(a,b)
          RETURN

        END IF

        w = zero
        IF (b<eight) GO TO 30

        betaln = gamln(a) + algdiv(a,b)
        RETURN

!                REDUCTION OF A WHEN B .LE. 1000

20      CONTINUE
        IF (b>thousand) GO TO 40
        n = a - one
        w = one

        DO i = 1, n
          a = a - one
          h = a/b
          w = w*(h/(one+h))
        END DO

        w = LOG(w)

        IF (b>=eight) THEN
          betaln = w + gamln(a) + algdiv(a,b)
          RETURN

!                 REDUCTION OF B WHEN B .LT. 8

        END IF

30      CONTINUE
        n = b - one
        z = one

        DO i = 1, n
          b = b - one
          z = z*(b/(a+b))
        END DO

        betaln = w + LOG(z) + (gamln(a)+(gamln(b)-gsumln(a,b)))
        RETURN

!                REDUCTION OF A WHEN B .GT. 1000

40      CONTINUE
        n = a - one
        w = one

        DO i = 1, n
          a = a - one
          w = w*(a/(one+a/b))
        END DO

        betaln = (LOG(w)-n*LOG(b)) + (gamln(a)+algdiv(a,b))
        RETURN

!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .GE. 8
!-----------------------------------------------------------------------
50      CONTINUE
        w = bcorr(a,b)
        h = a/b
        c = h/(one+h)
        u = -(a-half)*LOG(c)
        v = b*alnrel(h)

        IF (u>v) THEN
          betaln = ((((-half*LOG(b))+e)+w)-v) - u
          RETURN

        END IF

        betaln = ((((-half*LOG(b))+e)+w)-u) - v

      END FUNCTION betaln

!*********************************************************************

      FUNCTION bfrac(a,b,x,y,lambda,eps)
!-----------------------------------------------------------------------
!     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
!     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: bfrac
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, lambda, x, y
! ..
! .. Local Scalars ..
        REAL (dpkind) :: alpha, an, anp1, beta, bn, bnp1, c, c0, c1, e, &
          n, p, r, r0, s, t, w, yp1
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS
! ..
        bfrac = brcomp(a,b,x,y)

        IF (bfrac==zero) RETURN

        c = one + lambda
        c0 = b/a
        c1 = one + one/a
        yp1 = y + one

        n = zero
        p = one
        s = a + one
        an = zero
        bn = one
        anp1 = one
        bnp1 = c/c1
        r = c1/c

!        CONTINUED FRACTION CALCULATION

10      CONTINUE
        n = n + one
        t = n/a
        w = n*(b-n)*x
        e = a/s
        alpha = (p*(p+c0)*e*e)*(w*x)
        e = (one+t)/(c1+t+t)
        beta = n + w/s + e*(c+n*yp1)
        p = one + t
        s = s + two

!        UPDATE AN, BN, ANP1, AND BNP1

        t = alpha*an + beta*anp1
        an = anp1
        anp1 = t
        t = alpha*bn + beta*bnp1
        bn = bnp1
        bnp1 = t

        r0 = r
        r = anp1/bnp1
        IF (ABS(r-r0)<=eps*r) GO TO 20

!        RESCALE AN, BN, ANP1, AND BNP1

        an = an/bnp1
        bn = bn/bnp1
        anp1 = r
        bnp1 = one
        GO TO 10

20      CONTINUE
        bfrac = bfrac*r

        RETURN

      END FUNCTION bfrac

!*********************************************************************

      SUBROUTINE bgrat(a,b,x,y,w,eps,ierr)
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, x, y
        REAL (dpkind), INTENT (INOUT) :: w
        INTEGER, INTENT (OUT) :: ierr
! ..
! .. Local Scalars ..
        REAL (dpkind) :: bm1, bp2n, cn, coef, dj, j, l, lnx, n2, nu, p, &
          q, r, s, sum, t, t2, u, v, z
        INTEGER :: i, n, nm1
! ..
! .. Local Arrays ..
        REAL (dpkind) :: c(30), d(30)
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, LOG
! ..
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
!     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
!     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!-----------------------------------------------------------------------
        bm1 = (b-half) - half
        nu = a + half*bm1

        IF (y<=0.375_dpkind) THEN
          lnx = alnrel(-y)
        ELSE
          lnx = LOG(x)
        END IF

        z = -nu*lnx
        IF (b*z==zero) GO TO 10

!                 COMPUTATION OF THE EXPANSION
!                 SET R = EXP(-Z)*Z**B/GAMMA(B)

        r = b*(one+gam1(b))*EXP(b*LOG(z))
        r = r*EXP(a*lnx)*EXP(half*bm1*lnx)
        u = algdiv(b,a) + b*LOG(nu)
        u = r*EXP(-u)

        IF (u/=zero) THEN
          CALL grat1(b,z,r,p,q,eps)

          v = fourth*(one/nu)**2
          t2 = fourth*lnx*lnx
          l = w/u
          j = q/r
          sum = j
          t = one
          cn = one
          n2 = zero
          DO n = 1, 30
            bp2n = b + n2
            j = (bp2n*(bp2n+one)*j+(z+bp2n+one)*t)*v
            n2 = n2 + two
            t = t*t2
            cn = cn/(n2*(n2+one))
            c(n) = cn
            s = zero

            IF (n/=1) THEN
              nm1 = n - 1
              coef = b - n
              DO i = 1, nm1
                s = s + coef*c(i)*d(n-i)
                coef = coef + b
              END DO
            END IF

            d(n) = bm1*cn + s/n
            dj = d(n)*j
            sum = sum + dj
            IF (sum<=zero) GO TO 10
            IF (ABS(dj)<=eps*(sum+l)) EXIT
          END DO

! ADD THE RESULTS TO W

          ierr = 0
          w = w + u*sum
          RETURN

        END IF

! THE EXPANSION CANNOT BE COMPUTED

10      CONTINUE
        ierr = 1

        RETURN

      END SUBROUTINE bgrat

!*********************************************************************

      FUNCTION bpser(a,b,x,eps)
!-----------------------------------------------------------------------
!     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
!     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: bpser
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a0, apb, b0, c, n, sum, t, tol, u, w, z
        INTEGER :: i, m
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, LOG, MAX, MIN, REAL
! ..
        bpser = zero
        IF (x==zero) RETURN
!-----------------------------------------------------------------------
!            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
!-----------------------------------------------------------------------
        a0 = MIN(a,b)
        IF (a0<one) GO TO 10
        z = a*LOG(x) - betaln(a,b)
        bpser = EXP(z)/a
        GO TO 40

10      CONTINUE
        b0 = MAX(a,b)
        IF (b0>=eight) GO TO 30
        IF (b0>one) GO TO 20

!            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1

        bpser = x**a
        IF (bpser==zero) RETURN

        apb = a + b
        IF (apb<=one) THEN
          z = one + gam1(apb)
        ELSE
          u = REAL(a,kind=dpkind) + REAL(b,kind=dpkind) - one
          z = (one+gam1(u))/apb
        END IF

        c = (one+gam1(a))*(one+gam1(b))/z
        bpser = bpser*c*(b/apb)
        GO TO 40

!         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8

20      CONTINUE
        u = gamln1(a0)
        m = b0 - one

        IF (m>=1) THEN
          c = one
          DO i = 1, m
            b0 = b0 - one
            c = c*(b0/(a0+b0))
          END DO

          u = LOG(c) + u
        END IF

        z = a*LOG(x) - u
        b0 = b0 - one
        apb = a0 + b0

        IF (apb<=one) THEN
          t = one + gam1(apb)
        ELSE
          u = REAL(a0,kind=dpkind) + REAL(b0,kind=dpkind) - one
          t = (one+gam1(u))/apb
        END IF

        bpser = EXP(z)*(a0/a)*(one+gam1(b0))/t
        GO TO 40

!            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8

30      CONTINUE
        u = gamln1(a0) + algdiv(a0,b0)
        z = a*LOG(x) - u
        bpser = (a0/a)*EXP(z)

40      CONTINUE
        IF (bpser==zero .OR. a<=0.1_dpkind*eps) RETURN
!-----------------------------------------------------------------------
!                     COMPUTE THE SERIES
!-----------------------------------------------------------------------
        sum = zero
        n = zero
        c = one
        tol = eps/a

50      CONTINUE
        n = n + one
        c = c*(half+(half-b/n))*x
        w = c/(a+n)
        sum = sum + w
        IF (ABS(w)>tol) GO TO 50

        bpser = bpser*(one+a*sum)

        RETURN

      END FUNCTION bpser

!*********************************************************************

      SUBROUTINE bratio(a,b,x,y,w,w1,ierr)
!----------------------------------------------------------------------
!            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)
!                     --------------------
!     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
!     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES
!                      W  = IX(A,B)
!                      W1 = 1 - IX(A,B)
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
!     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
!     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
!     ONE OF THE FOLLOWING VALUES ...
!        IERR = 1  IF A OR B IS NEGATIVE
!        IERR = 2  IF A = B = 0
!        IERR = 3  IF X .LT. 0 OR X .GT. 1
!        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
!        IERR = 5  IF X + Y .NE. 1
!        IERR = 6  IF X = A = 0
!        IERR = 7  IF Y = B = 0
!--------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     REVISED ... NOV 1991
!-----------------------------------------------------------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, x, y
        REAL (dpkind), INTENT (OUT) :: w, w1
        INTEGER, INTENT (OUT) :: ierr
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a0, b0, eps, lambda, t, x0, y0, z
        INTEGER :: ierr1, ind, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EPSILON, MAX, MIN
! ..
        eps = EPSILON(one)

        w = zero
        w1 = zero

        IF (a<zero .OR. b<zero) THEN
          ierr = 1
          RETURN
        END IF

        IF (a==zero .AND. b==zero) THEN
          ierr = 2
          RETURN
        END IF

        IF (x<zero .OR. x>one) THEN
          ierr = 3
          RETURN
        END IF

        IF (y<zero .OR. y>one) THEN
          ierr = 4
          RETURN
        END IF

        z = ((x+y)-half) - half

        IF (ABS(z)>three*eps) THEN
          ierr = 5
          RETURN
        END IF

        ierr = 0
        IF (x==zero) GO TO 150
        IF (y==zero) GO TO 170
        IF (a==zero) GO TO 180
        IF (b==zero) GO TO 160

        eps = MAX(eps,1.E-15_dpkind)

        IF (MAX(a,b)<1.E-3_dpkind*eps) GO TO 200

        ind = 0
        a0 = a
        b0 = b
        x0 = x
        y0 = y
        IF (MIN(a0,b0)>one) GO TO 30

!             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1

        IF (x>half) THEN
          ind = 1
          a0 = b
          b0 = a
          x0 = y
          y0 = x
        END IF

        IF (b0<MIN(eps,eps*a0)) GO TO 50
        IF (a0<MIN(eps,eps*b0) .AND. b0*x0<=one) GO TO 60
        IF (MAX(a0,b0)>one) GO TO 10
        IF (a0>=MIN(0.2_dpkind,b0)) GO TO 70
        IF (x0**a0<=0.9_dpkind) GO TO 70
        IF (x0>=0.3_dpkind) GO TO 80
        n = 20
        GO TO 100

10      CONTINUE
        IF (b0<=one) GO TO 70
        IF (x0>=0.3_dpkind) GO TO 80
        IF (x0>=0.1_dpkind) GO TO 20
        IF ((x0*b0)**a0<=0.7_dpkind) GO TO 70

20      CONTINUE
        IF (b0>15.0_dpkind) GO TO 110
        n = 20
        GO TO 100

!             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1

30      CONTINUE
        IF (a<=b) THEN
          lambda = a - (a+b)*x
        ELSE
          lambda = (a+b)*y - b
        END IF

        IF (lambda<zero) THEN
          ind = 1
          a0 = b
          b0 = a
          x0 = y
          y0 = x
          lambda = ABS(lambda)
        END IF

        IF (b0<40.0_dpkind .AND. b0*x0<=0.7_dpkind) GO TO 70
        IF (b0<40.0_dpkind) GO TO 120
        IF (a0>b0) GO TO 40
        IF (a0<=hundred) GO TO 90
        IF (lambda>0.03_dpkind*a0) GO TO 90
        GO TO 140

40      CONTINUE
        IF (b0<=hundred) GO TO 90
        IF (lambda>0.03_dpkind*b0) GO TO 90
        GO TO 140

!            EVALUATION OF THE APPROPRIATE ALGORITHM

50      CONTINUE
        w = fpser(a0,b0,x0,eps)
        w1 = half + (half-w)
        GO TO 190

60      CONTINUE
        w1 = apser(a0,b0,x0,eps)
        w = half + (half-w1)
        GO TO 190

70      CONTINUE
        w = bpser(a0,b0,x0,eps)
        w1 = half + (half-w)
        GO TO 190

80      CONTINUE
        w1 = bpser(b0,a0,y0,eps)
        w = half + (half-w1)
        GO TO 190

90      CONTINUE
        w = bfrac(a0,b0,x0,y0,lambda,15.0_dpkind*eps)

        w1 = half + (half-w)
        GO TO 190

100     CONTINUE
        w1 = bup(b0,a0,y0,x0,n,eps)
        b0 = b0 + n

110     CONTINUE
        CALL bgrat(b0,a0,y0,x0,w1,15.0_dpkind*eps,ierr1)

        w = half + (half-w1)
        GO TO 190

120     CONTINUE
        n = b0
        b0 = b0 - n
        IF (b0==zero) THEN
          n = n - 1
          b0 = one
        END IF

        w = bup(b0,a0,y0,x0,n,eps)

        IF (x0>0.7_dpkind) GO TO 130

        w = w + bpser(a0,b0,x0,eps)
        w1 = half + (half-w)
        GO TO 190

130     CONTINUE
        IF (a0<=15.0_dpkind) THEN
          n = 20
          w = w + bup(a0,b0,x0,y0,n,eps)
          a0 = a0 + n
        END IF

        CALL bgrat(a0,b0,x0,y0,w,15.0_dpkind*eps,ierr1)

        w1 = half + (half-w)
        GO TO 190

140     CONTINUE

        w = basym(a0,b0,lambda,100.0_dpkind*eps)

        w1 = half + (half-w)
        GO TO 190

!               TERMINATION OF THE PROCEDURE

150     CONTINUE
        IF (a==zero) THEN
          ierr = 6
          RETURN
        END IF

160     CONTINUE
        w = zero
        w1 = one
        RETURN

170     CONTINUE
        IF (b==zero) THEN
          ierr = 7
          RETURN
        END IF

180     CONTINUE
        w = one
        w1 = zero
        RETURN

190     CONTINUE
        IF (ind==0) RETURN
        t = w
        w = w1
        w1 = t
        RETURN

!           PROCEDURE FOR A AND B .LT. 1.E-3*EPS

200     CONTINUE
        w = b/(a+b)
        w1 = a/(a+b)

        RETURN

      END SUBROUTINE bratio

!*********************************************************************

      FUNCTION brcmp1(mu,a,b,x,y)
!-----------------------------------------------------------------------
!          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
!-----------------------------------------------------------------------
!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------
! .. Function Return Value ..
        REAL (dpkind) :: brcmp1
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, x, y
        INTEGER, INTENT (IN) :: mu
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a0, apb, b0, c, e, h, lambda, lnx, lny, t, u, v, &
          x0, y0, z
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, LOG, MAX, MIN, SQRT
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: const = .398942280401433_dpkind
! ..
        a0 = MIN(a,b)

        IF (a0>=eight) GO TO 30

        IF (x<=0.375_dpkind) THEN
          lnx = LOG(x)
          lny = alnrel(-x)
        ELSE

          IF (y<=0.375_dpkind) THEN
            lnx = alnrel(-y)
            lny = LOG(y)
          ELSE
            lnx = LOG(x)
            lny = LOG(y)
          END IF
        END IF

        z = a*lnx + b*lny

        IF (a0>=one) THEN
          z = z - betaln(a,b)
          brcmp1 = esum(mu,z)

          RETURN
        END IF

!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
        b0 = MAX(a,b)
        IF (b0>=eight) GO TO 20
        IF (b0>one) GO TO 10

!                   ALGORITHM FOR B0 .LE. 1

        brcmp1 = esum(mu,z)

        IF (brcmp1==zero) RETURN

        apb = a + b
        IF (apb<=one) THEN
          z = one + gam1(apb)
        ELSE
          u = a + b - one
          z = (one+gam1(u))/apb
        END IF

        c = (one+gam1(a))*(one+gam1(b))/z
        brcmp1 = brcmp1*(a0*c)/(one+a0/b0)

        RETURN

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

10      CONTINUE
        u = gamln1(a0)
        n = b0 - one

        IF (n>=1) THEN
          c = one
          DO i = 1, n
            b0 = b0 - one
            c = c*(b0/(a0+b0))
          END DO
          u = LOG(c) + u

        END IF

        z = z - u
        b0 = b0 - one
        apb = a0 + b0

        IF (apb<=one) THEN
          t = one + gam1(apb)
        ELSE
          u = a0 + b0 - one
          t = (one+gam1(u))/apb
        END IF

        brcmp1 = a0*esum(mu,z)*(one+gam1(b0))/t

        RETURN

!                   ALGORITHM FOR B0 .GE. 8

20      CONTINUE
        u = gamln1(a0) + algdiv(a0,b0)
        brcmp1 = a0*esum(mu,z-u)

        RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .GE. 8 AND B .GE. 8
!-----------------------------------------------------------------------
30      CONTINUE
        IF (a<=b) THEN
          h = a/b
          x0 = h/(one+h)
          y0 = one/(one+h)
          lambda = a - (a+b)*x
        ELSE
          h = b/a
          x0 = one/(one+h)
          y0 = h/(one+h)
          lambda = (a+b)*y - b
        END IF

        e = -lambda/a
        IF (ABS(e)<=0.6_dpkind) THEN
          u = rlog1(e)
        ELSE
          u = e - LOG(x/x0)
        END IF

        e = lambda/b

        IF (ABS(e)<=0.6_dpkind) THEN
          v = rlog1(e)
        ELSE
          v = e - LOG(y/y0)
        END IF

        z = esum(mu,-(a*u+b*v))

        brcmp1 = const*SQRT(b*x0)*z*EXP(-bcorr(a,b))

      END FUNCTION brcmp1

!*********************************************************************

      FUNCTION brcomp(a,b,x,y)
!-----------------------------------------------------------------------
!               EVALUATION OF X**A*Y**B/BETA(A,B)
!-----------------------------------------------------------------------
!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------
! .. Function Return Value ..
        REAL (dpkind) :: brcomp
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, x, y
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a0, apb, b0, c, e, h, lambda, lnx, lny, t, u, v, &
          x0, y0, z
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, LOG, MAX, MIN, SQRT
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: const = .398942280401433_dpkind
! ..
        brcomp = zero

        IF (x<=zero .OR. y<=zero) RETURN

        a0 = MIN(a,b)

        IF (a0>=8.0_dpkind) GO TO 30

        IF (x<=0.375_dpkind) THEN
          lnx = LOG(x)
          lny = alnrel(-x)
        ELSE
          IF (y<=0.375_dpkind) THEN
            lnx = alnrel(-y)
            lny = LOG(y)
          ELSE
            lnx = LOG(x)
            lny = LOG(y)
          END IF
        END IF

        z = a*lnx + b*lny

        IF (a0>=one) THEN
          z = z - betaln(a,b)
          brcomp = EXP(z)

          RETURN
        END IF

!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------

        b0 = MAX(a,b)

        IF (b0>=eight) GO TO 20
        IF (b0>one) GO TO 10

!                   ALGORITHM FOR B0 .LE. 1

        brcomp = EXP(z)

        IF (brcomp==zero) RETURN

        apb = a + b

        IF (apb<=one) THEN
          z = one + gam1(apb)
        ELSE
          u = a + b - one
          z = (one+gam1(u))/apb
        END IF

        c = (one+gam1(a))*(one+gam1(b))/z
        brcomp = brcomp*(a0*c)/(one+a0/b0)
        RETURN

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

10      CONTINUE
        u = gamln1(a0)
        n = b0 - one

        IF (n>=1) THEN
          c = one
          DO i = 1, n
            b0 = b0 - one
            c = c*(b0/(a0+b0))
          END DO

          u = LOG(c) + u

        END IF

        z = z - u
        b0 = b0 - one
        apb = a0 + b0

        IF (apb<=one) THEN
          t = one + gam1(apb)
        ELSE
          u = a0 + b0 - one
          t = (one+gam1(u))/apb
        END IF

        brcomp = a0*EXP(z)*(one+gam1(b0))/t
        RETURN

!                   ALGORITHM FOR B0 .GE. 8

20      CONTINUE
        u = gamln1(a0) + algdiv(a0,b0)
        brcomp = a0*EXP(z-u)
        RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .GE. 8 AND B .GE. 8
!-----------------------------------------------------------------------
30      CONTINUE
        IF (a<=b) THEN
          h = a/b
          x0 = h/(one+h)
          y0 = one/(one+h)
          lambda = a - (a+b)*x
        ELSE
          h = b/a
          x0 = one/(one+h)
          y0 = h/(one+h)
          lambda = (a+b)*y - b
        END IF

        e = -lambda/a

        IF (ABS(e)<=0.6_dpkind) THEN
          u = rlog1(e)
        ELSE
          u = e - LOG(x/x0)
        END IF

        e = lambda/b

        IF (ABS(e)<=0.6_dpkind) THEN
          v = rlog1(e)
        ELSE
          v = e - LOG(y/y0)
        END IF
        z = EXP(-(a*u+b*v))

        brcomp = const*SQRT(b*x0)*z*EXP(-bcorr(a,b))

      END FUNCTION brcomp

!*********************************************************************

      FUNCTION bup(a,b,x,y,n,eps)
!-----------------------------------------------------------------------
!     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
!     EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: bup
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, x, y
        INTEGER, INTENT (IN) :: n
! ..
! .. Local Scalars ..
        REAL (dpkind) :: ap1, apb, d, l, r, t, w
        INTEGER :: i, k, kp1, mu, nm1
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, MIN
! ..
!          OBTAIN THE SCALING FACTOR EXP(-MU) AND
!             EXP(MU)*(X**A*Y**B/BETA(A,B))/A
        apb = a + b
        ap1 = a + one
        mu = 0
        d = one

        IF (n/=1 .AND. a>=one) THEN
          IF (apb>=1.1_dpkind*ap1) THEN
            mu = ABS(exparg(1))
            k = exparg(0)
            mu = MIN(k,mu)
            t = mu
            d = EXP(-t)

          END IF
        END IF

        bup = brcmp1(mu,a,b,x,y)/a

        IF (n==1 .OR. bup==zero) RETURN

        nm1 = n - 1
        w = d

!          LET K BE THE INDEX OF THE MAXIMUM TERM

        k = 0
        IF (b<=one) GO TO 30
        IF (y>1.E-4_dpkind) GO TO 10
        k = nm1
        GO TO 20

10      CONTINUE
        r = (b-one)*x/y - a
        IF (r<one) GO TO 30
        k = nm1
        t = nm1
        IF (r<t) k = r

!          ADD THE INCREASING TERMS OF THE SERIES

20      CONTINUE
        DO i = 1, k
          l = i - 1
          d = ((apb+l)/(ap1+l))*x*d
          w = w + d
        END DO

        IF (k==nm1) GO TO 40

!          ADD THE REMAINING TERMS OF THE SERIES

30      CONTINUE
        kp1 = k + 1

        DO i = kp1, nm1
          l = i - 1
          d = ((apb+l)/(ap1+l))*x*d
          w = w + d
          IF (d<=eps*w) EXIT
        END DO

40      CONTINUE
        bup = bup*w

        RETURN

      END FUNCTION bup

!*********************************************************************

      FUNCTION erf(x)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: erf
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: ax, bot, t, top, x2
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, SIGN
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: a(5) = (/ .128379167095513E+00_dpkind &
          , .479137145607681E-01_dpkind, .323076579225834E-01_dpkind, &
          -.133733772997339E-02_dpkind, .771058495001320E-04_dpkind/)
        REAL (dpkind), PARAMETER :: b(4) = (/ one, &
          .375795757275549E+00_dpkind, .538971687740286E-01_dpkind, &
          .301048631703895E-02_dpkind/)
        REAL (dpkind), PARAMETER :: p(8) = (/ 3.00459261020162E+02_dpkind &
          , 4.51918953711873E+02_dpkind, 3.39320816734344E+02_dpkind, &
          1.52989285046940E+02_dpkind, 4.31622272220567E+01_dpkind, &
          7.21175825088309E+00_dpkind, 5.64195517478974E-01_dpkind, &
          -1.36864857382717E-07_dpkind/)
        REAL (dpkind), PARAMETER :: q(8) = (/ 3.00459260956983E+02_dpkind &
          , 7.90950925327898E+02_dpkind, 9.31354094850610E+02_dpkind, &
          6.38980264465631E+02_dpkind, 2.77585444743988E+02_dpkind, &
          7.70001529352295E+01_dpkind, 1.27827273196294E+01_dpkind, one/)
        REAL (dpkind), PARAMETER :: r(5) = (/ 2.82094791773523E-01_dpkind &
          , 4.65807828718470E+00_dpkind, 2.13688200555087E+01_dpkind, &
          2.62370141675169E+01_dpkind, 2.10144126479064E+00_dpkind/)
        REAL (dpkind), PARAMETER :: s(5) = (/ one, &
          1.80124575948747E+01_dpkind, 9.90191814623914E+01_dpkind, &
          1.87114811799590E+02_dpkind, 9.41537750555460E+01_dpkind/)
        REAL (dpkind), PARAMETER :: c = .564189583547756_dpkind
! ..
        ax = ABS(x)

        IF (ax<=half) THEN
          t = x*x

          top = evaluate_polynomial(a,t) + one
          bot = evaluate_polynomial(b,t)

          erf = x*(top/bot)

        ELSE IF (ax<=four) THEN
          top = evaluate_polynomial(p,ax)
          bot = evaluate_polynomial(q,ax)

          erf = half + (half-EXP(-x*x)*top/bot)

          IF (x<zero) erf = -erf

        ELSE IF (ax<5.8_dpkind) THEN
          x2 = x*x
          t = one/x2

          top = evaluate_polynomial(r,t)
          bot = evaluate_polynomial(s,t)

          erf = (c-top/(x2*bot))/ax
          erf = half + (half-EXP(-x2)*erf)

          IF (x<zero) erf = -erf
        ELSE
          erf = SIGN(one,x)
        END IF

      END FUNCTION erf

!*********************************************************************

      FUNCTION erfc1(ind,x)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
!          ERFC1(IND,X) = ERFC(X)            IF IND = 0
!          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: erfc1
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x
        INTEGER, INTENT (IN) :: ind
! ..
! .. Local Scalars ..
        REAL (dpkind) :: ax, bot, e, t, top, w
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: a(5) = (/ .128379167095513E+00_dpkind &
          , .479137145607681E-01_dpkind, .323076579225834E-01_dpkind, &
          -.133733772997339E-02_dpkind, .771058495001320E-04_dpkind/)
        REAL (dpkind), PARAMETER :: b(4) = (/ one, &
          .375795757275549E+00_dpkind, .538971687740286E-01_dpkind, &
          .301048631703895E-02_dpkind/)
        REAL (dpkind), PARAMETER :: p(8) = (/ 3.00459261020162E+02_dpkind &
          , 4.51918953711873E+02_dpkind, 3.39320816734344E+02_dpkind, &
          1.52989285046940E+02_dpkind, 4.31622272220567E+01_dpkind, &
          7.21175825088309E+00_dpkind, 5.64195517478974E-01_dpkind, &
          -1.36864857382717E-07_dpkind/)
        REAL (dpkind), PARAMETER :: q(8) = (/ 3.00459260956983E+02_dpkind &
          , 7.90950925327898E+02_dpkind, 9.31354094850610E+02_dpkind, &
          6.38980264465631E+02_dpkind, 2.77585444743988E+02_dpkind, &
          7.70001529352295E+01_dpkind, 1.27827273196294E+01_dpkind, one/)
        REAL (dpkind), PARAMETER :: r(5) = (/ 2.82094791773523E-01_dpkind &
          , 4.65807828718470E+00_dpkind, 2.13688200555087E+01_dpkind, &
          2.62370141675169E+01_dpkind, 2.10144126479064E+00_dpkind/)
        REAL (dpkind), PARAMETER :: s(5) = (/ one, &
          1.80124575948747E+01_dpkind, 9.90191814623914E+01_dpkind, &
          1.87114811799590E+02_dpkind, 9.41537750555460E+01_dpkind/)
        REAL (dpkind), PARAMETER :: c = .564189583547756_dpkind
! ..
!-------------------------
!                     ABS(X) .LE. 0.5
        ax = ABS(x)

        IF (ax<=half) THEN
! x in [-0.5, 0.5]

          t = x*x

          top = evaluate_polynomial(a,t) + one
          bot = evaluate_polynomial(b,t)

          erfc1 = half + (half-x*(top/bot))

          IF (ind/=0) erfc1 = EXP(t)*erfc1

          RETURN

        END IF

        IF (ax>four) GO TO 10

        top = evaluate_polynomial(p,ax)
        bot = evaluate_polynomial(q,ax)

        erfc1 = top/bot

        GO TO 30

!                      ABS(X) .GT. 4

10      CONTINUE

        IF (x<=(-5.6_dpkind)) GO TO 40

        IF (ind/=0) GO TO 20

        IF (x>hundred) GO TO 50

        IF (x*x>-exparg(1)) GO TO 50

20      CONTINUE

        t = (one/x)**2

        top = evaluate_polynomial(r,t)
        bot = evaluate_polynomial(s,t)

        erfc1 = (c-t*top/bot)/ax

!                      FINAL ASSEMBLY

30      CONTINUE

        IF (ind/=0) THEN
          IF (x<zero) erfc1 = two*EXP(x*x) - erfc1
          RETURN

        END IF

        w = x*x
        t = w
        e = w - t
        erfc1 = ((half+(half-e))*EXP(-t))*erfc1

        IF (x<zero) erfc1 = two - erfc1

        RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

40      CONTINUE

        erfc1 = two

        IF (ind/=0) erfc1 = two*EXP(x*x)

        RETURN

!   LIMIT VALUE FOR LARGE POSITIVE X (WHEN IND = 0)

50      CONTINUE

        erfc1 = zero

      END FUNCTION erfc1

!*********************************************************************

      FUNCTION esum(mu,x)
!-----------------------------------------------------------------------
!                    EVALUATION OF EXP(MU + X)
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: esum
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: x
        INTEGER :: mu
! ..
! .. Local Scalars ..
        REAL (dpkind) :: w
! ..
! .. Intrinsic Functions ..
        INTRINSIC EXP
! ..
        IF (x>zero) GO TO 10

        IF (mu<0) GO TO 20

        w = mu + x

        IF (w>zero) GO TO 20
        esum = EXP(w)
        RETURN

10      CONTINUE
        IF (mu>0) GO TO 20

        w = mu + x

        IF (w>=zero) THEN
          esum = EXP(w)
          RETURN

        END IF

20      CONTINUE
        w = mu
        esum = EXP(w)*EXP(x)

        RETURN

      END FUNCTION esum

!*********************************************************************

      FUNCTION evaluate_polynomial(a,x)
!----------------------------------------------------------------------
!              Evaluate a PoLynomial at x
!                              Function
!     Returns:
!          A(1) + A(2)*X + ... + A(N)*X**(N-1)
!                              Arguments
!     A --> Array of coefficients of the polynomial.
!     X --> Point at which the polynomial is to be evaluated.
!----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: evaluate_polynomial
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x
! ..
! .. Array Arguments ..
        REAL (dpkind), INTENT (IN) :: a(:)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: term
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC SIZE
! ..
        n = SIZE(a)
        term = a(n)

        DO i = n - 1, 1, -1
          term = a(i) + term*x
        END DO

        evaluate_polynomial = term

      END FUNCTION evaluate_polynomial

!*********************************************************************

      FUNCTION exparg(l)
! .. Function Return Value ..
        REAL (dpkind) :: exparg
! ..
! .. Scalar Arguments ..
        INTEGER :: l
! ..
! .. Intrinsic Functions ..
        INTRINSIC HUGE, LOG, TINY
! ..
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
        IF (l==0) THEN
          exparg = LOG(HUGE(one))
        ELSE
          exparg = LOG(TINY(one))
        END IF

      END FUNCTION exparg

!*********************************************************************

      FUNCTION fpser(a,b,x,eps)
!-----------------------------------------------------------------------
!                 EVALUATION OF I (A,B)
!                                X
!          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: fpser
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, b, eps, x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: an, c, s, t, tol
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, LOG
! ..
!                  SET  FPSER = X**A
        fpser = one

        IF (a>1.E-3_dpkind*eps) THEN
          fpser = zero
          t = a*LOG(x)
          IF (t<exparg(1)) RETURN
          fpser = EXP(t)
        END IF

!                NOTE THAT 1/B(A,B) = B

        fpser = (b/a)*fpser
        tol = eps/a
        an = a + one
        t = x
        s = t/an

        DO
          an = an + one
          t = x*t
          c = t/an
          s = s + c

          IF (ABS(c)<=tol) EXIT
        END DO

        fpser = fpser*(one+a*s)

        RETURN

      END FUNCTION fpser

!*********************************************************************

      FUNCTION gam1(a)
!     ------------------------------------------------------------------
!     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
!     ------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: gam1
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a
! ..
! .. Local Scalars ..
        REAL (dpkind) :: d, t, w
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: p(7) = (/ .577215664901533E+00_dpkind &
          , -.409078193005776E+00_dpkind, -.230975380857675E+00_dpkind, &
          .597275330452234E-01_dpkind, .766968181649490E-02_dpkind, &
          -.514889771323592E-02_dpkind, .589597428611429E-03_dpkind/)
        REAL (dpkind), PARAMETER :: q(5) = (/ .100000000000000E+01_dpkind &
          , .427569613095214E+00_dpkind, .158451672430138E+00_dpkind, &
          .261132021441447E-01_dpkind, .423244297896961E-02_dpkind/)
        REAL (dpkind), PARAMETER :: r(9) = (/ - &
          .422784335098468E+00_dpkind, -.771330383816272E+00_dpkind, &
          -.244757765222226E+00_dpkind, .118378989872749E+00_dpkind, &
          .930357293360349E-03_dpkind, -.118290993445146E-01_dpkind, &
          .223047661158249E-02_dpkind, .266505979058923E-03_dpkind, &
          -.132674909766242E-03_dpkind/)
        REAL (dpkind), PARAMETER :: s1 = .273076135303957E+00_dpkind
        REAL (dpkind), PARAMETER :: s2 = .559398236957378E-01_dpkind
! ..
!     -------------------
        t = a
        d = a - half

        IF (d>zero) t = d - half

        IF (t==zero) THEN
          gam1 = zero
        ELSE IF (t>zero) THEN
          w = evaluate_polynomial(p,t)/evaluate_polynomial(q,t)

          IF (d<=zero) THEN
            gam1 = a*w
          ELSE
            gam1 = (t/a)*((w-half)-half)
          END IF

        ELSE
          w = evaluate_polynomial(r,t)/((s2*t+s1)*t+one)

          IF (d<=zero) THEN
            gam1 = a*((w+half)+half)
          ELSE
            gam1 = t*w/a
          END IF
        END IF

      END FUNCTION gam1

!*********************************************************************

      FUNCTION gamln(a)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
!--------------------------
! .. Function Return Value ..
        REAL (dpkind) :: gamln
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a
! ..
! .. Local Scalars ..
        REAL (dpkind) :: t, w
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: c(6) = (/ .833333333333333E-01_dpkind &
          , -.277777777760991E-02_dpkind, .793650666825390E-03_dpkind, &
          -.595202931351870E-03_dpkind, .837308034031215E-03_dpkind, &
          -.165322962780713E-02_dpkind/)
        REAL (dpkind), PARAMETER :: d = .418938533204673_dpkind
! ..
!-----------------------------------------------------------------------
        IF (a<=zero) THEN
          gamln = -one
        ELSE IF (a<=0.8_dpkind) THEN
          gamln = gamln1(a) - LOG(a)
        ELSE IF (a<=2.25_dpkind) THEN
          t = (a-half) - half
          gamln = gamln1(t)
        ELSE IF (a<10.0_dpkind) THEN
          n = a - 1.25_dpkind
          t = a
          w = one

          DO i = 1, n
            t = t - one
            w = t*w
          END DO

          gamln = gamln1(t-one) + LOG(w)
        ELSE

          t = (one/a)**2

          w = evaluate_polynomial(c,t)/a

          gamln = (d+w) + (a-half)*(LOG(a)-one)
        END IF

      END FUNCTION gamln

!*********************************************************************

      FUNCTION gamln1(a)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: gamln1
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a
! ..
! .. Local Scalars ..
        REAL (dpkind) :: w, x
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: p(7) = (/ .577215664901533E+00_dpkind &
          , .844203922187225E+00_dpkind, -.168860593646662E+00_dpkind, &
          -.780427615533591E+00_dpkind, -.402055799310489E+00_dpkind, &
          -.673562214325671E-01_dpkind, -.271935708322958E-02_dpkind/)
        REAL (dpkind), PARAMETER :: q(7) = (/ one, &
          .288743195473681E+01_dpkind, .312755088914843E+01_dpkind, &
          .156875193295039E+01_dpkind, .361951990101499E+00_dpkind, &
          .325038868253937E-01_dpkind, .667465618796164E-03_dpkind/)
        REAL (dpkind), PARAMETER :: r(6) = (/ .422784335098467E+00_dpkind &
          , .848044614534529E+00_dpkind, .565221050691933E+00_dpkind, &
          .156513060486551E+00_dpkind, .170502484022650E-01_dpkind, &
          .497958207639485E-03_dpkind/)
        REAL (dpkind), PARAMETER :: s(6) = (/ one, &
          .124313399877507E+01_dpkind, .548042109832463E+00_dpkind, &
          .101552187439830E+00_dpkind, .713309612391000E-02_dpkind, &
          .116165475989616E-03_dpkind/)
! ..
        IF (a<0.6_dpkind) THEN
          w = evaluate_polynomial(p,a)/evaluate_polynomial(q,a)

          gamln1 = -a*w

        ELSE

          x = (a-half) - half

          w = evaluate_polynomial(r,x)/evaluate_polynomial(s,x)

          gamln1 = x*w
        END IF

      END FUNCTION gamln1

!*********************************************************************

      FUNCTION gamma(a)
!-----------------------------------------------------------------------
!         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS
!                           -----------
!     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!     BE COMPUTED.
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
!--------------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a
! ..
! .. Local Scalars ..
        REAL (dpkind) :: g, lnx, s, t, w, x, z
        INTEGER :: j, m, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, HUGE, INT, LOG, MOD, SIN
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: d = .41893853320467274178E0_dpkind
        REAL (dpkind), PARAMETER :: pi = 3.1415926535898E0_dpkind
        REAL (dpkind), PARAMETER :: p(7) = (/ one, &
          .553413866010467E+00_dpkind, .279648642639792E+00_dpkind, &
          .730981088720487E-01_dpkind, .204493667594920E-01_dpkind, &
          .261939260042690E-02_dpkind, .539637273585445E-03_dpkind/)
        REAL (dpkind), PARAMETER :: q(7) = (/ one, &
          .113062953091122E+01_dpkind, -.567902761974940E-01_dpkind, &
          -.170458969313360E+00_dpkind, .225211131035340E-01_dpkind, &
          .470059485860584E-02_dpkind, -.832979206704073E-03_dpkind/)
        REAL (dpkind), PARAMETER :: r(5) = (/ .833333333333333E-01_dpkind &
          , -.277777777770481E-02_dpkind, .793650663183693E-03_dpkind, &
          -.595156336428591E-03_dpkind, .820756370353826E-03_dpkind/)
! ..
! .. Function Return Value ..
        REAL (dpkind) :: gamma
! ..
!--------------------------
        gamma = zero
        x = a

        IF (ABS(a)>=15.0E0_dpkind) GO TO 40

!-----------------------------------------------------------------------
!            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
!-----------------------------------------------------------------------
        t = one
        m = INT(a) - 1

        IF (m<0) GO TO 10
        IF (m/=0) THEN

!     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2

          DO j = 1, m
            x = x - one
            t = x*t
          END DO
        END IF

        x = x - one
        GO TO 30

!     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1

10      CONTINUE
        t = a
        IF (a>zero) GO TO 20

        m = (-m) - 1
        IF (m/=0) THEN
          DO j = 1, m
            x = x + one
            t = x*t
          END DO
        END IF

        x = (x+half) + half
        t = x*t
        IF (t==zero) RETURN

20      CONTINUE

!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

        IF (ABS(t)>=1.E-30_dpkind) GO TO 30
        IF (ABS(t)*HUGE(one)<=1.0001E0_dpkind) RETURN

        gamma = one/t
        RETURN

!     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1

30      CONTINUE

        gamma = evaluate_polynomial(p,x)/evaluate_polynomial(q,x)

!     TERMINATION

        IF (a>=one) THEN
          gamma = gamma*t
          RETURN

        END IF

        gamma = gamma/t
        RETURN

!-----------------------------------------------------------------------
!            EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
!-----------------------------------------------------------------------
40      CONTINUE
        IF (ABS(a)>=1.E3_dpkind) RETURN

        IF (a>zero) GO TO 50

        x = -a
        n = x
        t = x - n

        IF (t>0.9E0_dpkind) t = one - t

        s = SIN(pi*t)/pi

        IF (MOD(n,2)==0) s = -s

        IF (s==zero) RETURN

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

50      CONTINUE
        t = one/(x*x)

        g = evaluate_polynomial(r,t)/x

!     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
!     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.

        lnx = LOG(x)

!     FINAL ASSEMBLY

        z = x
        g = (d+g) + (z-half)*(lnx-1.E0_dpkind)
        w = g
        t = g - w

        IF (w>0.99999E0_dpkind*exparg(0)) RETURN

        gamma = EXP(w)*(one+t)

        IF (a<zero) gamma = (one/(gamma*s))/x

        RETURN

      END FUNCTION gamma

!*********************************************************************

      SUBROUTINE grat1(a,x,r,p,q,eps)
!-----------------------------------------------------------------------
!        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                      P(A,X) AND Q(A,X)
!     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
!     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
!-----------------------------------------------------------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, eps, r, x
        REAL (dpkind), INTENT (OUT) :: p, q
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a2n, a2nm1, am0, an, an0, b2n, b2nm1, c, cma, g, &
          h, j, l, sum, t, tol, w, z
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP, LOG, SQRT
! ..
        IF (a*x==zero) GO TO 110
        IF (a==half) GO TO 100
        IF (x<1.1E0_dpkind) GO TO 10
        GO TO 60

!             TAYLOR SERIES FOR P(A,X)/X**A

10      CONTINUE
        an = three
        c = x
        sum = x/(a+three)
        tol = tenth*eps/(a+one)

20      CONTINUE
        an = an + one
        c = -c*(x/an)
        t = c/(a+an)
        sum = sum + t

        IF (ABS(t)>tol) GO TO 20

        j = a*x*((sum/six-half/(a+two))*x+one/(a+one))

        z = a*LOG(x)
        h = gam1(a)
        g = one + h

        IF (x<fourth) GO TO 30
        IF (a<x/2.59E0_dpkind) GO TO 50
        GO TO 40

30      CONTINUE
        IF (z>(-.13394E0_dpkind)) GO TO 50

40      CONTINUE
        w = EXP(z)
        p = w*g*(half+(half-j))
        q = half + (half-p)
        RETURN

50      CONTINUE
        l = rexp(z)
        w = half + (half+l)
        q = (w*j-l)*g - h
        IF (q<zero) GO TO 90
        p = half + (half-q)
        RETURN

!              CONTINUED FRACTION EXPANSION

60      CONTINUE
        a2nm1 = one
        a2n = one
        b2nm1 = x
        b2n = x + (one-a)
        c = one

70      CONTINUE
        a2nm1 = x*a2n + c*a2nm1
        b2nm1 = x*b2n + c*b2nm1
        am0 = a2nm1/b2nm1
        c = c + one
        cma = c - a
        a2n = a2nm1 + cma*a2n
        b2n = b2nm1 + cma*b2n
        an0 = a2n/b2n

        IF (ABS(an0-am0)>=eps*an0) GO TO 70

        q = r*an0
        p = half + (half-q)
        RETURN

!                SPECIAL CASES

80      CONTINUE
        p = zero
        q = one
        RETURN

90      CONTINUE
        p = one
        q = zero
        RETURN

100     CONTINUE
        IF (x<fourth) THEN
          p = erf(SQRT(x))
          q = half + (half-p)
          RETURN

        END IF

        q = erfc1(0,SQRT(x))
        p = half + (half-q)
        RETURN

110     CONTINUE
        IF (x<=a) GO TO 80
        GO TO 90

        RETURN

      END SUBROUTINE grat1

!*********************************************************************

      SUBROUTINE gratio(a,x,ans,qans,ind)
! ----------------------------------------------------------------------
!        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                      P(A,X) AND Q(A,X)
!                        ----------
!     IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
!     ARE NOT BOTH 0.
!     ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
!     P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
!     IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
!     POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
!     IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
!     6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
!     IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.
!     ERROR RETURN ...
!        ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
!     WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
!     P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
!     X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.
! ----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!     --------------------
!     --------------------
!     ALOG10 = LN(10)
!     RT2PIN = 1/SQRT(2*PI)
!     RTPI   = SQRT(PI)
!     --------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, x
        REAL (dpkind), INTENT (OUT) :: ans, qans
        INTEGER, INTENT (IN) :: ind
! ..
! .. Local Arrays ..
        REAL (dpkind) :: c1(8), wk(20)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a2n, a2nm1, acc, am0, amn, an, an0, apn, b2n, &
          b2nm1, c, cma, e, e0, g, h, j, l, r, rta, rtx, s, sum, t, t1, &
          tol, twoa, u, w, x0, y, z
        INTEGER :: i, iop, m, max, n, n1
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, DMAX1, EPSILON, EXP, INT, LOG, REAL, SQRT
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: acc0(3) = (/ 5.E-15_dpkind, &
          5.E-7_dpkind, 5.E-4_dpkind/)
        REAL (dpkind), PARAMETER :: big(3) = (/ 20.0E0_dpkind, &
          14.0E0_dpkind, 10.0E0_dpkind/)
        REAL (dpkind), PARAMETER :: d1(13) = (/ - &
          .185185185185185E-02_dpkind, -.347222222222222E-02_dpkind, &
          .264550264550265E-02_dpkind, -.990226337448560E-03_dpkind, &
          .205761316872428E-03_dpkind, -.401877572016461E-06_dpkind, &
          -.180985503344900E-04_dpkind, .764916091608111E-05_dpkind, &
          -.161209008945634E-05_dpkind, .464712780280743E-08_dpkind, &
          .137863344691572E-06_dpkind, -.575254560351770E-07_dpkind, &
          .119516285997781E-07_dpkind/)
        REAL (dpkind), PARAMETER :: d2(11) = (/ &
          .413359788359788E-02_dpkind, -.268132716049383E-02_dpkind, &
          .771604938271605E-03_dpkind, .200938786008230E-05_dpkind, &
          -.107366532263652E-03_dpkind, .529234488291201E-04_dpkind, &
          -.127606351886187E-04_dpkind, .342357873409614E-07_dpkind, &
          .137219573090629E-05_dpkind, -.629899213838006E-06_dpkind, &
          .142806142060642E-06_dpkind/)
        REAL (dpkind), PARAMETER :: d3(9) = (/ &
          .649434156378601E-03_dpkind, .229472093621399E-03_dpkind, &
          -.469189494395256E-03_dpkind, .267720632062839E-03_dpkind, &
          -.756180167188398E-04_dpkind, -.239650511386730E-06_dpkind, &
          .110826541153473E-04_dpkind, -.567495282699160E-05_dpkind, &
          .142309007324359E-05_dpkind/)
        REAL (dpkind), PARAMETER :: d4(7) = (/ - &
          .861888290916712E-03_dpkind, .784039221720067E-03_dpkind, &
          -.299072480303190E-03_dpkind, -.146384525788434E-05_dpkind, &
          .664149821546512E-04_dpkind, -.396836504717943E-04_dpkind, &
          .113757269706784E-04_dpkind/)
        REAL (dpkind), PARAMETER :: d5(5) = (/ - &
          .336798553366358E-03_dpkind, -.697281375836586E-04_dpkind, &
          .277275324495939E-03_dpkind, -.199325705161888E-03_dpkind, &
          .679778047793721E-04_dpkind/)
        REAL (dpkind), PARAMETER :: d6(3) = (/ &
          .531307936463992E-03_dpkind, -.592166437353694E-03_dpkind, &
          .270878209671804E-03_dpkind/)
        REAL (dpkind), PARAMETER :: d7(4) = (/ 105.0E0_dpkind, &
          3.5E0_dpkind, -one, 0.75E0_dpkind/)
        REAL (dpkind), PARAMETER :: e00(3) = (/ .25E-3_dpkind, &
          .25E-1_dpkind, .14E0_dpkind/)
        REAL (dpkind), PARAMETER :: x00(3) = (/ 31.0E0_dpkind, &
          17.0E0_dpkind, 9.7E0_dpkind/)
        REAL (dpkind), PARAMETER :: alog10 = 2.30258509299405E0_dpkind
        REAL (dpkind), PARAMETER :: d70 = .344367606892378E-03_dpkind
        REAL (dpkind), PARAMETER :: rt2pin = .398942280401433E0_dpkind
        REAL (dpkind), PARAMETER :: rtpi = 1.77245385090552E0_dpkind
        REAL (dpkind), PARAMETER :: third = .333333333333333E0_dpkind
        REAL (dpkind), PARAMETER :: d0(14) = (/ -third, &
          .833333333333333E-01_dpkind, -.148148148148148E-01_dpkind, &
          .115740740740741E-02_dpkind, .352733686067019E-03_dpkind, &
          -.178755144032922E-03_dpkind, .391926317852244E-04_dpkind, &
          -.218544851067999E-05_dpkind, -.185406221071516E-05_dpkind, &
          .829671134095309E-06_dpkind, -.176659527368261E-06_dpkind, &
          .670785354340150E-08_dpkind, .102618097842403E-07_dpkind, &
          -.438203601845335E-08_dpkind/)
! ..
!     --------------------
!     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
        e = EPSILON(one)
!     --------------------

        IF (a<zero .OR. x<zero) GO TO 350

        IF (a==zero .AND. x==zero) GO TO 350

        IF (a*x==zero) GO TO 340

        iop = ind + 1
        c1(8) = d70

        IF (iop/=1 .AND. iop/=2) iop = 3

        acc = DMAX1(acc0(iop),e)

        e0 = e00(iop)
        x0 = x00(iop)

!            SELECT THE APPROPRIATE ALGORITHM

        IF (a>=one) GO TO 10
        IF (a==half) GO TO 320
        IF (x<1.1E0_dpkind) GO TO 120

        t1 = a*LOG(x) - x
        u = a*EXP(t1)

        IF (u==zero) GO TO 310

        r = u*(one+gam1(a))
        GO TO 200

10      CONTINUE
        IF (a>=big(iop)) GO TO 30

        IF (a>x .OR. x>=x0) GO TO 20

        twoa = a + a
        m = INT(twoa)

        IF (twoa/=REAL(m,kind=dpkind)) GO TO 20

        i = m/2

        IF (a==REAL(i,kind=dpkind)) GO TO 170
        GO TO 180

20      CONTINUE
        t1 = a*LOG(x) - x
        r = EXP(t1)/gamma(a)
        GO TO 40

30      CONTINUE
        l = x/a

        IF (l==zero) GO TO 300

        s = half + (half-l)
        z = rlog(l)

        IF (z>=700.0E0_dpkind/a) GO TO 330

        y = a*z
        rta = SQRT(a)

        IF (ABS(s)<=e0/rta) GO TO 260

        IF (ABS(s)<=0.4E0_dpkind) GO TO 220

        t = (one/a)**2

        t1 = evaluate_polynomial(d7,t)/(a*1260.0E0_dpkind)

        t1 = t1 - y
        r = rt2pin*rta*EXP(t1)

40      CONTINUE
        IF (r==zero) GO TO 340
        IF (x<=DMAX1(a,alog10)) GO TO 50
        IF (x<x0) GO TO 200
        GO TO 80

!                 TAYLOR SERIES FOR P/R

50      CONTINUE
        apn = a + one
        t = x/apn
        wk(1) = t

        DO n = 2, 20
          apn = apn + one
          t = t*(x/apn)
          IF (t<=1.E-3_dpkind) GO TO 60
          wk(n) = t
        END DO

        n = 20

60      CONTINUE
        sum = t
        tol = half*acc

70      CONTINUE
        apn = apn + one
        t = t*(x/apn)
        sum = sum + t

        IF (t>tol) GO TO 70

        max = n - 1
        DO m = 1, max
          n = n - 1
          sum = sum + wk(n)
        END DO

        ans = (r/a)*(one+sum)
        qans = half + (half-ans)

        RETURN

!                 ASYMPTOTIC EXPANSION

80      CONTINUE
        amn = a - one
        t = amn/x
        wk(1) = t

        DO n = 2, 20
          amn = amn - one
          t = t*(amn/x)
          IF (ABS(t)<=1.E-3_dpkind) GO TO 90
          wk(n) = t
        END DO
        n = 20

90      CONTINUE
        sum = t

100     CONTINUE
        IF (ABS(t)<=acc) GO TO 110

        amn = amn - one
        t = t*(amn/x)
        sum = sum + t
        GO TO 100

110     CONTINUE
        max = n - 1
        DO m = 1, max
          n = n - 1
          sum = sum + wk(n)
        END DO

        qans = (r/x)*(one+sum)
        ans = half + (half-qans)
        RETURN

!             TAYLOR SERIES FOR P(A,X)/X**A

120     CONTINUE
        an = three
        c = x
        sum = x/(a+three)
        tol = three*acc/(a+one)

130     CONTINUE
        an = an + one
        c = -c*(x/an)
        t = c/(a+an)
        sum = sum + t

        IF (ABS(t)>tol) GO TO 130

        j = a*x*((sum/six-half/(a+two))*x+one/(a+one))

        z = a*LOG(x)
        h = gam1(a)
        g = one + h

        IF (x<fourth) GO TO 140

        IF (a<x/2.59E0_dpkind) GO TO 160

        GO TO 150

140     CONTINUE
        IF (z>(-.13394E0_dpkind)) GO TO 160

150     CONTINUE
        w = EXP(z)
        ans = w*g*(half+(half-j))
        qans = half + (half-ans)

        RETURN

160     CONTINUE
        l = rexp(z)
        w = half + (half+l)
        qans = (w*j-l)*g - h

        IF (qans<zero) GO TO 310

        ans = half + (half-qans)

        RETURN

!             FINITE SUMS FOR Q WHEN A .GE. 1
!                 AND 2*A IS AN INTEGER

170     CONTINUE
        sum = EXP(-x)
        t = sum
        n = 1
        c = zero
        GO TO 190

180     CONTINUE
        rtx = SQRT(x)

        sum = erfc1(0,rtx)

        t = EXP(-x)/(rtpi*rtx)
        n = 0
        c = -half

190     CONTINUE
        n1 = n
        IF (i-n1>0) THEN
          DO n = n1, i - 1
            c = c + one
            t = (x*t)/c
            sum = sum + t
          END DO
        END IF

        qans = sum
        ans = half + (half-qans)
        RETURN

!              CONTINUED FRACTION EXPANSION

200     CONTINUE
        tol = DMAX1(five*e,acc)
        a2nm1 = one
        a2n = one
        b2nm1 = x
        b2n = x + (one-a)
        c = one
210     CONTINUE
        a2nm1 = x*a2n + c*a2nm1
        b2nm1 = x*b2n + c*b2nm1
        am0 = a2nm1/b2nm1
        c = c + one
        cma = c - a
        a2n = a2nm1 + cma*a2n
        b2n = b2nm1 + cma*b2n
        an0 = a2n/b2n
        IF (ABS(an0-am0)>=tol*an0) GO TO 210

        qans = r*an0
        ans = half + (half-qans)
        RETURN

!                GENERAL TEMME EXPANSION

220     CONTINUE
        IF (ABS(s)<=two*e .AND. a*e*e>3.28E-3_dpkind) GO TO 350
        c = EXP(-y)
        w = half*erfc1(1,SQRT(y))
        u = one/a
        z = SQRT(z+z)

        IF (l<one) z = -z
        IF (iop-2>0) GO TO 240
        IF (iop-2==0) GO TO 230

        IF (ABS(s)<=thousandth) GO TO 270

        c1(1) = evaluate_polynomial(d0,z)

        c1(2) = evaluate_polynomial(d1,z)

        c1(3) = evaluate_polynomial(d2,z)

        c1(4) = evaluate_polynomial(d3,z)

        c1(5) = evaluate_polynomial(d4,z)

        c1(6) = evaluate_polynomial(d5,z)

        c1(7) = evaluate_polynomial(d6,z)

        t = evaluate_polynomial(c1,u)

        GO TO 250

230     CONTINUE
        c1(1) = evaluate_polynomial(d0(1:7),z)

        c1(2) = evaluate_polynomial(d1(1:5),z)

        c1(3) = d2(2)*z + d2(1)

        t = evaluate_polynomial(c1(1:3),u)

        GO TO 250

240     CONTINUE
        t = evaluate_polynomial(d0(1:4),z)

250     CONTINUE
        IF (l>=one) THEN
          qans = c*(w+rt2pin*t/rta)
          ans = half + (half-qans)
          RETURN

        END IF

        ans = c*(w-rt2pin*t/rta)
        qans = half + (half-ans)
        RETURN

!               TEMME EXPANSION FOR L = 1

260     CONTINUE
        IF (a*e*e>3.28E-3_dpkind) GO TO 350

        c = half + (half-y)
        w = (half-SQRT(y)*(half+(half-y/three))/rtpi)/c
        u = one/a
        z = SQRT(z+z)

        IF (l<one) z = -z
        IF (iop-2>0) GO TO 290
        IF (iop-2==0) GO TO 280

270     CONTINUE
        c1(1) = evaluate_polynomial(d0(1:8),z)

        c1(2) = evaluate_polynomial(d1(1:7),z)

        c1(3) = evaluate_polynomial(d2(1:5),z)

        c1(4) = evaluate_polynomial(d3(1:5),z)

        c1(5) = evaluate_polynomial(d4(1:3),z)

        c1(6) = evaluate_polynomial(d5(1:3),z)

        c1(7) = d6(2)*z + d6(1)

        t = evaluate_polynomial(c1,u)
        GO TO 250

280     CONTINUE
        c1(1) = (d0(3)*z+d0(2))*z - d0(1)
        c1(2) = d1(2)*z + d1(1)

        t = (d2(1)*u+c1(2))*u + c1(1)
        GO TO 250

290     CONTINUE
        t = d0(2)*z - d0(1)
        GO TO 250

!                     SPECIAL CASES

300     CONTINUE
        ans = zero
        qans = one
        RETURN

310     CONTINUE
        ans = one
        qans = zero
        RETURN

320     CONTINUE
        IF (x<fourth) THEN
          ans = erf(SQRT(x))
          qans = half + (half-ans)
          RETURN

        END IF

        qans = erfc1(0,SQRT(x))
        ans = half + (half-qans)
        RETURN

330     CONTINUE
        IF (ABS(s)<=two*e) GO TO 350
340     CONTINUE
        IF (x<=a) GO TO 300
        GO TO 310

!                     ERROR RETURN

350     CONTINUE
        ans = two

        RETURN

      END SUBROUTINE gratio

!*********************************************************************

      FUNCTION gsumln(a,b)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: gsumln
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: a, b
! ..
! .. Local Scalars ..
        REAL (dpkind) :: x
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG, REAL
! ..
        x = REAL(a,kind=dpkind) + REAL(b,kind=dpkind) - 2._dpkind

        IF (x<=0.25_dpkind) THEN
          gsumln = gamln1(one+x)
          RETURN

        END IF

        IF (x<=1.25_dpkind) THEN
          gsumln = gamln1(x) + alnrel(x)
          RETURN

        END IF

        gsumln = gamln1(x-one) + LOG(x*(one+x))

      END FUNCTION gsumln

!*********************************************************************

      FUNCTION log_beta(a0,b0)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
! .. Function Return Value ..
        REAL (dpkind) :: log_beta
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a0, b0
! ..
! .. Local Scalars ..
        REAL (dpkind) :: a, b, c, h, u, v, w, z
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG, MAX, MIN
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: e = .918938533204673_dpkind
! ..
!--------------------------
        a = MIN(a0,b0)
        b = MAX(a0,b0)

        IF (a>=eight) GO TO 50
        IF (a>=one) GO TO 10
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .LT. 1
!-----------------------------------------------------------------------
        IF (b<eight) THEN
          log_beta = gamln(a) + (gamln(b)-gamln(a+b))
          RETURN

        END IF

        log_beta = gamln(a) + algdiv(a,b)
        RETURN
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 .LE. A .LT. 8
!-----------------------------------------------------------------------
10      CONTINUE
        IF (a>two) GO TO 20
        IF (b<=two) THEN
          log_beta = gamln(a) + gamln(b) - gsumln(a,b)
          RETURN

        END IF

        w = zero
        IF (b<eight) GO TO 30

        log_beta = gamln(a) + algdiv(a,b)
        RETURN

!                REDUCTION OF A WHEN B .LE. 1000

20      CONTINUE
        IF (b>thousand) GO TO 40
        n = a - one
        w = one

        DO i = 1, n
          a = a - one
          h = a/b
          w = w*(h/(one+h))
        END DO

        w = LOG(w)

        IF (b>=eight) THEN
          log_beta = w + gamln(a) + algdiv(a,b)
          RETURN

!                 REDUCTION OF B WHEN B .LT. 8

        END IF

30      CONTINUE
        n = b - one
        z = one

        DO i = 1, n
          b = b - one
          z = z*(b/(a+b))
        END DO

        log_beta = w + LOG(z) + (gamln(a)+(gamln(b)-gsumln(a,b)))
        RETURN

!                REDUCTION OF A WHEN B .GT. 1000

40      CONTINUE
        n = a - one
        w = one
        DO i = 1, n
          a = a - one
          w = w*(a/(one+a/b))
        END DO

        log_beta = (LOG(w)-n*LOG(b)) + (gamln(a)+algdiv(a,b))
        RETURN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .GE. 8
!-----------------------------------------------------------------------
50      CONTINUE
        w = bcorr(a,b)
        h = a/b
        c = h/(one+h)
        u = -(a-half)*LOG(c)
        v = b*alnrel(h)

        IF (u>v) THEN
          log_beta = ((((-half*LOG(b))+e)+w)-v) - u
          RETURN

        END IF

        log_beta = ((((-half*LOG(b))+e)+w)-u) - v

      END FUNCTION log_beta

!*********************************************************************

      FUNCTION log_bicoef(k,n)
! .. Function Return Value ..
        REAL (dpkind) :: log_bicoef
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: k, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG
! ..
        log_bicoef = -log_beta(k+one,n-k+one) - LOG(n+one)

      END FUNCTION log_bicoef

!*********************************************************************

      FUNCTION log_gamma(a)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
!--------------------------
! .. Function Return Value ..
        REAL (dpkind) :: log_gamma
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a
! ..
! .. Local Scalars ..
        REAL (dpkind) :: t, w
        INTEGER :: i, n
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: c(6) = (/ .833333333333333E-01_dpkind &
          , -.277777777760991E-02_dpkind, .793650666825390E-03_dpkind, &
          -.595202931351870E-03_dpkind, .837308034031215E-03_dpkind, &
          -.165322962780713E-02_dpkind/)
        REAL (dpkind), PARAMETER :: d = .418938533204673_dpkind
! ..
!-----------------------------------------------------------------------
        IF (a<=0.8_dpkind) THEN
          log_gamma = gamln1(a) - LOG(a)

        ELSE IF (a<=2.25_dpkind) THEN
          t = (a-half) - half
          log_gamma = gamln1(t)

        ELSE IF (a<10.0_dpkind) THEN
          n = a - 1.25_dpkind
          t = a
          w = one

          DO i = 1, n
            t = t - one
            w = t*w
          END DO

          log_gamma = gamln1(t-one) + LOG(w)

        ELSE

          t = (one/a)**2
          w = evaluate_polynomial(c,t)/a

          log_gamma = (d+w) + (a-half)*(LOG(a)-one)

        END IF

      END FUNCTION log_gamma

!*********************************************************************

      FUNCTION psi(xx)
!---------------------------------------------------------------------
!                 EVALUATION OF THE DIGAMMA FUNCTION
!                           -----------
!     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
!     BE COMPUTED.
!     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
!     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
!     CODY, STRECOK AND THACHER.
!---------------------------------------------------------------------
!     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
!     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
!     A.H. MORRIS (NSWC).
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     PIOV4 = PI/4
!     DX0 = ZERO OF PSI TO EXTENDED PRECISION
!---------------------------------------------------------------------
!     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
!     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0
!---------------------------------------------------------------------
!     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
!     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0
!---------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: psi
! ..
! .. Scalar Arguments ..
        REAL (dpkind) :: xx
! ..
! .. Local Scalars ..
        REAL (dpkind) :: aug, den, sgn, w, x, xmax1, xmx0, xsmall, z
        INTEGER :: m, n, nq
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, COS, EPSILON, HUGE, INT, LOG, MIN, REAL, SIN
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: dx0 = &
          1.461632144968362341262659542325721325_dpkind
        REAL (dpkind), PARAMETER :: piov4 = .785398163397448_dpkind
        REAL (dpkind), PARAMETER :: p1(7) = (/ &
          .130560269827897E+04_dpkind, .413810161269013E+04_dpkind, &
          .363351846806499E+04_dpkind, .118645200713425E+04_dpkind, &
          .142441585084029E+03_dpkind, .477762828042627E+01_dpkind, &
          .895385022981970E-02_dpkind/)
        REAL (dpkind), PARAMETER :: p2(5) = (/ zero, &
          -.648157123766197E+00_dpkind, -.448616543918019E+01_dpkind, &
          -.701677227766759E+01_dpkind, -.212940445131011E+01_dpkind/)
        REAL (dpkind), PARAMETER :: q1(7) = (/ &
          .691091682714533E-05_dpkind, .190831076596300E+04_dpkind, &
          .364127349079381E+04_dpkind, .221000799247830E+04_dpkind, &
          .520752771467162E+03_dpkind, .448452573429826E+02_dpkind, one/)
        REAL (dpkind), PARAMETER :: q2(5) = (/ &
          .777788548522962E+01_dpkind, .546117738103215E+02_dpkind, &
          .892920700481861E+02_dpkind, .322703493791143E+02_dpkind, one/)
! ..
!---------------------------------------------------------------------
!     MACHINE DEPENDENT CONSTANTS ...
!        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
!                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
!                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
!                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
!                 PSI MAY BE REPRESENTED AS ALOG(X).
!        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
!                 MAY BE REPRESENTED BY 1/X.
!---------------------------------------------------------------------
        xmax1 = HUGE(1)
        xmax1 = MIN(xmax1,one/EPSILON(one))
        xsmall = 1.E-9_dpkind
!---------------------------------------------------------------------
        x = xx
        aug = zero

        IF (x>=half) GO TO 40

!---------------------------------------------------------------------
!     X .LT. 0.5,  USE REFLECTION FORMULA
!     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
        IF (ABS(x)>xsmall) GO TO 10
        IF (x==zero) GO TO 50

!---------------------------------------------------------------------
!     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
!     FOR  PI*COTAN(PI*X)
!---------------------------------------------------------------------
        aug = -one/x
        GO TO 30

!---------------------------------------------------------------------
!     REDUCTION OF ARGUMENT FOR COTAN
!---------------------------------------------------------------------
10      CONTINUE
        w = -x
        sgn = piov4
        IF (w<=zero) THEN
          w = -w
          sgn = -sgn
!---------------------------------------------------------------------
!     MAKE AN ERROR EXIT IF X .LE. -XMAX1
!---------------------------------------------------------------------
        END IF

        IF (w>=xmax1) GO TO 50

        nq = INT(w)
        w = w - REAL(nq,kind=dpkind)
        nq = INT(w*4.0_dpkind)
        w = 4.0_dpkind*(w-REAL(nq,kind=dpkind)*.25_dpkind)

!---------------------------------------------------------------------
!     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
!     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
!     QUADRANT AND DETERMINE SIGN
!---------------------------------------------------------------------
        n = nq/2
        IF (n+n/=nq) w = one - w

        z = piov4*w
        m = n/2

        IF (m+m/=n) sgn = -sgn

!---------------------------------------------------------------------
!     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
!---------------------------------------------------------------------
        n = (nq+1)/2
        m = n/2
        m = m + m
        IF (m/=n) GO TO 20

!---------------------------------------------------------------------
!     CHECK FOR SINGULARITY
!---------------------------------------------------------------------
        IF (z==zero) GO TO 50
!---------------------------------------------------------------------
!     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
!     SIN/COS AS A SUBSTITUTE FOR TAN
!---------------------------------------------------------------------
        aug = sgn*((COS(z)/SIN(z))*four)
        GO TO 30

20      CONTINUE
        aug = sgn*((SIN(z)/COS(z))*four)

30      CONTINUE
        x = one - x

40      CONTINUE
        IF (x<=three) THEN
!---------------------------------------------------------------------
!     0.5 .LE. X .LE. 3.0
!---------------------------------------------------------------------

          den = evaluate_polynomial(p1,x)/evaluate_polynomial(q1,x)

          xmx0 = x - dx0
          psi = den*xmx0 + aug
          RETURN

!---------------------------------------------------------------------
!     IF X .GE. XMAX1, PSI = LN(X)
!---------------------------------------------------------------------
        END IF

        IF (x<xmax1) THEN
!---------------------------------------------------------------------
!     3.0 .LT. X .LT. XMAX1
!---------------------------------------------------------------------
          w = one/(x*x)

          aug = evaluate_polynomial(p2,w)/evaluate_polynomial(q2,w) - &
            half/x + aug
        END IF

        psi = aug + LOG(x)

        RETURN

!---------------------------------------------------------------------
!     ERROR RETURN
!---------------------------------------------------------------------
50      CONTINUE

        psi = zero

      END FUNCTION psi

!*********************************************************************

      FUNCTION rcomp(a,x)
!     -------------------
!     EVALUATION OF EXP(-X)*X**A/GAMMA(A)
!     -------------------
!     RT2PIN = 1/SQRT(2*PI)
!     -------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: a, x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: t, t1, u
! ..
! .. Intrinsic Functions ..
        INTRINSIC EXP, LOG, SQRT
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: c(4) = (/ -105.0E0_dpkind, &
          3.5E0_dpkind, -one, 0.75E0_dpkind/)
        REAL (dpkind), PARAMETER :: rt2pin = .398942280401433E0_dpkind
! ..
! .. Function Return Value ..
        REAL (dpkind) :: rcomp
! ..
        rcomp = zero

        IF (x<=zero) RETURN

        IF (a<20.0E0_dpkind) THEN
          t = a*LOG(x) - x

          IF (a<one) THEN
            rcomp = (a*EXP(t))*(one+gam1(a))
          ELSE
            rcomp = EXP(t)/gamma(a)
          END IF

        ELSE
          u = x/a

          t = (one/a)**2
          t1 = evaluate_polynomial(c,t)/(a*1260.0E0_dpkind) - a*rlog(u)

          rcomp = rt2pin*SQRT(a)*EXP(t1)
        END IF

      END FUNCTION rcomp

!*********************************************************************

      FUNCTION rexp(x)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
! .. Function Return Value ..
        REAL (dpkind) :: rexp
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: w
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EXP
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: p1 = .914041914819518E-09_dpkind
        REAL (dpkind), PARAMETER :: p2 = .238082361044469E-01_dpkind
        REAL (dpkind), PARAMETER :: q(5) = (/ one, &
          -.499999999085958E+00_dpkind, .107141568980644E+00_dpkind, &
          -.119041179760821E-01_dpkind, .595130811860248E-03_dpkind/)
! ..
        IF (ABS(x)<=0.15_dpkind) THEN
          rexp = x*(((p2*x+p1)*x+one)/evaluate_polynomial(q,x))
        ELSE
          w = EXP(x)

          IF (x<=zero) THEN
            rexp = (w-half) - half
          ELSE
            rexp = w*(half+(half-one/w))
          END IF
        END IF
      END FUNCTION rexp

!*********************************************************************

      FUNCTION rlog(x)
!     -------------------
!     COMPUTATION OF  X - 1 - LN(X)
!     -------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: r, t, u, w, w1
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: a = .566749439387324E-01_dpkind
        REAL (dpkind), PARAMETER :: b = .456512608815524E-01_dpkind
        REAL (dpkind), PARAMETER :: p(3) = (/ .333333333333333E+00_dpkind &
          , -.224696413112536E+00_dpkind, .620886815375787E-02_dpkind/)
        REAL (dpkind), PARAMETER :: q(3) = (/ one, &
          -.127408923933623E+01_dpkind, .354508718369557E+00_dpkind/)
! ..
! .. Function Return Value ..
        REAL (dpkind) :: rlog
! ..
!     -------------------
        IF (x>=0.61E0_dpkind .AND. x<=1.57E0_dpkind) THEN
          IF (x<0.82E0_dpkind) GO TO 10
          IF (x>1.18E0_dpkind) GO TO 20

!              ARGUMENT REDUCTION

          u = (x-half) - half
          w1 = zero
          GO TO 30

10        CONTINUE
          u = x - 0.7E0_dpkind
          u = u/0.7E0_dpkind
          w1 = a - u*0.3E0_dpkind
          GO TO 30

20        CONTINUE
          u = 0.75E0_dpkind*x - one
          w1 = b + u/three

!               SERIES EXPANSION

30        CONTINUE
          r = u/(u+two)
          t = r*r
          w = evaluate_polynomial(p,t)/evaluate_polynomial(q,t)

          rlog = two*t*(one/(one-r)-r*w) + w1

          RETURN

        END IF

        r = (x-half) - half

        rlog = r - LOG(x)

        RETURN

      END FUNCTION rlog

!*********************************************************************

      FUNCTION rlog1(x)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - LN(1 + X)
!-----------------------------------------------------------------------
!------------------------
! .. Function Return Value ..
        REAL (dpkind) :: rlog1
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x
! ..
! .. Local Scalars ..
        REAL (dpkind) :: h, r, t, w, w1
! ..
! .. Intrinsic Functions ..
        INTRINSIC LOG
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: a = .566749439387324E-01_dpkind
        REAL (dpkind), PARAMETER :: b = .456512608815524E-01_dpkind
        REAL (dpkind), PARAMETER :: p(3) = (/ .333333333333333E+00_dpkind &
          , -.224696413112536E+00_dpkind, .620886815375787E-02_dpkind/)
        REAL (dpkind), PARAMETER :: q(3) = (/ one, &
          -.127408923933623E+01_dpkind, .354508718369557E+00_dpkind/)
! ..
!------------------------
        IF (x<=-one) THEN
          rlog1 = -one
        ELSE IF (x>=(-0.39_dpkind) .AND. x<=0.57_dpkind) THEN
          IF (x<-0.18_dpkind) THEN
            h = x + 0.3_dpkind
            h = h/0.7_dpkind
            w1 = a - h*0.3_dpkind
          ELSE IF (x>0.18_dpkind) THEN
            h = 0.75_dpkind*x - 0.25_dpkind
            w1 = b + h/three
          ELSE
! ARGUMENT REDUCTION

            h = x
            w1 = zero
          END IF

! SERIES EXPANSION

          r = h/(h+two)
          t = r*r

          w = evaluate_polynomial(p,t)/evaluate_polynomial(q,t)

          rlog1 = two*t*(one/(one-r)-r*w) + w1

        ELSE
          w = (x+half) + half
          rlog1 = x - LOG(w)
        END IF

      END FUNCTION rlog1

!*********************************************************************

    END MODULE biomath_mathlib_mod
