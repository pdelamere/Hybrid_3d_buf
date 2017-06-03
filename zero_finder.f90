    MODULE zero_finder
!!!**********************************************************************
!!!   These variables provide the  status of the zero_finder  routine either
!!!   currently or most recently in use.
!!!
!!!   zf_current_status:
!!!      1 - A reverse communication zero-finder is in use
!!!      0 - The most recently used zero-finder completed successfully
!!!     -1 - The most recently used zero-finder completed unsuccessfully
!!!          (NO roots between the bounds)
!!!
!!!   zf_crash_left:
!!!           .TRUE. iff the zero finder tried to exceed the left bound
!!!                  of the interval.  Root is set to the left bound.
!!!                  
!!!           .FALSE. iff the zero finder tried to exceed the right bound
!!!                  of the interval.  Root is set to the right bound.
!!!
!!!   zf_crash_hi:
!!!           .TRUE. iff f(x) > 0 at BOTH left and right bounds of the
!!!                  interval.
!!!           .FALSE. iff f(x) < 0 at BOTH left and right bounds of the
!!!                  interval.
!!!
!!!   NOTE:
!!!   The two logical status variables described above are defined ONLY if:
!!!             zf_current_status == -1.
!!!
!!!**********************************************************************
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: final_zf_state, interval_zf, rc_interval_zf, rc_step_zf, &
        set_zero_finder, step_zf, zf_bound_high, zf_bound_low, &
        zf_crash_hi, zf_crash_left, zf_current_status, zf_locals
! ..
! .. Parameters ..
      REAL (dpkind), PARAMETER :: converge = 1.0E-6_dpkind
      REAL (dpkind), PARAMETER :: default_abs_tol = converge
      REAL (dpkind), PARAMETER :: default_rel_step = 1.0E0_dpkind
      REAL (dpkind), PARAMETER :: default_rel_tol = converge
      REAL (dpkind), PARAMETER :: default_step_multiplier = 2.0_dpkind
      REAL (dpkind), PARAMETER :: infinity = 1.0E35_dpkind
      REAL (dpkind), PARAMETER :: interval_init_f = 1.0E5_dpkind
      REAL (dpkind), PARAMETER :: interval_lambda = 0.7_dpkind
      REAL (dpkind), PARAMETER :: interval_mu = 0.5_dpkind
      REAL (dpkind), PARAMETER :: small_step = 1.0E-4_dpkind
      REAL (dpkind), PARAMETER :: default_abs_step = small_step
      REAL (dpkind), PARAMETER :: default_hi_limit = infinity
      REAL (dpkind), PARAMETER :: default_low_limit = -infinity
      INTEGER, PARAMETER :: common_tolerances_are_set = -123456, &
        interval_neps = 100
! ..
! .. Derived Type Declarations ..
      TYPE :: zf_locals
        PRIVATE
        REAL (dpkind) :: step_absstp
        REAL (dpkind) :: common_abstol
        REAL (dpkind) :: common_big
        REAL (dpkind) :: step_fbig
        REAL (dpkind) :: step_fsmall
        REAL (dpkind) :: step_relstp
        REAL (dpkind) :: common_reltol
        REAL (dpkind) :: common_small
        REAL (dpkind) :: step_step
        REAL (dpkind) :: step_stpmul
        REAL (dpkind) :: step_xhi
        REAL (dpkind) :: step_xlb
        REAL (dpkind) :: step_xlo
        REAL (dpkind) :: step_xsave
        REAL (dpkind) :: step_xub
        REAL (dpkind) :: step_yy
        INTEGER :: step_from
        INTEGER :: common_tolerances_set = 0
        LOGICAL :: step_qbdd
        LOGICAL :: step_qcond
        LOGICAL :: step_qincr
        LOGICAL :: step_qlim
        LOGICAL :: step_qok
        LOGICAL :: step_qup
        LOGICAL :: common_called_from_step
        REAL (dpkind) :: interval_a
        REAL (dpkind) :: interval_a0
        REAL (dpkind) :: interval_b
        REAL (dpkind) :: interval_b0
        REAL (dpkind) :: interval_c
        REAL (dpkind) :: interval_d
        REAL (dpkind) :: interval_e
        REAL (dpkind) :: interval_u
        REAL (dpkind) :: interval_fa
        REAL (dpkind) :: interval_fb
        REAL (dpkind) :: interval_fc
        REAL (dpkind) :: interval_fd
        REAL (dpkind) :: interval_fe
        REAL (dpkind) :: interval_fu
        REAL (dpkind) :: interval_tol0
        REAL (dpkind) :: interval_tol
        REAL (dpkind) :: interval_eps
        REAL (dpkind) :: interval_lenab
        INTEGER :: interval_from
        INTEGER :: interval_itnum
        INTEGER :: common_current_status
        LOGICAL :: common_crash_hi
        LOGICAL :: common_crash_left
        REAL (dpkind) :: common_zf_bound_low
        REAL (dpkind) :: common_zf_bound_high
      END TYPE zf_locals
! ..
! .. Local Scalars ..
      REAL (dpkind) :: zf_bound_high, zf_bound_low
      INTEGER :: zf_current_status = 1
      LOGICAL :: zf_crash_hi, zf_crash_left
! ..
! .. Dependents ..
      TYPE (zf_locals), SAVE :: default_locals_zf
! ..
    CONTAINS

!*********************************************************************

      SUBROUTINE final_zf_state(local,status,crash_left,crash_hi, &
          left_end,right_end)
! .. Structure Arguments ..
        TYPE (zf_locals), INTENT (IN) :: local
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: left_end, right_end
        INTEGER, INTENT (OUT) :: status
        LOGICAL, INTENT (OUT) :: crash_hi, crash_left
! ..
        status = local%common_current_status

        IF (status==0) RETURN

        crash_left = local%common_crash_left
        crash_hi = local%common_crash_hi

        left_end = local%common_small
        right_end = local%common_big

        RETURN

      END SUBROUTINE final_zf_state

!*********************************************************************

      SUBROUTINE interval_zf(f,y,answer,status,local)
! .. Structure Arguments ..
        TYPE (zf_locals), OPTIONAL :: local
! ..
! .. Non-Generic Interface Blocks ..
        INTERFACE
          FUNCTION f(x)
! .. Use Statements ..
            USE biomath_constants_mod
! ..
! .. Function Return Value ..
            REAL (dpkind) :: f
! ..
! .. Scalar Arguments ..
            REAL (dpkind), INTENT (IN) :: x
! ..
          END FUNCTION f
        END INTERFACE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: answer
        REAL (dpkind), INTENT (IN) :: y
        INTEGER, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: fx, x
        LOGICAL :: locals_present
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        locals_present = PRESENT(local)

        status = 0

        DO
          IF (locals_present) THEN
            CALL rc_interval_local_zf(status,x,fx,local)
          ELSE
            CALL rc_interval_local_zf(status,x,fx,default_locals_zf)
          END IF

          IF (status/=1) EXIT

          fx = f(x) - y

        END DO

        IF (status==0) answer = x

        RETURN

      END SUBROUTINE interval_zf

!*********************************************************************

      SUBROUTINE rc_interval_local_zf(status,x,fx,local)
!---------------------------------------------------------------------
! Finds either an EXACT solution or an APPROXIMATE solution of the
! equation f(x)=0 in the interval [a,b].
! At the begining of each iteration, the current enclosing interval is
! recorded as [a0,b0].
! The FIRST iteration is simply a secant step.
! Starting with the SECOND iteration, THREE STEPS are taken in EACH
! iteration.  The first two steps are either QUADRATIC interpolation
! or CUBIC inverse interpolation. The third step is a double-size
! secant step.  If the diameter of the enclosing interval obtained
! after those three steps is larger than 0.5*(b0-a0), then
! an additional bisection step will be taken.
!  A,B   <-> INPUT:  the initial interval
!            OUTPUT: the enclosing interval at the termination;
!  ROOT  <-- Solution of the equation.
! DMS Each time f(x) is needed, the routine sets status=1 and
!     returns to the caller, then comes back and continues where it
!     left off.
! The algorithm is described in the article:
! Algorithm 748: Enclosing Zeros of Continuous Functions,
! by G. E. Alefeld, F. A. Potra, YiXun Shi,
! ACM Transactions on Mathematical Software, Vol. 21, No. 3, Sep. 1995
! pages 327-344
!---------------------------------------------------------------------
! .. Intrinsic Functions ..
        INTRINSIC ABS, EPSILON
! ..
! .. Structure Arguments ..
! .. Local Arguments ..
        TYPE (zf_locals) :: local
! ..
! .. Local Scalars ..
        INTEGER :: iprof
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: fx
        REAL (dpkind), INTENT (OUT) :: x
        INTEGER, INTENT (INOUT) :: status
! ..
        IF (status==0) THEN
!!!     Make sure that local has been set
          IF (local%common_tolerances_set/=common_tolerances_are_set) &
              THEN
            STOP 'Must call set_zero_finder before a zero finding &
              &routine -- ABORT'
          END IF

! Some initializations are done

          IF (local%common_called_from_step) THEN
            local%interval_a = local%step_xlb
            local%interval_b = local%step_xub
          ELSE
            local%interval_a = local%common_small
            local%interval_b = local%common_big
          END IF

! ***********************
          local%interval_tol = zero
! ***********************

          local%interval_itnum = 0
          local%interval_e = interval_init_f
          local%interval_fe = interval_init_f
          local%interval_eps = EPSILON(one)

          IF (interval_neps>300) THEN
            local%interval_tol0 = zero
          ELSE
            local%interval_tol0 = one/(ten**interval_neps)
          END IF


! Go back to the caller and compute f(a)

          x = local%interval_a
          status = 1
          local%interval_from = 1
          RETURN
        END IF

        IF (status==1) THEN
          SELECT CASE (local%interval_from)
          CASE (1)
! Called only once to compute: fa = f(a)
            GO TO 10

          CASE (2)
! Called only once to compute: fb = f(b)
            GO TO 20

          CASE (3)
! Called only for ITNUM=1 to compute: fc = f(c)
            GO TO 40

          CASE (4)
            GO TO 60

          CASE (5)
            GO TO 70

          CASE (6)
            GO TO 80

          CASE (7)
            GO TO 90

          END SELECT
        END IF

        IF (status/=0) THEN
          WRITE (*,*) 'Value is: ', status
          STOP 'Illegal input value of status in routine rc_interval_local&
            &_zf'
        END IF

10      CONTINUE

! Function value at the LEFT bound
        local%interval_fa = fx

! Go back to the caller and compute f(b)

        x = local%interval_b
        status = 1
        local%interval_from = 2
        RETURN

20      CONTINUE

! Function value at the RIGHT bound
        local%interval_fb = fx

! Check whether there is a root in [a, b]

        IF ((local%interval_fa<zero) .AND. (local%interval_fb<zero)) THEN
! Function values at BOTH bounds are NEGATIVE.
! NO root in [a, b]

          status = -1
          zf_crash_hi = .FALSE.
          zf_crash_left = (local%interval_fa>local%interval_fb)
        ELSE IF ((local%interval_fa>zero) .AND. (local%interval_fb>zero)) &
            THEN
! Function values at BOTH bounds are POSITIVE.
! NO root in [a, b]

          status = -1
          zf_crash_hi = .TRUE.
          zf_crash_left = (local%interval_fa<local%interval_fb)
        END IF

        IF (status==-1) THEN
! NO root in [a, b]

          zf_current_status = -1
          local%common_current_status = status
          local%common_crash_left = zf_crash_left
          local%common_crash_hi = zf_crash_hi

! Set the root to one of the bounds.

          IF (zf_crash_left) THEN
            x = local%interval_a
          ELSE
            x = local%interval_b
          END IF

          RETURN
        END IF

! f(a) * f(b) <= 0 ==> there is a root x in [a, b] | f(x) = 0

! Iteration starts.
! The enclosing interval BEFORE executing the iteration is recorded
! as [a0, b0].

30      local%interval_a0 = local%interval_a
        local%interval_b0 = local%interval_b

! Length of the interval [a, b]
        local%interval_lenab = local%interval_b - local%interval_a

        local%interval_itnum = local%interval_itnum + 1

! Calculates the TERMINATION criterion.
! EXITs the procedure if the criterion is satisfied.

        CALL ftol

        IF (local%interval_lenab<=local%interval_tol) THEN
          GO TO 100
        END IF

! For the FIRST iteration, SECANT STEP is taken.

        IF (local%interval_itnum/=1) GO TO 50
        local%interval_c = local%interval_a - (local%interval_fa/(local% &
          interval_fb-local%interval_fa))*local%interval_lenab

! Call subroutine "BRACKT" to get a shrinked enclosing interval
! as well as to update the termination criterion.
! EXITs the procedure if the criterion is satisfied OR
! the EXACT SOLUTION is obtained.

        CALL adjust_c

        x = local%interval_c
        local%interval_from = 3
        status = 1
        RETURN

40      CONTINUE
        local%interval_fc = fx

        CALL bracket

        IF (has_converged()) THEN
          GO TO 100
        ELSE
          GO TO 30
        END IF

50      CONTINUE

! Starting with the SECOND iteration, in the FIRST TWO STEPS, either
! QUADRATIC INTERPOLATION is used by calling the subroutine "Newton_Quadratic"
! or THE CUBIC INVERSE INTERPOLATION is used by calling the subroutine
! "PZERO".
! In the following, if "IPROF" is not equal to 0, then the
! FOUR function values "fa", "fb", "fd", and "fe" are DISTINCT,
! and hence "PZERO" will be called.

        IF ((local%interval_itnum==2) .OR. (ifabde()==0)) THEN
          CALL newton_quadratic(2)
        ELSE
          CALL pzero

          IF (iisign(local%interval_c-local%interval_a)* &
              iisign(local%interval_c-local%interval_b)>=0) THEN
            CALL newton_quadratic(2)
          END IF
        END IF

        local%interval_e = local%interval_d
        local%interval_fe = local%interval_fd

        CALL adjust_c

        x = local%interval_c
        local%interval_from = 4
        status = 1
        RETURN

60      CONTINUE
        local%interval_fc = fx

        CALL bracket

        IF (has_converged()) GO TO 100

        IF (ifabde()==0) THEN
          CALL newton_quadratic(3)
        ELSE
          CALL pzero

          iprof = iisign(local%interval_c-local%interval_a)* &
            iisign(local%interval_c-local%interval_b)

          IF (iprof>=0) THEN
            CALL newton_quadratic(3)
          END IF
        END IF

        CALL adjust_c

        x = local%interval_c
        local%interval_from = 5
        status = 1
        RETURN

70      CONTINUE
        local%interval_fc = fx

        CALL bracket

        IF (has_converged()) GO TO 100

        local%interval_e = local%interval_d
        local%interval_fe = local%interval_fd

! Takes the DOUBLE-SIZE SECANT STEP.

        IF (ABS(local%interval_fa)<ABS(local%interval_fb)) THEN
          local%interval_u = local%interval_a
          local%interval_fu = local%interval_fa
        ELSE
          local%interval_u = local%interval_b
          local%interval_fu = local%interval_fb
        END IF

        local%interval_c = local%interval_u - two*(local%interval_fu/( &
          local%interval_fb-local%interval_fa))*local%interval_lenab

        IF (ABS(local%interval_c-local%interval_u)> &
            half*local%interval_lenab) THEN
          local%interval_c = half*(local%interval_a+local%interval_b)
        END IF

        CALL adjust_c

        x = local%interval_c
        local%interval_from = 6
        status = 1
        RETURN

80      CONTINUE
        local%interval_fc = fx

        CALL bracket

        IF (has_converged()) GO TO 100

! Determines whether an additional BISECTION STEP is needed.
! Takes it if necessary.

        IF (local%interval_lenab<interval_mu*(local%interval_b0-local% &
            interval_a0)) THEN
          GO TO 30
        END IF

        local%interval_e = local%interval_d
        local%interval_fe = local%interval_fd

        local%interval_c = half*(local%interval_a+local%interval_b)

        CALL adjust_c

        x = local%interval_c
        local%interval_from = 7
        status = 1
        RETURN

90      CONTINUE
        local%interval_fc = fx

        CALL bracket

        IF (has_converged()) THEN
          GO TO 100
        ELSE
          GO TO 30
        END IF

! Terminates the procedure and return the "root" in x

100     CONTINUE

        status = 0
        zf_current_status = 0
        zf_bound_low = local%interval_a
        zf_bound_high = local%interval_b
        local%common_current_status = 0
        local%common_zf_bound_low = zf_bound_low
        local%common_zf_bound_high = zf_bound_high
        x = half*(zf_bound_low+zf_bound_high)

        RETURN

!-----------------------------------------------------------------------

      CONTAINS

!.....................................................................

        SUBROUTINE adjust_c
! Adjust c if (b-a) is very small or if c is very close to a or b.
          local%interval_tol = interval_lambda*local%interval_tol

          IF ((local%interval_lenab)<=two*local%interval_tol) THEN
            local%interval_c = half*(local%interval_a+local%interval_b)
          ELSE IF (local%interval_c<=local%interval_a+local%interval_tol) &
              THEN
            local%interval_c = local%interval_a + local%interval_tol
          ELSE
            IF (local%interval_c>=local%interval_b-local%interval_tol) &
                THEN
              local%interval_c = local%interval_b - local%interval_tol
            END IF
          END IF

        END SUBROUTINE adjust_c

!.....................................................................

        SUBROUTINE bracket
! Given current enclosing interval [a,b] and a number c in (a,b), if
! f(c)=0 then sets the output a=c. Otherwise determines the new
! enclosing interval: [a,b]=[a,c] or [a,b]=[c,b]. Also updates the
! termination criterion corresponding to the new enclosing interval.
!  a, b     <-> INPUT: [a,b] is the current enclosing interval
!              OUTPUT: [a,b] is the NEW shrinked enclosing interval;
!  c       --> Used to determine the new enclosing interval;
!  d       <-- If the new enclosing interval is [a,c] then d=b,
!              otherwise d=a;
!  fa,fb,fd -- fa=f(a), fb=f(b), fd=f(d);
!  tol      <-> INPUT: the current termination criterion
!               OUTPUT: the updated termination criterion,
!                       according to the new enclosing interval;
! If f(c)=0, then set a=c and return. This will terminate the
! procedure in subroutine "rroot" and give the EXACT solution of
! the equation f(x)=0.
          IF (local%interval_fc==zero) THEN
            local%interval_a = local%interval_c
            local%interval_fa = zero
            local%interval_d = zero
            local%interval_fd = zero

            RETURN
          END IF

! If f(c) is not zero, then determine the new enclosing interval.

          IF ((iisign(local%interval_fa)*iisign(local%interval_fc))<0) &
              THEN
! ... Root is in (a, c)
            local%interval_d = local%interval_b
            local%interval_fd = local%interval_fb
            local%interval_b = local%interval_c
            local%interval_fb = local%interval_fc
          ELSE
! ... Root is in (c, b)
            local%interval_d = local%interval_a
            local%interval_fd = local%interval_fa
            local%interval_a = local%interval_c
            local%interval_fa = local%interval_fc
          END IF

          local%interval_lenab = local%interval_b - local%interval_a

! UPDATE the TERMINATION CRITERION according to
! the new enclosing interval.

          CALL ftol

          RETURN

        END SUBROUTINE bracket

!.....................................................................

        SUBROUTINE ftol

          IF (ABS(local%interval_fb)<=ABS(local%interval_fa)) THEN
            local%interval_tol = two*(local%interval_tol0+two*ABS(local% &
              interval_b)*local%interval_eps)
          ELSE
            local%interval_tol = two*(local%interval_tol0+two*ABS(local% &
              interval_a)*local%interval_eps)
          END IF

          RETURN

        END SUBROUTINE ftol

!.....................................................................

        FUNCTION has_converged()

! .. Function Return Value ..
          LOGICAL :: has_converged
! ..
          has_converged = (local%interval_fa==zero) .OR. &
            (local%interval_lenab<=local%interval_tol)

          RETURN

        END FUNCTION has_converged

!.....................................................................

        FUNCTION ifabde()
! .. Function Return Value ..
          INTEGER :: ifabde
! ..
          ifabde = iisign(local%interval_fa-local%interval_fb)* &
            iisign(local%interval_fa-local%interval_fd)* &
            iisign(local%interval_fa-local%interval_fe)* &
            iisign(local%interval_fb-local%interval_fd)* &
            iisign(local%interval_fb-local%interval_fe)* &
            iisign(local%interval_fd-local%interval_fe)

          RETURN

        END FUNCTION ifabde

!.....................................................................

        FUNCTION iisign(x)
!    Returns the SIGN of variable "x".
! .. Function Return Value ..
          INTEGER :: iisign
! ..
! .. Scalar Arguments ..
          REAL (dpkind), INTENT (IN) :: x
! ..
          IF (x>zero) THEN
            iisign = 1
          ELSE IF (x==zero) THEN
            iisign = 0
          ELSE
            iisign = -1
          END IF

          RETURN

        END FUNCTION iisign

!.....................................................................

        SUBROUTINE newton_quadratic(k)
! Uses K Newton steps to approximate the zero in (a,b) of the
! QUADRATIC POLYNOMIAL interpolating f(x) at a, b, and d.
! Safeguard is used to AVOID OVERFLOW.
!  a,b,d,fa,fb,fd --> d lies outside [a,b].
!                     fa=f(a), fb=f(b), and fd=f(d). f(a)f(b)<0.
!  c              <-- Approximate zero in (a,b) of the QUADRATIC POLYNOMIAL.
!  K              --> INTEGER. The number of Newton steps to take.
! .. Scalar Arguments ..
          INTEGER :: k
! ..
! .. Local Scalars ..
          REAL (dpkind) :: a0, a1, a2, pc, pdc
          INTEGER :: i, ierror
! ..
! Find the coefficients of the QUADRATIC POLYNOMIAL.
          ierror = 0

          a0 = local%interval_fa
          a1 = (local%interval_fb-local%interval_fa)/local%interval_lenab
          a2 = ((local%interval_fd-local%interval_fb)/ &
            (local%interval_d-local%interval_b)-a1)/ &
            (local%interval_d-local%interval_a)

! Safeguard to AVOID OVERFLOW.

10        IF (a2==zero .OR. ierror==1) THEN
            local%interval_c = local%interval_a - a0/a1

            RETURN
          END IF

! Determine the STARTING POINT of Newton steps.

          IF (iisign(a2)*iisign(local%interval_fa)>0) THEN
            local%interval_c = local%interval_a
          ELSE
            local%interval_c = local%interval_b
          END IF

! Start the SAFEGUARDED Newton steps.

          DO i = 1, k
            IF (ierror==0) THEN
              pc = a0 + (a1+a2*(local%interval_c-local%interval_b))*( &
                local%interval_c-local%interval_a)
              pdc = a1 + a2*((two*local%interval_c)-(local%interval_a+ &
                local%interval_b))
              IF (pdc==zero) THEN
                ierror = 1
              ELSE
                local%interval_c = local%interval_c - pc/pdc
              END IF
            END IF
          END DO

          IF (ierror==1) GO TO 10

          RETURN

        END SUBROUTINE newton_quadratic

!.....................................................................

        SUBROUTINE pzero
! Uses CUBIC INVERSE INTERPOLATION of f(x) at a, b, d, and e
! to get an approximate root of f(x).
! This procedure is a slight modification of Aitken-Neville algorithm
! for interpolation described by:
! Stoer and Bulirsch in "Introduction to Numerical Analysis",
! Springer-Verlag. New York (1980).
!  a,b,d,e,fa,fb,fd,fe -- d and e lie outside the interval [a,b].
!                         fa=f(a), fb=f(b), fd=f(d), fe=f(e).
!  c                   -- Output of the subroutine.
! NOTE: The parameters described above are members of the "local" structure
!       and are prepended by local%interval_ in the code
! .. Local Scalars ..
          REAL (dpkind) :: const, d21, d31, d32, q11, q21, q22, q31, q32, &
            q33
! ..
          q11 = (local%interval_d-local%interval_e)* &
            (local%interval_fd/(local%interval_fe-local%interval_fd))
          q21 = (local%interval_b-local%interval_d)* &
            (local%interval_fb/(local%interval_fd-local%interval_fb))

          const = -local%interval_lenab/(local%interval_fb-local% &
            interval_fa)
          q31 = const*local%interval_fa
          d21 = (local%interval_b-local%interval_d)* &
            (local%interval_fd/(local%interval_fd-local%interval_fb))
          d31 = const*local%interval_fb

          q22 = (d21-q11)*(local%interval_fb/(local%interval_fe-local% &
            interval_fb))
          const = (d31-q21)/(local%interval_fd-local%interval_fa)
          q32 = const*local%interval_fa
          d32 = const*local%interval_fd
          q33 = (d32-q22)*(local%interval_fa/(local%interval_fe-local% &
            interval_fa))

! Calculate the output c.

          local%interval_c = local%interval_a + q31 + q32 + q33

          RETURN

        END SUBROUTINE pzero

!.....................................................................

      END SUBROUTINE rc_interval_local_zf

!*********************************************************************

      SUBROUTINE rc_interval_zf(status,x,fx,local)
! .. Structure Arguments ..
        TYPE (zf_locals), OPTIONAL :: local
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: fx
        REAL (dpkind), INTENT (OUT) :: x
        INTEGER, INTENT (INOUT) :: status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        IF (PRESENT(local)) THEN
          CALL rc_interval_local_zf(status,x,fx,local)
        ELSE
          CALL rc_interval_local_zf(status,x,fx,default_locals_zf)
        END IF

        RETURN

      END SUBROUTINE rc_interval_zf

!*********************************************************************

      SUBROUTINE rc_step_local_zf(status,x,fx,local)
! .. Intrinsic Functions ..
        INTRINSIC ABS, MAX, MIN
! ..
! .. Structure Arguments ..
        TYPE (zf_locals) :: local
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: fx
        REAL (dpkind), INTENT (OUT) :: x
        INTEGER, INTENT (INOUT) :: status
! ..
        IF (status==1) THEN
          SELECT CASE (local%step_from)
          CASE (1)
            GO TO 10

          CASE (2)
            GO TO 20

          CASE (3)
            GO TO 30

          CASE (4)
            GO TO 50

          CASE (5)
            GO TO 80

          CASE (6)
            GO TO 100
          END SELECT
        END IF

        IF (status/=0) THEN
          WRITE (*,*) 'Value is: ', status
          STOP &
            'Illegal input value of status in routine rc_step_local_zf'
        END IF

        local%common_current_status = status

        IF ( .NOT. qxmon(local%common_small,x,local%common_big)) THEN
          WRITE (*,*) 'Small, x, big not monotone in rc_step_local_zf'
          WRITE (*,*) 'Value of small: ', local%common_small
          WRITE (*,*) 'Value of x: ', x
          WRITE (*,*) 'Value of big: ', local%common_big
          STOP ' x, big not monotone in rc_step_local_zf'
        END IF

        local%step_xsave = x

!     See that SMALL and BIG bound the zero and set QINCR

        x = local%common_small
        local%step_from = 1
        status = 1
        local%common_current_status = status

        RETURN

10      CONTINUE

! Function value at the LEFT bound
        local%step_fsmall = fx

        x = local%common_big
        local%step_from = 2
        status = 1
        local%common_current_status = status

        RETURN

20      CONTINUE

! Function value at the RIGHT bound
        local%step_fbig = fx

! Decide the MONOTONICITY (whether the function is) :
! INCREASING - when f(left) < f(right)
! DECREASING - when f(left) > f(right)

        local%step_qincr = (local%step_fbig>local%step_fsmall)

        IF (local%step_qincr) THEN
! The function is INCREASING
! To have a root in [left, right]
!    * f(left)  MUST be NEGATIVE
!    * f(right) MUST be POSITIVE

          IF (local%step_fsmall>zero) THEN
! BOTH functions are POSITIVE.
! Condition is NOT met and as such:
! there is NO root in [left, right]

            status = -1
            zf_crash_left = .TRUE.
            zf_crash_hi = .TRUE.
          ELSE IF (local%step_fbig<zero) THEN
! BOTH functions are NEGATIVE.
! NO root in [left, right]

            status = -1
            zf_crash_left = .FALSE.
            zf_crash_hi = .FALSE.
          END IF
        ELSE
! Function is DECREASING
! To have a root in [left, right], these conditions MUST be met:
!    * f(left) > 0
!    * f(right) < 0

          IF (local%step_fsmall<zero) THEN
! f(left) < 0, but because f(x) is DECREASING ==> f(right) < 0
! BOTH functions are NEGATIVE at the bounds.
! NO root in f(left, right)

            status = -1
            zf_crash_left = .TRUE.
            zf_crash_hi = .FALSE.
          ELSE IF (local%step_fbig>zero) THEN
! f(right) > 0, but because f(x) is DECREASING ==> f(left) > 0
! BOTH functions are POSITIVE at the bounds.
! NO root in f(left, right)

            status = -1
            zf_crash_left = .FALSE.
            zf_crash_hi = .TRUE.
          END IF
        END IF

        IF (status==-1) THEN
! NO root in [left, right]

          zf_current_status = -1
          local%common_current_status = status
          local%common_crash_left = zf_crash_left
          local%common_crash_hi = zf_crash_hi

          IF (zf_crash_left) THEN
! Set root to the LEFT bound

            x = local%common_small
          ELSE
! Set root to the RIGHT bound

            x = local%common_big
          END IF

          RETURN
        END IF

! There is a root in [left, right]

        x = local%step_xsave
        local%step_step = MAX(local%step_absstp,local%step_relstp*ABS(x))
        local%step_from = 3
        status = 1

        RETURN

30      CONTINUE

        local%step_yy = fx
        IF (local%step_yy==zero) THEN
          status = 0
          zf_current_status = 0
          local%step_qok = .TRUE.
          RETURN
        END IF
        local%step_qup = (local%step_qincr .AND. (local%step_yy<zero)) &
          .OR. ( .NOT. local%step_qincr .AND. (local%step_yy>zero))

!      WE MUST STEP HIGHER

        IF ( .NOT. local%step_qup) GO TO 60
        local%step_xlb = local%step_xsave
        local%step_xub = MIN(local%step_xlb+local%step_step, &
          local%common_big)

40      CONTINUE
        x = local%step_xub
        local%step_from = 4
        status = 1
        RETURN

50      CONTINUE
        local%step_yy = fx
        local%step_qbdd = (local%step_qincr .AND. (local%step_yy>=zero)) &
          .OR. ( .NOT. local%step_qincr .AND. (local%step_yy<=zero))
        local%step_qlim = local%step_xub >= local%common_big
        local%step_qcond = local%step_qbdd .OR. local%step_qlim
        IF ( .NOT. (local%step_qcond)) THEN
          local%step_step = local%step_stpmul*local%step_step
          local%step_xlb = local%step_xub
          local%step_xub = MIN(local%step_xlb+local%step_step, &
            local%common_big)
          GO TO 40
        END IF

        IF (local%step_qlim .AND. .NOT. local%step_qbdd) THEN
          status = -1
          zf_current_status = -1
          local%common_current_status = -1

          zf_crash_left = .FALSE.
          zf_crash_hi = .NOT. local%step_qincr

          local%common_crash_left = zf_crash_left
          local%common_crash_hi = zf_crash_hi

          x = local%common_big

          RETURN
        END IF
        GO TO 90

60      CONTINUE
        local%step_xub = local%step_xsave
        local%step_xlb = MAX(local%step_xub-local%step_step, &
          local%common_small)

70      CONTINUE
        x = local%step_xlb
        local%step_from = 5
        status = 1
        RETURN

80      CONTINUE
        local%step_yy = fx
        local%step_qbdd = (local%step_qincr .AND. (local%step_yy<=zero)) &
          .OR. ( .NOT. local%step_qincr .AND. (local%step_yy>=zero))
        local%step_qlim = local%step_xlb <= local%common_small
        local%step_qcond = local%step_qbdd .OR. local%step_qlim

        IF ( .NOT. (local%step_qcond)) THEN
          local%step_step = local%step_stpmul*local%step_step
          local%step_xub = local%step_xlb
          local%step_xlb = MAX(local%step_xub-local%step_step, &
            local%common_small)
          GO TO 70
        END IF

        IF (local%step_qlim .AND. .NOT. local%step_qbdd) THEN
          status = -1
          zf_current_status = -1
          local%common_current_status = -1

          zf_crash_left = .TRUE.
          zf_crash_hi = local%step_qincr

          local%common_crash_left = zf_crash_left
          local%common_crash_hi = zf_crash_hi

          x = local%common_small

          RETURN
        END IF

!     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.

90      CONTINUE
        status = 0
        zf_current_status = 0
        local%common_current_status = 0
        local%common_called_from_step = .TRUE.

100     CALL rc_interval_local_zf(status,x,fx,local)

        local%step_from = 6

        IF (status/=1) status = 0
        RETURN

      CONTAINS

!.....................................................................

        FUNCTION qxmon(zx,zy,zz)
! .. Function Return Value ..
          LOGICAL :: qxmon
! ..
! .. Scalar Arguments ..
          REAL (dpkind), INTENT (IN) :: zx, zy, zz
! ..
          qxmon = zx <= zy .AND. zy <= zz
        END FUNCTION qxmon

!.....................................................................

      END SUBROUTINE rc_step_local_zf

!*********************************************************************

      SUBROUTINE rc_step_zf(status,x,fx,local)
! .. Structure Arguments ..
        TYPE (zf_locals), OPTIONAL :: local
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: fx
        REAL (dpkind), INTENT (OUT) :: x
        INTEGER, INTENT (INOUT) :: status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        IF (PRESENT(local)) THEN
          CALL rc_step_local_zf(status,x,fx,local)
        ELSE
          CALL rc_step_local_zf(status,x,fx,default_locals_zf)
        END IF


        RETURN

      END SUBROUTINE rc_step_zf

!*********************************************************************

      SUBROUTINE set_zero_finder(low_limit,hi_limit,abs_tol,rel_tol, &
          abs_step,rel_step,step_multiplier,local)
! .. Structure Arguments ..
        TYPE (zf_locals), OPTIONAL :: local
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: abs_step, abs_tol, &
          hi_limit, low_limit, rel_step, rel_tol, step_multiplier
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! .. Local Scalars ..
        LOGICAL :: locals_present
! ..
!!! Set flag to show that this routine has been called
!!! and that the interval routine was not called from the step routine
!!! If the interval routine was so callled, the flag will be reset
!!! by the step routine
        locals_present = PRESENT(local)

        IF (locals_present) THEN
          local%common_tolerances_set = common_tolerances_are_set
          local%common_called_from_step = .FALSE.
        ELSE
          default_locals_zf%common_tolerances_set = &
            common_tolerances_are_set
          default_locals_zf%common_called_from_step = .FALSE.
        END IF

!!!     Low Limit

        IF (PRESENT(low_limit)) THEN
          IF (locals_present) THEN
            local%common_small = low_limit
          ELSE
            default_locals_zf%common_small = low_limit
          END IF
        ELSE
          IF (locals_present) THEN
            local%common_small = default_low_limit
          ELSE
            default_locals_zf%common_small = default_low_limit
          END IF
        END IF

!!!     Hi Limit

        IF (PRESENT(hi_limit)) THEN
          IF (locals_present) THEN
            local%common_big = hi_limit
          ELSE
            default_locals_zf%common_big = hi_limit
          END IF
        ELSE
          IF (locals_present) THEN
            local%common_big = default_hi_limit
          ELSE
            default_locals_zf%common_big = default_hi_limit
          END IF
        END IF

!!!     Absolute tolerance

        IF (PRESENT(abs_tol)) THEN
          IF (locals_present) THEN
            local%common_abstol = abs_tol
          ELSE
            default_locals_zf%common_abstol = abs_tol
          END IF
        ELSE
          IF (locals_present) THEN
            local%common_abstol = default_abs_tol
          ELSE
            default_locals_zf%common_abstol = default_abs_tol
          END IF
        END IF

!     Relative tolerance

        IF (PRESENT(rel_tol)) THEN
          IF (locals_present) THEN
            local%common_reltol = rel_tol
          ELSE
            default_locals_zf%common_reltol = rel_tol
          END IF
        ELSE
          IF (locals_present) THEN
            local%common_reltol = default_rel_tol
          ELSE
            default_locals_zf%common_reltol = default_rel_tol
          END IF
        END IF

!     Absolute step size

        IF (PRESENT(abs_step)) THEN
          IF (locals_present) THEN
            local%step_absstp = abs_step
          ELSE
            default_locals_zf%step_absstp = abs_step
          END IF
        ELSE
          IF (locals_present) THEN
            local%step_absstp = default_abs_step
          ELSE
            default_locals_zf%step_absstp = default_abs_step
          END IF
        END IF

!     Relative step size

        IF (PRESENT(rel_step)) THEN
          IF (locals_present) THEN
            local%step_relstp = rel_step
          ELSE
            default_locals_zf%step_relstp = rel_step
          END IF
        ELSE
          IF (locals_present) THEN
            local%step_relstp = default_rel_step
          ELSE
            default_locals_zf%step_relstp = default_rel_step
          END IF
        END IF

!     Step size multiplier

        IF (PRESENT(step_multiplier)) THEN
          IF (locals_present) THEN
            local%step_stpmul = step_multiplier
          ELSE
            default_locals_zf%step_stpmul = step_multiplier
          END IF
        ELSE
          IF (locals_present) THEN
            local%step_stpmul = default_step_multiplier
          ELSE
            default_locals_zf%step_stpmul = default_step_multiplier
          END IF
        END IF

      END SUBROUTINE set_zero_finder

!*********************************************************************

      SUBROUTINE step_zf(f,y,answer,status,local)
! .. Structure Arguments ..
        TYPE (zf_locals), OPTIONAL :: local
! ..
! .. Non-Generic Interface Blocks ..
        INTERFACE
          FUNCTION f(x)
! .. Use Statements ..
            USE biomath_constants_mod
! ..
! .. Function Return Value ..
            REAL (dpkind) :: f
! ..
! .. Scalar Arguments ..
            REAL (dpkind), INTENT (IN) :: x
! ..
          END FUNCTION f
        END INTERFACE
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (INOUT) :: answer
        REAL (dpkind), INTENT (IN) :: y
        INTEGER, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: fx, x
        LOGICAL :: locals_present
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        locals_present = PRESENT(local)

        status = 0
        x = answer

        DO
          IF (locals_present) THEN
            CALL rc_step_local_zf(status,x,fx,local)
          ELSE
            CALL rc_step_local_zf(status,x,fx,default_locals_zf)
          END IF

          IF (status/=1) EXIT

          fx = f(x) - y

        END DO

        IF (status==0) answer = x

        RETURN

      END SUBROUTINE step_zf

!*********************************************************************

    END MODULE zero_finder
