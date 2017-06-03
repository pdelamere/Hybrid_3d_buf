    MODULE cdf_aux_mod
! .. Use Statements ..
      USE biomath_constants_mod
      USE biomath_interface_mod
      USE zero_finder
! ..
! .. Default Accessibility ..
      PUBLIC
! ..
! .. Private Statements ..
      PRIVATE :: bdf, bpnonc, bprob, immense, large, r_immense, sdf, &
        spnonc, sprob
! ..
! .. Parameters ..
      REAL (dpkind), PARAMETER :: bdf = 1.0E10_dpkind
      REAL (dpkind), PARAMETER :: bpnonc = 1.0E4_dpkind
      REAL (dpkind), PARAMETER :: immense = 1.0E100_dpkind
      REAL (dpkind), PARAMETER :: large = 1.0E10_dpkind
      REAL (dpkind), PARAMETER :: r_immense = 1.0E-100_dpkind
      REAL (dpkind), PARAMETER :: sdf = 1.0E-3_dpkind
      REAL (dpkind), PARAMETER :: spnonc = zero
      REAL (dpkind), PARAMETER :: sprob = 1.0E-10_dpkind
      REAL (dpkind), PARAMETER :: bprob = one - sprob
! ..
! .. Generic Interface Blocks ..
      INTERFACE in_range
        MODULE PROCEDURE dbl_in_range, int_in_range
      END INTERFACE
! ..
! .. Derived Type Declarations ..
! .. Local Structures
! Each parameter has this associated fields:
!   * name
!   * no_check: for this value of 'which'
!   * low_bound: minimum value for the parameter
!   * high_bound: high value for the parameter
      TYPE :: one_parameter
        CHARACTER (31) :: name
        INTEGER :: no_check
        REAL (dpkind) :: low_bound
        REAL (dpkind) :: high_bound
      END TYPE one_parameter
      TYPE :: the_distribution
        CHARACTER (31) :: name
        INTEGER :: max_which
        INTEGER :: nparam
        TYPE (one_parameter) :: parameters(6)
      END TYPE the_distribution
! ..
! .. Dependents ..
      TYPE (the_distribution), PARAMETER :: the_beta = the_distribution( &
        'cdf_beta',4,6,(/one_parameter('cum',1,zero,one), &
        one_parameter('ccum',1,zero,one),one_parameter('x',2,zero,one), &
        one_parameter('cx',2,zero,one),one_parameter('a',3,sprob,bdf), &
        one_parameter('b',4,sprob,bdf)/))
      TYPE (the_distribution), PARAMETER :: the_binomial = &
        the_distribution('cdf_binomial',4,6,(/one_parameter('cum',1,zero, &
        bprob),one_parameter('ccum',1,sprob,one), &
        one_parameter('s',2,zero,large),one_parameter('n',3,zero,large), &
        one_parameter('pr',4,zero,one),one_parameter('cpr',4,zero,one)/))
      TYPE (the_distribution), PARAMETER :: the_chi_square = &
        the_distribution('cdf_chi_square',3,4,(/one_parameter('cum',1, &
        zero,bprob),one_parameter('ccum',1,sprob,one), &
        one_parameter('x',2,zero,immense),one_parameter('df',3,sdf,bdf), &
        one_parameter(' ',0,zero,zero),one_parameter(' ',0,zero,zero)/))
      TYPE (the_distribution), PARAMETER :: the_dummy_binomial = &
        the_distribution('cdf_binomial',0,2,(/one_parameter('pr',0,zero, &
        half),one_parameter('cpr',0,zero,half),one_parameter(' ',0,zero, &
        zero),one_parameter(' ',0,zero,zero),one_parameter(' ',0,zero, &
        zero),one_parameter(' ',0,zero,zero)/))
      TYPE (the_distribution), PARAMETER :: the_f = the_distribution( &
        'cdf_f',2,5,(/one_parameter('cum',1,zero,bprob), &
        one_parameter('ccum',1,sprob,one),one_parameter('f',2,zero, &
        immense),one_parameter('dfn',3,sdf,bdf), &
        one_parameter('dfd',4,sdf,bdf),one_parameter(' ',0,zero,zero)/))
      TYPE (the_distribution), PARAMETER :: the_gamma = the_distribution( &
        'cdf_gamma',4,5,(/one_parameter('cum',1,zero,bprob), &
        one_parameter('ccum',1,sprob,one),one_parameter('x',2,zero, &
        immense),one_parameter('shape',3,sprob,immense), &
        one_parameter('scale',4,sprob,immense),one_parameter(' ',0,zero, &
        zero)/))
      TYPE (the_distribution), PARAMETER :: the_negative_binomial = &
        the_distribution('cdf_negative_binomial',4,6, &
        (/one_parameter('cum',1,zero,bprob),one_parameter('ccum',1,sprob, &
        one),one_parameter('f',2,zero,large),one_parameter('s',3,zero, &
        large),one_parameter('pr',4,zero,one),one_parameter('cpr',4,zero, &
        one)/))
      TYPE (the_distribution), PARAMETER :: the_non_central_chi_square = &
        the_distribution('cdf_nc_chisq',4,5,(/one_parameter('cum',1,zero, &
        bprob),one_parameter('ccum',1,sprob,one), &
        one_parameter('x',2,zero,immense),one_parameter('df',3,sdf,bdf), &
        one_parameter('pnonc',4,spnonc,bpnonc),one_parameter(' ',0,zero, &
        zero)/))
      TYPE (the_distribution), PARAMETER :: the_non_central_f = &
        the_distribution('cdf_nc_f',3,6,(/one_parameter('cum',1,zero, &
        bprob),one_parameter('ccum',1,sprob,one), &
        one_parameter('f',2,zero,immense),one_parameter('dfn',3,sdf,bdf), &
        one_parameter('dfd',4,sdf,bdf),one_parameter('pnonc',5,spnonc, &
        bpnonc)/))
      TYPE (the_distribution), PARAMETER :: the_non_central_t = &
        the_distribution('cdf_non_central_t',4,5, &
        (/one_parameter('cum',1,sprob,bprob),one_parameter('ccum',1,sprob &
        ,bprob),one_parameter('t',2,-immense,immense), &
        one_parameter('df',3,sdf,bdf),one_parameter('pnonc',4,spnonc, &
        bpnonc),one_parameter(' ',0,zero,zero)/))
      TYPE (the_distribution), PARAMETER :: the_normal = the_distribution &
        ('cdf_normal',4,5,(/one_parameter('cum',1,sprob,bprob), &
        one_parameter('ccum',1,sprob,bprob),one_parameter('x',2,-immense, &
        immense),one_parameter('mean',3,-immense,immense), &
        one_parameter('sd',4,sprob,immense),one_parameter(' ',0,zero,zero &
        )/))
      TYPE (the_distribution), PARAMETER :: the_poisson = &
        the_distribution('cdf_poisson',3,4,(/one_parameter('cum',1,zero, &
        bprob),one_parameter('ccum',1,sprob,one), &
        one_parameter('s',2,zero,immense),one_parameter('lambda',3,sprob, &
        immense),one_parameter(' ',0,zero,zero), &
        one_parameter(' ',0,zero,zero)/))
      TYPE (the_distribution), PARAMETER :: the_t = the_distribution( &
        'cdf_t',3,4,(/one_parameter('cum',1,sprob,bprob), &
        one_parameter('ccum',1,sprob,bprob),one_parameter('t',2,-immense, &
        immense),one_parameter('df',3,sdf,bdf),one_parameter(' ',0,zero, &
        zero),one_parameter(' ',0,zero,zero)/))
! ..
    CONTAINS

!*********************************************************************

      FUNCTION add_to_one(x,y,routine_name,x_name,y_name,bad_status, &
          status)
!------------------------------------------------------------------------
! The sum of complementary parameters x and y MUST be 1.0
! (or very close to one.)

! If sum is 1.0, then function returns .TRUE.

! If it is NOT, then:
!   * if STATUS is present (the caller routine had status in its
!     argument list), then status is set to bad_status and the function
!     returns .FALSE.  The user MUST check the status after the return
!     from this function to make sure that x+y is indeed 1.0
!   * if STATUS is NOT present, in order to prevent following nonsensical
!     computations, the application BREAKS here with an error message.
!------------------------------------------------------------------------
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: x, y
        INTEGER, INTENT (IN) :: bad_status
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        CHARACTER (*), INTENT (IN) :: routine_name, x_name, y_name
! ..
! .. Function Return Value ..
        LOGICAL :: add_to_one
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, EPSILON, PRESENT
! ..
        IF (ABS(((x+y)-half)-half)>three*EPSILON(one)) THEN
          add_to_one = .FALSE.
          IF (PRESENT(status)) THEN
            status = bad_status
          ELSE
            WRITE (*,90000) routine_name, x_name, y_name, x_name, x, &
              y_name, y
90000       FORMAT (' Error in routine: ',A/A,' + ',A, &
              ' must add to one'/'Value of ',A,' is: ',G15.8/'Value of ', &
              A,' is: ',G15.8/'Program aborting'/)
            STOP 'Error in cdf routine -- see output'
          END IF
        ELSE
          add_to_one = .TRUE.
        END IF

        RETURN

      END FUNCTION add_to_one

!*********************************************************************

      SUBROUTINE cdf_finalize_status(local,status)
!------------------------------------------------------------------------
! This routine should be called at the end of EACH root finding in the
! cdf_ routines.  It mainly sets the status of the zero finding.
! Again, the user is URGED to check the return value of this routine.
!------------------------------------------------------------------------
! .. Scalar Arguments ..
        INTEGER, INTENT (OUT) :: status
! ..
! .. Structure Arguments ..
        TYPE (zf_locals), INTENT (IN) :: local
! ..
! .. Local Scalars ..
        REAL (dpkind) :: left_end, right_end
        INTEGER :: zf_status
        LOGICAL :: crash_hi, crash_left
! ..
        CALL final_zf_state(local,zf_status,crash_left,crash_hi,left_end, &
          right_end)

        IF (zf_status==0) THEN
          status = 0
        ELSE
          IF (crash_left) THEN
! The root was set in zero_finder to: x=left_bound
            status = -50
          ELSE
! The root was set in zero_finder to: x=right_bound
            status = 50
          END IF
        END IF

        RETURN

      END SUBROUTINE cdf_finalize_status

!*********************************************************************

      SUBROUTINE cdf_set_zero_finder(distrib,which,local)
! Called before solving a problem.
! Sets the bounds for the parameter to be found.
! The bounds are in the structure distrib, corresponding to the
! "which-th" parameter
! .. Parameters ..
        REAL (dpkind), PARAMETER :: absstp = half
        REAL (dpkind), PARAMETER :: atol = 1.0E-50_dpkind
        REAL (dpkind), PARAMETER :: relstp = half
        REAL (dpkind), PARAMETER :: reltol = 1.0E-8_dpkind
        REAL (dpkind), PARAMETER :: stpmul = five
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: which
! ..
! .. Structure Arguments ..
        TYPE (the_distribution), INTENT (IN) :: distrib
        TYPE (zf_locals), INTENT (OUT) :: local
! ..
        CALL set_zero_finder(low_limit=distrib%parameters(which) &
          %low_bound,hi_limit=distrib%parameters(which)%high_bound, &
          abs_step=absstp,rel_step=relstp,step_multiplier=stpmul, &
          abs_tol=atol,rel_tol=reltol,local=local)

      END SUBROUTINE cdf_set_zero_finder

!*********************************************************************

      SUBROUTINE check_complements(x,y,routine_name,x_name,y_name, &
          local_x,local_y,set_values,bad_status,status)
! At least one of the optional parameters x and y must be defined in
! the caller.
! If NONE is defined, the program STOPS.
! If ONE parameter is defined:
!   * if:   set_values is FALSE, return to the caller
!     else: set the non-defined parameter to 1 minus the value of the
!          defined parameter.
! If BOTH parameters are defined:
!   * if   : set_values is FALSE, return to the caller
!     else : set local_x=x and local_y=y
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: local_x, local_y
        REAL (dpkind), OPTIONAL, INTENT (IN) :: x, y
        INTEGER, INTENT (IN) :: bad_status
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL :: set_values
        CHARACTER (*), INTENT (IN) :: routine_name, x_name, y_name
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! .. Local Scalars ..
        INTEGER :: count
! ..
        count = 0

        IF (PRESENT(x)) count = count + 1
        IF (PRESENT(y)) count = count + 1

        SELECT CASE (count)

        CASE (0)

! None of the complementary parameters are present in the call

          IF (PRESENT(status)) THEN
            status = bad_status
          ELSE
            WRITE (*,90000) routine_name, x_name, y_name
            STOP 'Complement error -- see output'
          END IF

90000     FORMAT (' In routine: ',A/' Neither argument: ',A,' or: ',A, &
            ' is present'/' At least one is required'/ &
            ' Program aborting')

        CASE (1)

          IF ( .NOT. set_values) RETURN

          IF (PRESENT(x)) THEN
            local_x = x
            local_y = one - x
          ELSE
            local_y = y
            local_x = one - y
          END IF

        CASE (2)

          IF ( .NOT. set_values) RETURN

          local_x = x
          local_y = y

        END SELECT

        RETURN

      END SUBROUTINE check_complements

!*********************************************************************

      FUNCTION dbl_in_range(value,lo,hi,routine_name,arg_name,bad_status, &
          status)
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: hi, lo, value
        INTEGER, INTENT (IN) :: bad_status
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        CHARACTER (*) :: arg_name, routine_name
! ..
! .. Local Scalars ..
        LOGICAL :: has_status
! ..
! .. Function Return Value ..
        LOGICAL :: dbl_in_range
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        dbl_in_range = .TRUE.
        has_status = PRESENT(status)

        IF (value>hi) THEN
          dbl_in_range = .FALSE.
          IF (has_status) THEN
            status = bad_status
            RETURN
          ELSE
            WRITE (*,90000) routine_name, arg_name, arg_name, value, hi
90000       FORMAT (' Input argument out of range in routine: ', &
              A/' Argument: ',A,' above upper bound'/' The value of ',A, &
              ' is: ',G15.8/' The upper bound on the value is: ', &
              G15.8/' The program is aborting')
            STOP 'Bound error on cdf routine -- see output'
          END IF
        END IF

        IF (value<lo) THEN
          dbl_in_range = .FALSE.
          IF (has_status) THEN
            status = bad_status
            RETURN
          ELSE
            WRITE (*,90100) routine_name, arg_name, arg_name, value, lo
90100       FORMAT (' Input argument out of range in routine: ', &
              A/' Argument: ',A,' below lower bound'/' The value of ',A, &
              ' is: ',G15.8/' The lower bound on the value is: ', &
              G15.8/' The program is aborting')
            STOP 'Bound error on cdf routine -- see output'
          END IF
        END IF

      END FUNCTION dbl_in_range

!*********************************************************************

      FUNCTION int_in_range(value,lo,hi,routine_name,arg_name,bad_status, &
          status)
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: bad_status, hi, lo, value
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        CHARACTER (*), INTENT (IN) :: arg_name, routine_name
! ..
! .. Local Scalars ..
        LOGICAL :: has_status
! ..
! .. Function Return Value ..
        LOGICAL :: int_in_range
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        int_in_range = .TRUE.
        has_status = PRESENT(status)

        IF (value>hi) THEN
          IF (has_status) THEN
            int_in_range = .FALSE.
            status = bad_status
            RETURN
          ELSE
            WRITE (*,*) ' Input argument error in: ', routine_name
            WRITE (*,*) ' Argument ', arg_name, ' above upper bound'
            WRITE (*,*) arg_name, ' value: ', value
            WRITE (*,*) ' Upper bound: ', hi
            STOP &
              'Stop: uncaught bound error on cdf routine -- see output'
          END IF
        END IF

        IF (value<lo) THEN
          IF (has_status) THEN
            int_in_range = .FALSE.
            status = bad_status
            RETURN
          ELSE
            WRITE (*,*) ' Input argument error in: ', routine_name
            WRITE (*,*) ' Argument ', arg_name, ' below lower bound'
            WRITE (*,*) arg_name, ' value: ', value
            WRITE (*,*) ' Lower bound: ', lo
            STOP &
              'Stop: uncaught bound error on cdf routine -- see output'
          END IF
        END IF

        RETURN

      END FUNCTION int_in_range

!*********************************************************************

      SUBROUTINE validate_parameters(distrib,which,params,status)
!        Generic subroutine which performs validation for
!        ALL the parameters of ALL distributions.
! .. Local Scalars ..
        INTEGER :: i, nparam
        CHARACTER (31) :: tn1, tn2
! ..
! .. Structure Arguments ..
        TYPE (the_distribution) :: distrib
! ..
! .. Scalar Arguments ..
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        INTEGER :: which
! ..
! .. Array Arguments ..
        REAL (dpkind), INTENT (IN) :: params(6)
! ..
        nparam = distrib%nparam

! Check for which to be within limits

        tn1 = 'which'

        IF ( .NOT. in_range(value=which,lo=1,hi=distrib%max_which, &
          routine_name=distrib%name,arg_name=tn1,bad_status=-1, &
          status=status)) RETURN

! If the unknown is not cdf (which=1), check that cum and ccum add to one

        IF (which/=1) THEN
          tn1 = 'cum'
          tn2 = 'ccum'
          IF ( .NOT. add_to_one(params(1),params( &
            2),distrib%name,tn1,tn2,3,status)) RETURN
        END IF

        DO i = 1, nparam

          IF (which/=distrib%parameters(i)%no_check) THEN
            IF ( .NOT. in_range(value=params(i),lo=distrib%parameters( &
              i)%low_bound,hi=distrib%parameters( &
              i)%high_bound,routine_name=distrib%name, &
              arg_name=distrib%parameters(i)%name,bad_status=-(i+ &
              1),status=status)) RETURN
          END IF

        END DO

      END SUBROUTINE validate_parameters

!*********************************************************************

      SUBROUTINE which_miss(which,routine_name,arg_name)
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: which
        CHARACTER (*), INTENT (IN) :: arg_name, routine_name
! ..
        WRITE (*,90000) routine_name, which, arg_name, arg_name

90000   FORMAT (' In routine: ',A/' which has value: ',I0, &
          ' indicating that ',A,'is to be calculated'/'But argument: ',A, &
          ' is missing'/'The program is aborting')

        STOP 'Which missing_arg error - see output'

        RETURN

      END SUBROUTINE which_miss

!*********************************************************************

    END MODULE cdf_aux_mod
