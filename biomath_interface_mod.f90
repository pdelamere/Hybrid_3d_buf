    MODULE biomath_interface_mod
! This module contains procedures (methods) for
! interaction with the user.
! .. Use Statements ..
      USE biomath_constants_mod
      USE biomath_strings_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: always_print, clear_screen, get_character, &
        get_list_double, get_string, get_yn, hold, message_format, &
        num_subs, print_message_format, print_off, prompt, report_unit, &
        write_array, write_error, write_message
      PUBLIC :: get_numbers
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: max_try = 3, mx_format_length = 5000, &
        mx_subs = 100, mx_sub_length = 100
! ..
! .. Local Scalars ..
! Refer to :Local Scalars
! Refer to :Local Scalars
      INTEGER :: num_subs, report_unit
      INTEGER, SAVE :: print_level = 1
      LOGICAL :: always_print, print_off
      LOGICAL, SAVE :: clear_screen_before_print = .FALSE., &
        format_printed
      CHARACTER (mx_format_length) :: message_format
! ..
! .. Local Arrays ..
      INTEGER :: sub_pos(2*mx_subs)
      CHARACTER (mx_sub_length) :: sub_strings(mx_subs)
! ..
! .. Generic Interface Blocks ..
      INTERFACE get_numbers
        MODULE PROCEDURE get_double_array, get_integer_array, &
          get_real_array, get_double_scalar, get_integer_scalar, &
          get_real_scalar
      END INTERFACE
! ..
    CONTAINS

!*********************************************************************

      SUBROUTINE clear_screen
! writes out WINDOW_SIZE (in this module) blank lines
! to STDOUT
! .. Local Scalars ..
        INTEGER :: i
        INTEGER, SAVE :: window_size = 24
! ..
        DO i = 1, window_size

          WRITE (stdout,*)

        END DO

      END SUBROUTINE clear_screen

!*********************************************************************

      SUBROUTINE hold
! displays hold message to module unit STDOUT,
! waits for user input from module unit STDIN,
! and returns - all this should take place on one line.
! .. Local Scalars ..
        INTEGER :: ios
        CHARACTER (1) :: dummy_string
! ..
        WRITE (stdout, &
          '(/1X,''Press the Return or Enter key to continue ...'')')

        READ (stdin,'(A)',iostat=ios) dummy_string

        RETURN

      END SUBROUTINE hold

!*********************************************************************

      SUBROUTINE prompt
! Writes a prompt to module unit STDOUT, without advancing the line.
        WRITE (stdout,'(1X,''> '')',advance='NO')

        RETURN

      END SUBROUTINE prompt

!*********************************************************************

      SUBROUTINE get_character(mssg,chars,num_chars,which,qerr)
! Displays mssg to user (explaining character choices), gets response
! from user, evaluates first non-blank character, returns index of
! response character in chars(1:num_chars) as which.
! chars should all be lower case.  Response letter is translated to
! lower case before evaluation.
! If no good answer in max_bad tries,
! qerr = .TRUE., otherwise qerr = .FALSE.
! Note: if num_chars is not positive and at most the length of chars,
! undefined behavior may result.
! .. Parameters ..
        INTEGER, PARAMETER :: line_len = 80, max_bad = 3
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: num_chars
        INTEGER, INTENT (OUT) :: which
        LOGICAL, INTENT (OUT) :: qerr
        CHARACTER (*), INTENT (IN) :: chars, mssg
! ..
! .. Local Scalars ..
        INTEGER :: i, ios, j
        CHARACTER (1) :: answer
        CHARACTER (line_len) :: line
! ..
! .. Intrinsic Functions ..
        INTRINSIC INDEX
! ..
        qerr = .FALSE.

GET_ANSWER: DO i = 1, max_bad

          IF (i>1) WRITE (stdout,'(/5X,''Please try again ...'')')

          WRITE (stdout,mssg)

          WRITE (stdout,'(1X,''Please enter one of ['',A,'']:'')', &
            advance='NO') chars(1:num_chars)

          CALL prompt

          READ (stdin,'(A)',iostat=ios) line

          IF (ios/=0) THEN

            WRITE (stdout,'(/5X,''Could not read the response.'')')

            CYCLE GET_ANSWER

          END IF

          DO j = 1, line_len

            IF (line(j:j)/=' ') EXIT

          END DO

          IF (j>line_len) THEN

            WRITE (stdout,'(/5X''A blank line is not allowed.'')')

            CYCLE GET_ANSWER

          END IF

          answer = lower_case_char(line(j:j))

          which = INDEX(chars(1:num_chars),answer)

          IF (which>0) RETURN

          WRITE (stdout,'(/5X,A,'' is not allowed.'')') answer

        END DO GET_ANSWER

        qerr = .TRUE.

        RETURN

      END SUBROUTINE get_character

!*********************************************************************

      SUBROUTINE get_yn(mssg,qyes,qerr)
! Uses mssg to elicit a Y/N response from the user with get_character.
! If get_character fails, qerr == .TRUE., otherwise qerr == .FALSE.
! The value for qgetyn is undefined if qerr = .TRUE.
! otherwise is .TRUE. if user enters Y
!              .FALSE. if user enters N
! .. Local Scalars ..
        INTEGER :: which
! ..
! .. Scalar Arguments ..
        LOGICAL, INTENT (OUT) :: qerr, qyes
        CHARACTER (*), INTENT (IN) :: mssg
! ..
! .. Executable code ..
        CALL get_character(mssg,'yn',2,which,qerr)

        IF ( .NOT. qerr) qyes = (which==1)

        RETURN

      END SUBROUTINE get_yn

!*********************************************************************

      SUBROUTINE get_string(fmt,nline,string,qnulok,qfail)
!-----------------------------------------------------------------------
!     Name: Get String from User
!     Description:
!          Gets a character string from the user. Ignores lines
!     starting with a "#". Ignores any text entered after a "#".
!          If QNULOK is FALSE then blank lines are not allowed.
!     The user is warned if he enteres a blank line and asked
!     to try again. The user gets 3 tries to enter an acceptable
!     line.
!     Arguments:
!     FMT -> Array of character strings defining the format for the
!            query statements. Query will be made as write(*,fmt(i))
!            for every line in FMT.
!            CHARACTER*(*) FMT(*)
!     NLINE -> Number of lines in the query message.
!              If NLINE is 0 then no query message will be written.
!              INTEGER NLINE
!     STRING <- Characters entered by user.
!               CHARACTER*(*) STRING
!     QNULOK -> If TRUE then a null line (line with all blanks or
!               comments) is OK. If FALSE then a null line
!               will give the user an error message.
!               LOGICAL QNULOK
!     QFAIL <- TRUE if read failed 3 times. Note that get_string
!              does not inform the user of this state.
!              LOGICAL QFAIL
!-----------------------------------------------------------------------
! .. Parameters ..
        INTEGER, PARAMETER :: itry = 3
! ..
! .. Scalar Arguments ..
        INTEGER :: nline
        LOGICAL :: qfail, qnulok
        CHARACTER (*) :: string
! ..
! .. Array Arguments ..
        CHARACTER (*) :: fmt(*)
! ..
! .. Local Scalars ..
        INTEGER :: i, nbad
        LOGICAL :: qblank, qcomnt
! ..
! .. Intrinsic Functions ..
        INTRINSIC LEN
! ..
        nbad = 0

!     Query user.

10      CONTINUE
        IF (nline>=1) THEN
          DO i = 1, nline
            WRITE (*,fmt(i))
          END DO
        END IF
        CALL prompt

!     Read in the string.

20      CONTINUE
        READ (*,'(A)',err=30) string

!     Check to see if the string is a comment.
!     If so read in another line.

        IF (string(1:1)=='#') GO TO 20

!     Loop through the string.
!     Blank out any characters past a "#".
!     Determine if the string is blank.

        qcomnt = .FALSE.
        qblank = .TRUE.
        DO i = 1, LEN(string)
          IF (qcomnt) THEN
            string(i:i) = ' '
          ELSE
            IF (string(i:i)=='#') THEN
              string(i:i) = ' '
              qcomnt = .TRUE.
            ELSE IF (qblank) THEN
              IF (string(i:i)/=' ') qblank = .FALSE.
            END IF
          END IF
        END DO

!     If the string is blank and that is not allowed.
!     Then send the user an error message. Otherwise,
!     return with the string.

        IF ( .NOT. qnulok .AND. qblank) THEN
          WRITE (*,90000) 'Blank line not allowed.'
        ELSE
          qfail = .FALSE.
          RETURN
        END IF

!     If we get here there was an error

30      CONTINUE
        nbad = nbad + 1

!     If there were too many bad entries return with QFAIL TRUE.

        IF (nbad>=itry) THEN
          qfail = .TRUE.
          RETURN
        END IF

!     Otherwise print an error message and try again.

        WRITE (*,90000) 'Bad Entry. Try again.'
        GO TO 10
!     ..
!     .. Format Statements ..
90000   FORMAT (1X,1A)

      END SUBROUTINE get_string

!!!**********************************************************************

      SUBROUTINE get_double_array(format,lo,lo_eq_ok,hi,hi_eq_ok,x, &
          number_wanted,failed)
!!!
!!!SUMMARY
!!!
!!!get_numbers is a  generic Fortran 90 routine  that prompts a user then
!!!reads   a  specified number  of  values   into  array  x.  Optionally,
!!!get_numbers checks  high and low bounds  on the numbers.   The type of
!!!each element of x and of  lo and hi must be  the same and either
!!!double precision, real, or integer.  x, lo, and hi can all be scalars
!!!or all can be arrays of rank 1.
!!!
!!!DETAILS
!!!
!!!    SUBROUTINE get_numbers(format,lo,lo_eq_ok,hi,hi_eq_ok,x,
!!!               number_wanted, failed)
!!!
!!!    CHARACTER (*), OPTIONAL, INTENT (IN) :: format
!!!
!!!If present,   a  character string  containing   a Fortran format.  The
!!!format will be printed to standard output (*)  to prompt the user.  If
!!!the argument is missing, no prompt message will be printed.
!!!
!!!      <TYPE>, OPTIONAL, INTENT (IN) :: lo(:)     or
!!!      <TYPE>, OPTIONAL, INTENT (IN) :: lo
!!!
!!!<TYPE> should be the same as the type of x: double precision, integer,
!!!or real.
!!!
!!!Low bounds for the values entered.  If lo has size 1, the single value
!!!is used to check all input values.   If lo does not  have size 1, then
!!!it should have size  of at least the  number of values to be obtained.
!!!In this case the i'th value obtained  from the user is checked against
!!!the i'th value  of lo.   If this argument  is missing,  no lower bound
!!!checking is performed.
!!!
!!!      LOGICAL, OPTIONAL, INTENT (IN) :: lo_eq_ok
!!!
!!!If .TRUE., a value can equal its low bound without triggering an error
!!!message; if .FALSE., it can't.  Default if missing is .TRUE.
!!!
!!!      <TYPE>, OPTIONAL, INTENT (IN) :: hi(:)     or
!!!      <TYPE>, OPTIONAL, INTENT (IN) :: hi
!!!
!!!High bounds for the values entered.  See lo.
!!!
!!!      LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok
!!!
!!!If   .TRUE., a value  can equal  its high bound  without triggering an
!!!error message; if .FALSE., it can't.  Default if missing is .TRUE.
!!!
!!!      <TYPE>, INTENT (OUT) :: x(:)     or
!!!      <TYPE>, INTENT (OUT) :: x
!!!
!!!<TYPE> should be one of: double precision, integer, or real. The array
!!!into which the values obtained from the user should be placed.
!!!
!!!      INTEGER, OPTIONAL, INTENT (IN) :: number_wanted
!!!
!!!The number of  values to be obtained from  the user.  If missing, this
!!!value is SIZE( x ).
!!!
!!!      LOGICAL, OPTIONAL, INTENT (OUT) :: failed
!!!
!!!If .TRUE. on exit, three attempts failed to obtain the values from the
!!!user.  If argument is missing,  this subroutine will abort (STOP) with
!!!an error message.
!!!
!!!**********************************************************************
! .. Scalar Arguments ..
        INTEGER, OPTIONAL, INTENT (IN) :: number_wanted
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Array Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: hi(:), lo(:)
        REAL (dpkind), INTENT (OUT) :: x(:)
! ..
! .. Local Scalars ..
        INTEGER :: check_hi_status, check_lo_status, i, ios, i_try, &
          use_n_wanted
        LOGICAL :: use_hi_eq_ok, use_lo_eq_ok
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT, SIZE
! ..
!     Query user.
        use_lo_eq_ok = .TRUE.
        IF (PRESENT(lo_eq_ok)) use_lo_eq_ok = lo_eq_ok

        use_hi_eq_ok = .TRUE.
        IF (PRESENT(hi_eq_ok)) use_hi_eq_ok = hi_eq_ok

!!!   CHECK_[LO|HI]_STATUS
!!!   -1 No check
!!!    1 Check against first and only element in bounds array
!!!   >1 Check against corresponding element in bounds array

        check_lo_status = -1
        IF (PRESENT(lo)) THEN
          check_lo_status = 1
          IF (SIZE(lo)/=1) check_lo_status = 2
        END IF

        check_hi_status = -1
        IF (PRESENT(hi)) THEN
          check_hi_status = 1
          IF (SIZE(hi)/=1) check_hi_status = 2
        END IF

        IF (PRESENT(number_wanted)) THEN
          use_n_wanted = number_wanted
        ELSE
          use_n_wanted = SIZE(x)
        END IF

TRY_LOOP: DO i_try = 1, max_try

          IF (i_try>1) WRITE (*,*) 'Please try again from the beginning.'

          IF (PRESENT(format)) WRITE (*,format)

          CALL prompt

!     Read in value.

          READ (*,*,iostat=ios) (x(i),i=1,use_n_wanted)
          IF (ios/=0) THEN
            WRITE (*,*) &
              'Some entered value not accepted as a legal number.'
            CYCLE TRY_LOOP
          END IF

!     Check bound on each entry in X.

CHECK_VALUES: DO i = 1, use_n_wanted

!        Check lower bound if required.

90000       FORMAT ('The ',I3, &
              '''th value entered < the lower bound allowed.'/ &
              'The value entered was: ',G20.7/'The lower bound is: ', &
              G20.7)
90100       FORMAT ('The ',I3, &
              '''th value entered <= the lower bound allowed.'/ &
              'The value entered was: ',G20.7/'The lower bound is: ', &
              G20.7)

LO_CHECK:   IF (check_lo_status>0) THEN
              IF (check_lo_status==1) THEN
                IF (use_lo_eq_ok .AND. x(i)<lo(1)) THEN
                  WRITE (*,90000) i, x(i), lo(1)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_lo_eq_ok .AND. x(i)<=lo(1)) THEN
                  WRITE (*,90100) i, x(i), lo(1)
                  CYCLE TRY_LOOP
                END IF
              ELSE
                IF (use_lo_eq_ok .AND. x(i)<lo(i)) THEN
                  WRITE (*,90000) i, x(i), lo(i)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_lo_eq_ok .AND. x(i)<=lo(i)) THEN
                  WRITE (*,90100) i, x(i), lo(i)
                  CYCLE TRY_LOOP
                END IF
              END IF

            END IF LO_CHECK

!        Check upper bound if required.

90200       FORMAT ('The ',I3, &
              '''th value entered > the upper bound allowed.'/ &
              'The value entered was: ',G20.7/'The upper bound is: ', &
              G20.7)
90300       FORMAT ('The ',I3, &
              '''th value entered >= the lower bound allowed.'/ &
              'The value entered was: ',G20.7/'The upper bound is: ', &
              G20.7)

HI_CHECK:   IF (check_hi_status>0) THEN
              IF (check_hi_status==1) THEN
                IF (use_hi_eq_ok .AND. x(i)>hi(1)) THEN
                  WRITE (*,90200) i, x(i), hi(1)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_hi_eq_ok .AND. x(i)>=hi(1)) THEN
                  WRITE (*,90300) i, x(i), hi(1)
                  CYCLE TRY_LOOP
                END IF
              ELSE
                IF (use_hi_eq_ok .AND. x(i)>hi(i)) THEN
                  WRITE (*,90200) i, x(i), hi(i)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_hi_eq_ok .AND. x(i)>=hi(i)) THEN
                  WRITE (*,90300) i, x(i), hi(i)
                  CYCLE TRY_LOOP
                END IF
              END IF

            END IF HI_CHECK

          END DO CHECK_VALUES

          IF (PRESENT(failed)) failed = .FALSE.
          RETURN

        END DO TRY_LOOP

!!! Failure if we reach here

        IF (PRESENT(failed)) THEN
          failed = .TRUE.
        ELSE
          WRITE (*,*) 'Too many failed attempts at entry in &
            &get_double_array. ABORT'
          STOP 'Too many failed attempts at entry in get_double_array. &
            &ABORT'
        END IF

      END SUBROUTINE get_double_array

!!!======================================================================

      SUBROUTINE get_integer_array(format,lo,lo_eq_ok,hi,hi_eq_ok,x, &
          number_wanted,failed)
! .. Scalar Arguments ..
        INTEGER, OPTIONAL, INTENT (IN) :: number_wanted
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Array Arguments ..
        INTEGER, OPTIONAL, INTENT (IN) :: hi(:), lo(:)
        INTEGER, INTENT (OUT) :: x(:)
! ..
! .. Local Scalars ..
        INTEGER :: check_hi_status, check_lo_status, i, ios, i_try, &
          use_n_wanted
        LOGICAL :: use_hi_eq_ok, use_lo_eq_ok
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT, SIZE
! ..
!     Query user.
        use_lo_eq_ok = .TRUE.
        IF (PRESENT(lo_eq_ok)) use_lo_eq_ok = lo_eq_ok

        use_hi_eq_ok = .TRUE.
        IF (PRESENT(hi_eq_ok)) use_hi_eq_ok = hi_eq_ok

!!!   CHECK_[LO|HI]_STATUS
!!!   -1 No check
!!!    1 Check against first and only element in bounds array
!!!   >1 Check against corresponding element in bounds array

        check_lo_status = -1
        IF (PRESENT(lo)) THEN
          check_lo_status = 1
          IF (SIZE(lo)/=1) check_lo_status = 2
        END IF

        check_hi_status = -1
        IF (PRESENT(hi)) THEN
          check_hi_status = 1
          IF (SIZE(hi)/=1) check_hi_status = 2
        END IF

        IF (PRESENT(number_wanted)) THEN
          use_n_wanted = number_wanted
        ELSE
          use_n_wanted = SIZE(x)
        END IF

TRY_LOOP: DO i_try = 1, max_try

          IF (i_try>1) WRITE (*,*) 'Please try again from the beginning.'

          IF (PRESENT(format)) WRITE (*,format)

          CALL prompt

!     Read in value.

          READ (*,*,iostat=ios) (x(i),i=1,use_n_wanted)
          IF (ios/=0) THEN
            WRITE (*,*) &
              'Some entered value not accepted as a legal number.'
            CYCLE TRY_LOOP
          END IF

!     Check bound on each entry in X.

CHECK_VALUES: DO i = 1, use_n_wanted

!        Check lower bound if required.

90000       FORMAT ('The ',I3, &
              '''th value entered < the lower bound allowed.'/ &
              'The value entered was: ',I20/'The lower bound is: ',I20)
90100       FORMAT ('The ',I3, &
              '''th value entered <= the lower bound allowed.'/ &
              'The value entered was: ',I20/'The lower bound is: ',I20)

LO_CHECK:   IF (check_lo_status>0) THEN
              IF (check_lo_status==1) THEN
                IF (use_lo_eq_ok .AND. x(i)<lo(1)) THEN
                  WRITE (*,90000) i, x(i), lo(1)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_lo_eq_ok .AND. x(i)<=lo(1)) THEN
                  WRITE (*,90100) i, x(i), lo(1)
                  CYCLE TRY_LOOP
                END IF
              ELSE
                IF (use_lo_eq_ok .AND. x(i)<lo(i)) THEN
                  WRITE (*,90000) i, x(i), lo(i)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_lo_eq_ok .AND. x(i)<=lo(i)) THEN
                  WRITE (*,90100) i, x(i), lo(i)
                  CYCLE TRY_LOOP
                END IF
              END IF

            END IF LO_CHECK

!        Check upper bound if required.

90200       FORMAT ('The ',I3, &
              '''th value entered > the upper bound allowed.'/ &
              'The value entered was: ',I20/'The upper bound is: ',I20)
90300       FORMAT ('The ',I3, &
              '''th value entered >= the upper bound allowed.'/ &
              'The value entered was: ',I20/'The upper bound is: ',I20)

HI_CHECK:   IF (check_hi_status>0) THEN
              IF (check_hi_status==1) THEN
                IF (use_hi_eq_ok .AND. x(i)>hi(1)) THEN
                  WRITE (*,90200) i, x(i), hi(1)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_hi_eq_ok .AND. x(i)>=hi(1)) THEN
                  WRITE (*,90300) i, x(i), hi(1)
                  CYCLE TRY_LOOP
                END IF
              ELSE
                IF (use_hi_eq_ok .AND. x(i)>hi(i)) THEN
                  WRITE (*,90200) i, x(i), hi(i)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_hi_eq_ok .AND. x(i)>=hi(i)) THEN
                  WRITE (*,90300) i, x(i), hi(i)
                  CYCLE TRY_LOOP
                END IF
              END IF

            END IF HI_CHECK

          END DO CHECK_VALUES

          IF (PRESENT(failed)) failed = .FALSE.
          RETURN

        END DO TRY_LOOP

!!! Failure if we reach here

        IF (PRESENT(failed)) THEN
          failed = .TRUE.
          RETURN
        ELSE
          WRITE (*,*) 'Too many failed attempts at entry in &
            &get_integer_array. ABORT'
          STOP 'Too many failed attempts at entry in get_integer_array. &
            &ABORT'
        END IF

      END SUBROUTINE get_integer_array

!!!======================================================================

      SUBROUTINE get_real_array(format,lo,lo_eq_ok,hi,hi_eq_ok,x, &
          number_wanted,failed)
! .. Scalar Arguments ..
        INTEGER, OPTIONAL, INTENT (IN) :: number_wanted
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Array Arguments ..
        REAL (spkind), OPTIONAL, INTENT (IN) :: hi(:), lo(:)
        REAL (spkind), INTENT (OUT) :: x(:)
! ..
! .. Local Scalars ..
        INTEGER :: check_hi_status, check_lo_status, i, ios, i_try, &
          use_n_wanted
        LOGICAL :: use_hi_eq_ok, use_lo_eq_ok
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT, SIZE
! ..
!     Query user.
        use_lo_eq_ok = .TRUE.
        IF (PRESENT(lo_eq_ok)) use_lo_eq_ok = lo_eq_ok

        use_hi_eq_ok = .TRUE.
        IF (PRESENT(hi_eq_ok)) use_hi_eq_ok = hi_eq_ok

!!!   CHECK_[LO|HI]_STATUS
!!!   -1 No check
!!!    1 Check against first and only element in bounds array
!!!   >1 Check against corresponding element in bounds array

        check_lo_status = -1
        IF (PRESENT(lo)) THEN
          check_lo_status = 1
          IF (SIZE(lo)/=1) check_lo_status = 2
        END IF

        check_hi_status = -1
        IF (PRESENT(hi)) THEN
          check_hi_status = 1
          IF (SIZE(hi)/=1) check_hi_status = 2
        END IF

        IF (PRESENT(number_wanted)) THEN
          use_n_wanted = number_wanted
        ELSE
          use_n_wanted = SIZE(x)
        END IF

TRY_LOOP: DO i_try = 1, max_try

          IF (i_try>1) WRITE (*,*) 'Please try again from the beginning.'

          IF (PRESENT(format)) WRITE (*,format)

          CALL prompt

!     Read in value.

          READ (*,*,iostat=ios) (x(i),i=1,use_n_wanted)
          IF (ios/=0) THEN
            WRITE (*,*) &
              'Some entered value not accepted as a legal number.'
            CYCLE TRY_LOOP
          END IF

!     Check bound on each entry in X.

CHECK_VALUES: DO i = 1, use_n_wanted

!        Check lower bound if required.

90000       FORMAT ('The ',I3, &
              '''th value entered < the lower bound allowed.'/ &
              'The value entered was: ',G20.7/'The lower bound is: ', &
              G20.7)
90100       FORMAT ('The ',I3, &
              '''th value entered <= the lower bound allowed.'/ &
              'The value entered was: ',G20.7/'The lower bound is: ', &
              G20.7)

LO_CHECK:   IF (check_lo_status>0) THEN
              IF (check_lo_status==1) THEN
                IF (use_lo_eq_ok .AND. x(i)<lo(1)) THEN
                  WRITE (*,90000) i, x(i), lo(1)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_lo_eq_ok .AND. x(i)<=lo(1)) THEN
                  WRITE (*,90100) i, x(i), lo(1)
                  CYCLE TRY_LOOP
                END IF
              ELSE
                IF (use_lo_eq_ok .AND. x(i)<lo(i)) THEN
                  WRITE (*,90000) i, x(i), lo(i)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_lo_eq_ok .AND. x(i)<=lo(i)) THEN
                  WRITE (*,90100) i, x(i), lo(i)
                  CYCLE TRY_LOOP
                END IF
              END IF
            END IF LO_CHECK

!        Check upper bound if required.

90200       FORMAT ('The ',I3, &
              '''th value entered > the upper bound allowed.'/ &
              'The value entered was: ',G20.7/'The upper bound is: ', &
              G20.7)
90300       FORMAT ('The ',I3, &
              '''th value entered >= the upper bound allowed.'/ &
              'The value entered was: ',G20.7/'The upper bound is: ', &
              G20.7)

HI_CHECK:   IF (check_hi_status>0) THEN
              IF (check_hi_status==1) THEN
                IF (use_hi_eq_ok .AND. x(i)>hi(1)) THEN
                  WRITE (*,90200) i, x(i), hi(1)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_hi_eq_ok .AND. x(i)>=hi(1)) THEN
                  WRITE (*,90300) i, x(i), hi(1)
                  CYCLE TRY_LOOP
                END IF
              ELSE
                IF (use_hi_eq_ok .AND. x(i)>hi(i)) THEN
                  WRITE (*,90200) i, x(i), hi(i)
                  CYCLE TRY_LOOP
                ELSE IF ( .NOT. use_hi_eq_ok .AND. x(i)>=hi(i)) THEN
                  WRITE (*,90300) i, x(i), hi(i)
                  CYCLE TRY_LOOP
                END IF
              END IF
            END IF HI_CHECK

          END DO CHECK_VALUES

          IF (PRESENT(failed)) failed = .FALSE.
          RETURN

        END DO TRY_LOOP

!!! Failure if we reach here

        IF (PRESENT(failed)) THEN
          failed = .TRUE.
          RETURN
        ELSE
          WRITE (*,*) &
            'Too many failed attempts at entry in get_real_array. ABORT'
          STOP &
            'Too many failed attempts at entry in get_real_array. ABORT'
        END IF

      END SUBROUTINE get_real_array

!!!======================================================================

      SUBROUTINE get_double_scalar(format,lo,lo_eq_ok,hi,hi_eq_ok,x, &
          failed)
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: hi, lo
        REAL (dpkind), INTENT (OUT) :: x
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Local Scalars ..
        INTEGER :: ios, i_try
        LOGICAL :: hi_present, lo_present, use_hi_eq_ok, use_lo_eq_ok
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
!     Query user.
        use_lo_eq_ok = .TRUE.
        IF (PRESENT(lo_eq_ok)) use_lo_eq_ok = lo_eq_ok

        use_hi_eq_ok = .TRUE.
        IF (PRESENT(hi_eq_ok)) use_hi_eq_ok = hi_eq_ok

        lo_present = PRESENT(lo)
        hi_present = PRESENT(hi)

TRY_LOOP: DO i_try = 1, max_try

          IF (i_try>1) WRITE (*,*) 'Please try again.'

          IF (PRESENT(format)) THEN
            WRITE (*,format)
          END IF

          CALL prompt

!     Read in value.

          READ (*,*,iostat=ios) x
          IF (ios/=0) THEN
            WRITE (*,*) 'Entered value not accepted as a legal number.'
            CYCLE TRY_LOOP
          END IF

!     Check bound on each entry in X.

!        Check lower bound if required.

90000     FORMAT ('The value entered < the lower bound allowed.'/ &
            'The value entered was: ',G20.7/'The lower bound is: ',G20.7)
90100     FORMAT ('The value entered <= the lower bound allowed.'/ &
            'The value entered was: ',G20.7/'The lower bound is: ',G20.7)

LO_CHECK: IF (lo_present) THEN
            IF (use_lo_eq_ok .AND. x<lo) THEN
              WRITE (*,90000) x, lo
              CYCLE TRY_LOOP
            ELSE IF ( .NOT. use_lo_eq_ok .AND. x<=lo) THEN
              WRITE (*,90100) x, lo
              CYCLE TRY_LOOP
            END IF
          END IF LO_CHECK

!        Check upper bound if required.

90200     FORMAT ('The value entered > the upper bound allowed.'/ &
            'The value entered was: ',G20.7/'The upper bound is: ',G20.7)
90300     FORMAT ('The value entered >= the upper bound allowed.'/ &
            'The value entered was: ',G20.7/'The upper bound is: ',G20.7)

HI_CHECK: IF (hi_present) THEN
            IF (use_hi_eq_ok .AND. x>hi) THEN
              WRITE (*,90200) x, hi
              CYCLE TRY_LOOP
            ELSE IF ( .NOT. use_hi_eq_ok .AND. x>=hi) THEN
              WRITE (*,90300) x, hi
              CYCLE TRY_LOOP
            END IF
          END IF HI_CHECK

          IF (PRESENT(failed)) failed = .FALSE.
          RETURN

        END DO TRY_LOOP

!!! Failure if we reach here

        IF (PRESENT(failed)) THEN
          failed = .TRUE.
          RETURN
        ELSE
          WRITE (*,*) 'Too many failed attempts at entry in &
            &get_double_scalar. ABORT'
          STOP 'Too many failed attempts at entry in get_double_scalar. &
            &ABORT'
        END IF

      END SUBROUTINE get_double_scalar

!!!======================================================================

      SUBROUTINE get_real_scalar(format,lo,lo_eq_ok,hi,hi_eq_ok,x,failed)
! .. Scalar Arguments ..
        REAL (spkind), OPTIONAL, INTENT (IN) :: hi, lo
        REAL (spkind), INTENT (OUT) :: x
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Local Scalars ..
        INTEGER :: ios, i_try
        LOGICAL :: hi_present, lo_present, use_hi_eq_ok, use_lo_eq_ok
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
!     Query user.
        use_lo_eq_ok = .TRUE.
        IF (PRESENT(lo_eq_ok)) use_lo_eq_ok = lo_eq_ok

        use_hi_eq_ok = .TRUE.
        IF (PRESENT(hi_eq_ok)) use_hi_eq_ok = hi_eq_ok

        lo_present = PRESENT(lo)
        hi_present = PRESENT(hi)

TRY_LOOP: DO i_try = 1, max_try

          IF (i_try>1) WRITE (*,*) 'Please try again.'

          IF (PRESENT(format)) WRITE (*,format)

          CALL prompt

!     Read in value.

          READ (*,*,iostat=ios) x
          IF (ios/=0) THEN
            WRITE (*,*) 'Entered value not accepted as a legal number.'
            CYCLE TRY_LOOP
          END IF

!     Check bound on each entry in X.

!        Check lower bound if required.

90000     FORMAT ('The value entered < the lower bound allowed.'/ &
            'The value entered was: ',G20.7/'The lower bound is: ',G20.7)
90100     FORMAT ('The value entered <= the lower bound allowed.'/ &
            'The value entered was: ',G20.7/'The lower bound is: ',G20.7)

LO_CHECK: IF (lo_present) THEN
            IF (use_lo_eq_ok .AND. x<lo) THEN
              WRITE (*,90000) x, lo
              CYCLE TRY_LOOP
            ELSE IF ( .NOT. use_lo_eq_ok .AND. x<=lo) THEN
              WRITE (*,90100) x, lo
              CYCLE TRY_LOOP
            END IF
          END IF LO_CHECK

!        Check upper bound if required.

90200     FORMAT ('The value entered > the upper bound allowed.'/ &
            'The value entered was: ',G20.7/'The upper bound is: ',G20.7)
90300     FORMAT ('The value entered >= the upper bound allowed.'/ &
            'The value entered was: ',G20.7/'The upper bound is: ',G20.7)

HI_CHECK: IF (hi_present) THEN
            IF (use_hi_eq_ok .AND. x>hi) THEN
              WRITE (*,90200) x, hi
              CYCLE TRY_LOOP
            ELSE IF ( .NOT. use_hi_eq_ok .AND. x>=hi) THEN
              WRITE (*,90300) x, hi
              CYCLE TRY_LOOP
            END IF
          END IF HI_CHECK

          IF (PRESENT(failed)) failed = .FALSE.
          RETURN

        END DO TRY_LOOP

!!! Failure if we reach here

        IF (PRESENT(failed)) THEN
          failed = .TRUE.
          RETURN
        ELSE
          WRITE (*,*) 'Too many failed attempts at entry in &
            &get_real_scalar. ABORT'
          STOP 'Too many failed attempts at entry in get_real_scalar. &
            &ABORT'
        END IF

      END SUBROUTINE get_real_scalar

!!!======================================================================

      SUBROUTINE get_integer_scalar(format,lo,lo_eq_ok,hi,hi_eq_ok,x, &
          failed)
! .. Scalar Arguments ..
        INTEGER, OPTIONAL, INTENT (IN) :: hi, lo
        INTEGER, INTENT (OUT) :: x
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Local Scalars ..
        INTEGER :: ios, i_try
        LOGICAL :: hi_present, lo_present, use_hi_eq_ok, use_lo_eq_ok
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
!     Query user.
        use_lo_eq_ok = .TRUE.
        IF (PRESENT(lo_eq_ok)) use_lo_eq_ok = lo_eq_ok

        use_hi_eq_ok = .TRUE.
        IF (PRESENT(hi_eq_ok)) use_hi_eq_ok = hi_eq_ok

        lo_present = PRESENT(lo)

        hi_present = PRESENT(hi)


TRY_LOOP: DO i_try = 1, max_try

          IF (i_try>1) WRITE (*,*) 'Please try again.'

          IF (PRESENT(format)) WRITE (*,format)

          CALL prompt

!     Read in value.

          READ (*,*,iostat=ios) x
          IF (ios/=0) THEN
            WRITE (*,*) 'Entered value not accepted as a legal number.'
            CYCLE TRY_LOOP
          END IF

!     Check bound on each entry in X.

!        Check lower bound if required.

90000     FORMAT ('The value entered < the lower bound allowed.'/ &
            'The value entered was: ',I20/'The lower bound is: ',I20)
90100     FORMAT ('The value entered <= the lower bound allowed.'/ &
            'The value entered was: ',I20/'The lower bound is: ',I20)

LO_CHECK: IF (lo_present) THEN
            IF (use_lo_eq_ok .AND. x<lo) THEN
              WRITE (*,90000) x, lo
              CYCLE TRY_LOOP
            ELSE IF ( .NOT. use_lo_eq_ok .AND. x<=lo) THEN
              WRITE (*,90100) x, lo
              CYCLE TRY_LOOP
            END IF
          END IF LO_CHECK

!        Check upper bound if required.

90200     FORMAT ('The value entered > the upper bound allowed.'/ &
            'The value entered was: ',I20/'The upper bound is: ',I20)
90300     FORMAT ('The value entered >= the upper bound allowed.'/ &
            'The value entered was: ',I20/'The upper bound is: ',I20)

HI_CHECK: IF (hi_present) THEN
            IF (use_hi_eq_ok .AND. x>hi) THEN
              WRITE (*,90200) x, hi
              CYCLE TRY_LOOP
            ELSE IF ( .NOT. use_hi_eq_ok .AND. x>=hi) THEN
              WRITE (*,90300) x, hi
              CYCLE TRY_LOOP
            END IF
          END IF HI_CHECK

          IF (PRESENT(failed)) failed = .FALSE.
          RETURN

        END DO TRY_LOOP

!!! Failure if we reach here

        IF (PRESENT(failed)) THEN
          failed = .TRUE.
          RETURN
        ELSE
          WRITE (*,*) 'Too many failed attempts at entry in &
            &get_integer_scalar. ABORT'
          STOP 'Too many failed attempts at entry in get_integer_scalar. &
            &ABORT'
        END IF

      END SUBROUTINE get_integer_scalar

!*********************************************************************

      SUBROUTINE get_list_double(format,list,n_list,max_n_list,lo, &
          lo_eq_ok,hi,hi_eq_ok,failed)
!-----------------------------------------------------------------------
!     Name: get_list_double (Double precision LIST)
!     Description:
!     Creates or modifies a list of double precision values via a
!     dialogue with the user.  The list values input by the user may
!     optionally be required to be greater than or equal to (or
!     strictly greater than) a supplied lower bound, less than or
!     equal to (or strictly less than) a supplied upper bound, or
!     both.  See below for details.  Will optionally STOP or return
!     with an error indicator if the user makes too many mistakes
!     while engaged in the dialogue.
!     Arguments:
!     FORMAT      -> If present, a character string containing a
!                    Fortran format.  The format will be printed to
!                    standard output (*) to prompt the user.  If the
!                    argument is missing, no prompt message will be
!                    printed.
!                    CHARACTER(*), OPTIONAL, INTENT(IN) :: format
!     LIST       <-> The list of values created or modified by the
!                    user.
!                    REAL(dpkind), INTENT(INOUT) :: list(:)
!     N_LIST     <-> Length of the list of numbers in LIST.
!                    INTEGER, INTENT(INOUT) :: n_list
!     MAX_N_LIST  -> If present, the maximum length of the list
!                    created or modified by the user.  If the argument
!                    is missing, the maximum will be the size of LIST.
!                    INTEGER, OPTIONAL, INTENT(IN) :: max_n_list
!     LO          -> If present, the lower bound on list values which
!                    may be entered by the user.  If absent, no lower
!                    bound checking is performed.  If checking is
!                    performed, the input of values lower than LO
!                    counts as a failed attempt to obtain list input
!                    from the user.  For the case of values which
!                    equal LO, see the LO_EQ_OK entry, and the Details
!                    and Notes sections below.
!                    REAL (dpkind), OPTIONAL, INTENT (IN) :: lo
!     LO_EQ_OK    -> If LO is present, this value if present
!                    determines if a value equal to LO is acceptable
!                    as a list value input by the user.  The default
!                    value for LO_EQ_OK if missing is .TRUE., by which
!                    it is meant that if LO is present and LO_EQ_OK is
!                    absent, the default behavior is to allow list
!                    entries which equal LO.  If LO is present and
!                    LO_EQ_OK is present and .TRUE., list entries
!                    equal to LO are allowed.  If LO is present and
!                    LO_EQ_OK is present and .FALSE., list entries
!                    equal to LO are not allowed.  If LO is not
!                    present, the presence and value of LO_EQ_OK has
!                    no effect.  For more information on exactly what
!                    equality means in this case, see the Notes
!                    section below.
!                    LOGICAL, OPTIONAL, INTENT (IN) :: lo_eq_ok
!     HI          -> If present, the upper bound on list values which
!                    may be entered by the user.  If absent, no upper
!                    bound checking is performed.  If checking is
!                    performed, the input of values greater than HI
!                    counts as a failed attempt to obtain list input
!                    from the user.  For the case of values which
!                    equal HI, see the HI_EQ_OK entry, and the Details
!                    and Notes sections below.
!                    REAL (dpkind), OPTIONAL, INTENT (IN) :: hi
!     HI_EQ_OK    -> If HI is present, this value if present
!                    determines if a value equal to HI is acceptable
!                    as a list value input by the user.  The default
!                    value for HI_EQ_OK if missing is .TRUE., by which
!                    it is meant that if HI is present and HI_EQ_OK is
!                    absent, the default behavior is to allow list
!                    entries which equal HI.  If HI is present and
!                    HI_EQ_OK is present and .TRUE., list entries
!                    equal to HI are allowed.  If HI is present and
!                    HI_EQ_OK is present and .FALSE., list entries
!                    equal to HI are not allowed.  If HI is not
!                    present, the presence and value of HI_EQ_OK has
!                    no effect.  For more information on exactly what
!                    equality means in this case, see the Notes
!                    section below.
!                    LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok
!     FAILED     <-  If present and .TRUE. on exit, MAX_FAILURES
!                    attempts failed to obtain list input from the
!                    user.  If present and .FALSE. on exit, list was
!                    obtained successfully.  If argument is missing,
!                    this subroutine will abort (STOP) with an error
!                    message if MAX_FAILURES attempts fail to obtain
!                    list input from the user.  Note that regardless
!                    of the presence of this optional argument, the
!                    user is notified of errors by messages written to
!                    the * output unit.
!                    LOGICAL, OPTIONAL, INTENT(OUT) :: failed
!     Details:
!     Conducts a dialog with the user to create or modify a list of
!     values in LIST.  An arbitrary number of actions can be taken to
!     make the list.  Each action is one of the following:
!          (1) Enter a user specified number of values
!          (2) Enter low value, high value, and number of intervals.
!              The values added will be equally spaced between the low
!              and high values specified.  Note that there will be
!              number of intervals plus one values added because both
!              the low and high values are on the list.
!          (3) Like (2), but values are equally spaced on a
!              logarithmic scale.
!          (4) Print the list.
!          (5) Delete a single element of the list.
!          (6) Delete consecutive elements of the list.
!          (7) Sort the list of values in ascending order and delete
!              duplicates.
!          (8) Quit modifying the list.
!     The values input by the user may be bounded from below by the
!     optional argument LO, and from above by the optional argument
!     HI, according to the presence and value of the corresponding
!     optional logical argument LO_EQ_OK or HI_EQ_OK.
!     If LO is absent no lower bound checking is performed.
!     If LO is present and LO_EQ_OK is present and .FALSE., the user
!     input values must be stricly greater than LO.  Otherwise the
!     user input must be greater than or equal to LO.  I.e. the
!     default value for LO_EQ_OK if it is missing is .TRUE., meaning
!     that user inputs may equal the lower bound supplied by LO.
!     If HI is absent no upper bound checking is performed.
!     If HI is present and HI_EQ_OK is present and .FALSE., the user
!     input values must be stricly less than HI.  Otherwise the user
!     input must be less than or equal to HI.  I.e. the default value
!     for HI_EQ_OK if it is missing is .TRUE., meaning that user
!     inputs may equal the upper bound supplied by HI.
!     Notes:
!     The user is allowed to quit without adding any values to the
!     list.
!     The user is notified by a message printed to the * output unit
!     when a mistake is made, regardless of the presence of the
!     optional argument FAILED.
!     For the purposes of eliminating duplicates in action (7)
!     described above, duplicates are defined as inputs to the
!     function APPROX_EQUAL which yield a .TRUE. return value.  This
!     subroutine CONTAINS the FUNCTION APPROX_EQUAL; see below for
!     more details.
!     For the purposes of bounds checking, equality is ordinary
!     floating point equality.  I.e. X is less than or equal to HI if
!     X <= HI evaluates to .TRUE. in your implementation of Fortran.
!     Note that this is NOT the approximate equality discussed above
!     and the appropriate care should be taken when using this
!     feature.
!-----------------------------------------------------------------------
! .. Use Statements ..
        USE biomath_sort_mod, ONLY : sort_list
! ..
! .. Parameters ..
! the allowed number of mistakes before an abort
! INTEGER, PARAMETER :: max_failures = 3
! number of lines in terminal or window; 0 = no limit
! integer, parameter :: window_lines = 24
! number of lines used by hold
! integer, parameter :: hold_lines = 3
! lines used by list; 0 = no limit
! integer, parameter :: list_lines = max(0,window_lines-hold_lines)
! maximum length of character representation of a default integer
! INTEGER, PARAMETER :: integer_length = 10
! the format corresponding to integer_length
! character(len=*), parameter :: integer_format = '(I10)'
! the main menu as a Fortran FORMAT string
        INTEGER, PARAMETER :: hold_lines = 3, integer_length = 10, &
          max_failures = 3, window_lines = 24
        INTEGER, PARAMETER :: list_lines = MAX(0,window_lines-hold_lines)
        CHARACTER (*), PARAMETER :: integer_format = '(I10)'
        CHARACTER (*), PARAMETER :: menu_format = '(/5X,''Select &
          &an option:''/5X,''   1 - ADD individually specified &
          &values to the list''/5X,''   2 - ADD equally spaced &
          &values to the list''/5X,''   3 - ADD logarithmically &
          &spaced values to the list''/5X,''   4 - PRINT the &
          &list''/5X,''   5 - DELETE a single element of the &
          &list''/5X,''   6 - DELETE consecutive elements of &
          &the list''/5X,''   7 - SORT the list in ascending &
          &order and eliminate duplicate values''/5X,''   8 &
          &- QUIT modifying the list'')'
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL, INTENT (IN) :: hi, lo
        INTEGER, OPTIONAL, INTENT (IN) :: max_n_list
        INTEGER, INTENT (INOUT) :: n_list
        LOGICAL, OPTIONAL, INTENT (OUT) :: failed
        LOGICAL, OPTIONAL, INTENT (IN) :: hi_eq_ok, lo_eq_ok
        CHARACTER (*), OPTIONAL, INTENT (IN) :: format
! ..
! .. Array Arguments ..
        REAL (dpkind), INTENT (INOUT) :: list(:)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: xhi, xint, xlo
        INTEGER :: i, idiff, idum, ihi, ilo, imenu, local_max_n_list, &
          nin, nleft, n_failures
        LOGICAL :: check_hi, check_lo, local_failed, local_hi_eq_ok, &
          local_lo_eq_ok
        CHARACTER (integer_length) :: cint, cmax
        CHARACTER (120) :: string
! ..
! .. Local Arrays ..
        REAL (dpkind) :: temp(3)
        INTEGER :: itemp(2)
! ..
! .. Intrinsic Functions ..
        INTRINSIC ADJUSTL, EXP, INT, LOG, MAX, MIN, MOD, PRESENT, REAL, &
          SIZE, TRIM
! ..
!***********************************************************************
!     Initalize Variables
!***********************************************************************
        IF (PRESENT(max_n_list)) THEN
          local_max_n_list = max_n_list
        ELSE
          local_max_n_list = SIZE(list)
        END IF

        check_lo = PRESENT(lo)

        IF (PRESENT(lo_eq_ok)) THEN
          local_lo_eq_ok = lo_eq_ok
        ELSE
          local_lo_eq_ok = .TRUE.
        END IF

        check_hi = PRESENT(hi)

        IF (PRESENT(hi_eq_ok)) THEN
          local_hi_eq_ok = hi_eq_ok
        ELSE
          local_hi_eq_ok = .TRUE.
        END IF

        WRITE (cmax,integer_format) local_max_n_list
        cmax = ADJUSTL(cmax)

        n_failures = 0

10      CONTINUE

!***********************************************************************
!     Outer loop allows user to select an option for modifying a list
!***********************************************************************

        IF (PRESENT(failed)) failed = .FALSE.

!     Put up the query statement

        IF (PRESENT(format)) WRITE (*,format)

!     Put put the statement of how many items are in the list
!     and the maximum items allowed.

        WRITE (cint,integer_format) n_list
        cint = ADJUSTL(cint)

        WRITE (*,90000) ' '
        WRITE (*,90200) 'The list currently contains ', TRIM(cint), &
          ' values.'
        WRITE (*,90200) 'The maximum number allowed is ', TRIM(cmax), '.'

!     Put up main menu and read in selection

        CALL get_numbers(format=menu_format,lo=1,hi=8,x=imenu, &
          failed=failed)

!     If FAILED present and .TRUE., print error message and return.
!     Note: if FAILED not present, above call STOPs on failure.

        IF (PRESENT(failed)) THEN
          IF (failed) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Too many bad entries.'
            RETURN
          END IF
        END IF

!     Select an action based on the menu selection

        IF (imenu==1) THEN
!***********************************************************************
!     Enter a specified number of values.
!***********************************************************************

!        Determine the maximum number of values that the user
!        can input.

          nleft = local_max_n_list - n_list

!        If there is no space left in the list echo an error message
!        and return to the main menu.

          IF (nleft<1) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'The list is full.'
            WRITE (*,90000) &
              'Attempted to add too many values to the list.'
            GO TO 20
          END IF

!        Ask the user how many values are being input.

          WRITE (*,90000) ' '

          WRITE (cint,integer_format) nleft
          cint = ADJUSTL(cint)

          WRITE (string,90200) &
            '(1X,''How many values will be entered? [0,', TRIM(cint), &
            ']'')'
          CALL get_numbers(format=string,lo=0,hi=nleft,x=nin, &
            failed=local_failed)

!     If LOCAL_FAILED is .TRUE., print error message and branch to
!     error handler

          IF (local_failed) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Too many bad entries.'
            GO TO 20
          END IF

!        If the user enters 0 return to the main menu.

          IF (nin==0) GO TO 10

!        Ask the user to input the requested number of values.

          WRITE (*,90000) ' '
          IF (nin==1) THEN
            string = '(1X,''Enter 1 value.'')'
          ELSE

            WRITE (cint,integer_format) nin
            cint = ADJUSTL(cint)

            WRITE (string,90200) '(1X,''Enter ', TRIM(cint), &
              ' values.'')'
          END IF
          IF (check_lo) THEN
            IF (check_hi) THEN
              CALL get_numbers(format=string,lo=(/lo/),lo_eq_ok=lo_eq_ok, &
                hi=(/hi/),hi_eq_ok=hi_eq_ok,x=list(n_list+1:), &
                number_wanted=nin,failed=local_failed)
            ELSE
              CALL get_numbers(format=string,lo=(/lo/),lo_eq_ok=lo_eq_ok, &
                x=list(n_list+1:),number_wanted=nin,failed=local_failed)
            END IF
          ELSE
            IF (check_hi) THEN
              CALL get_numbers(format=string,hi=(/hi/),hi_eq_ok=hi_eq_ok, &
                x=list(n_list+1:),number_wanted=nin,failed=local_failed)
            ELSE
              CALL get_numbers(format=string,x=list(n_list+1:), &
                number_wanted=nin,failed=local_failed)
            END IF
          END IF

!     If LOCAL_FAILED is .TRUE., print error message and branch to
!     error handler

          IF (local_failed) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Too many bad entries.'
            GO TO 20
          END IF

!        Update the length of the list

          n_list = n_list + nin
!***********************************************************************
!     Add values specified between bounds with either linear
!     or logrithmic spacing.
!***********************************************************************
        ELSE IF (imenu==2 .OR. imenu==3) THEN

!        Determine the maximum number of values that the user
!        can input.

          nleft = local_max_n_list - n_list

!        If there is no space left in the list echo an error message
!        and return to the main menu.

          IF (nleft<1) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'The list is full.'
            WRITE (*,90000) &
              'Attempted to add too many values to the list.'
            GO TO 20
          END IF

!        Get both bounds and the number of
!        intervals from the user.

          WRITE (*,90000) ' '
          string = '(1X,''Enter beginning bound, ending bound &
            &and number of intervals.'')'
          CALL get_numbers(format=string,x=temp,failed=local_failed)

!     If LOCAL_FAILED is .TRUE., print error message and branch to
!     error handler

          IF (local_failed) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Too many bad entries.'
            GO TO 20
          END IF

!        Define xlo, xhi and nin
!        NOTE that nin, the number of points being added,
!        is the number of intervals + 1.

          xlo = temp(1)
          xhi = temp(2)
          nin = INT(temp(3))

!        Check that NIN, the number of intervals is greater than or
!        equal to 1.

          IF (nin<1) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Number of intervals must be at least 1'
            GO TO 20
          END IF

!        Check that NIN+1, the number begin added, is less than NLEFT,
!        the number of positions left in the list.

          IF (nin+1>nleft) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Attempted to add too many values to list.'
            GO TO 20
          END IF

!        If Log spacing is being used check that XLO and XHI
!        are greater than 0.

          IF (imenu==3) THEN
            IF (xlo<=0.0E0_dpkind .OR. xhi<=0.0E0_dpkind) THEN
              WRITE (*,90000) ' '
              WRITE (*,90000) &
                'Beginning and ending bounds must be greater than 0.'
              GO TO 20
            END IF
          END IF

!        If bounds are being checked, check that XLO and XHI
!        are within the bounds.

          IF (check_lo) THEN
            IF (local_lo_eq_ok .AND. MIN(xlo,xhi)<lo) THEN
              WRITE (*,90000) ' '
              WRITE (*,90300) 'Beginning and ending bounds &
                &must be greater than or equal to ', lo, '.'
              GO TO 20
            ELSE IF ( .NOT. local_lo_eq_ok .AND. MIN(xlo,xhi)<=lo) THEN
              WRITE (*,90000) ' '
              WRITE (*,90300) &
                'Beginning and ending bounds must be greater than ', lo, &
                '.'
              GO TO 20
            END IF
          END IF
          IF (check_hi) THEN
            IF (local_hi_eq_ok .AND. MAX(xlo,xhi)>hi) THEN
              WRITE (*,90000) ' '
              WRITE (*,90300) 'Beginning and ending bounds &
                &must be less than or equal to ', hi, '.'
              GO TO 20
            ELSE IF ( .NOT. local_hi_eq_ok .AND. MAX(xlo,xhi)>=hi) THEN
              WRITE (*,90000) ' '
              WRITE (*,90300) &
                'Beginning and ending bounds must be less than ', hi, '.'
              GO TO 20
            END IF
          END IF

!        Create a linear spaced list.

          IF (imenu==2) THEN
            xint = (xhi-xlo)/REAL(nin,kind=dpkind)
            DO i = 1, nin + 1
              n_list = n_list + 1
              list(n_list) = xlo + REAL(i-1,kind=dpkind)*xint
            END DO

!        Create a log spaced list.

          ELSE IF (imenu==3) THEN
            xint = LOG(xhi/xlo)/REAL(nin,kind=dpkind)
            xlo = LOG(xlo)
            DO i = 1, nin + 1
              n_list = n_list + 1
              list(n_list) = EXP(xlo+REAL(i-1,kind=dpkind)*xint)
            END DO
          END IF
!***********************************************************************
!     Print the list.
!***********************************************************************
        ELSE IF (imenu==4) THEN
          IF (n_list>0) THEN
            DO i = 1, n_list
              WRITE (*,90400) i, list(i)
              IF (i==n_list .OR. list_lines==0 .OR. MOD(i,list_lines)/=0) &
                CYCLE
              CALL hold
            END DO
            IF (list_lines>0) THEN
              CALL hold
            END IF
          END IF
!***********************************************************************
!     Delete an single element or consecutive elements from the array.
!***********************************************************************
        ELSE IF (imenu==5 .OR. imenu==6) THEN

!        Check to see that there is at least 1 element in the list.

          IF (n_list<1) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'There are no elements in the list.'
            WRITE (*,90000) 'You must create a list before changing it.'
            GO TO 20
          END IF

!        Ask the user to enter the array position of elements to
!        be deleted.

          WRITE (cint,integer_format) n_list
          cint = ADJUSTL(cint)

          IF (imenu==5) THEN
            i = 1
            WRITE (string,90200) '(1X,''Enter position in array &
              &for value to be deleted. [1,', TRIM(cint), ']'')'
          ELSE IF (imenu==6) THEN
            i = 2
            WRITE (string,90200) '(1X,''Enter low and high &
              &positions in array for values to be deleted. &
              &[1,', TRIM(cint), ']'')'
          END IF

          WRITE (*,90000) ' '
          CALL get_numbers(format=string,lo=(/1/),hi=(/n_list/),x=itemp, &
            number_wanted=i,failed=local_failed)

!     If LOCAL_FAILED is .TRUE., print error message and branch to
!     error handler

          IF (local_failed) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'Too many bad entries.'
            GO TO 20
          END IF

!        Define ILO and IHI, the starting and stopping array positions
!        to be deleted.

          IF (imenu==5) THEN
            ilo = itemp(1)
            ihi = ilo
          ELSE IF (imenu==6) THEN
            IF (itemp(1)<=itemp(2)) THEN
              ilo = itemp(1)
              ihi = itemp(2)
            ELSE
              ilo = itemp(2)
              ihi = itemp(1)
            END IF
          END IF

!        Delete the selected positions from the list.

          idiff = (ihi-ilo) + 1
          i = ilo
          IF (n_list-idiff-ilo+1>0) THEN
            list(ilo:n_list-idiff) = list(ilo+idiff:n_list)
            i = n_list - idiff + 1
          END IF
          n_list = n_list - idiff
!***********************************************************************
!     Sort list in ascending order and eliminate duplicate values.
!***********************************************************************
        ELSE IF (imenu==7) THEN

!        Check to see that there is at least 1 element in the list.

          IF (n_list<1) THEN
            WRITE (*,90000) ' '
            WRITE (*,90000) 'There are no elements in the list.'
            WRITE (*,90000) 'You must create a list before sorting it.'
            GO TO 20
          END IF

!        Sort the list in ascending order.

          CALL sort_list(list,n_list)

!        Delete any duplicate values.

          idum = 0
          DO i = 2, n_list
            IF (approx_equal(list(i),list(i-idum-1))) THEN
              idum = idum + 1
            ELSE
              list(i-idum) = list(i)
            END IF
          END DO

          n_list = n_list - idum
!***********************************************************************
!     Quit if requested.
!***********************************************************************
        ELSE IF (imenu==8) THEN
          RETURN
        END IF
!***********************************************************************
!     Return to the main menu
!***********************************************************************
        GO TO 10
!***********************************************************************
!     Deal with an error condition.
!***********************************************************************
20      CONTINUE
        n_failures = n_failures + 1

!     If too many errors have occured:
!        if FAILED present, set FAILED to TRUE and return to the calling
!        program
!        if FAILED not present, WRITE abort message and STOP
!     Otherwise, return to the main menu.

        IF (n_failures>max_failures) THEN
          IF (PRESENT(failed)) THEN
            failed = .TRUE.
            RETURN
          ELSE
            WRITE (*,90000) 'Too many failed attempts at list &
              &creation or modification in get_list_double. &
              &ABORT'
            STOP ' Too many failed attempts at list creation &
              &or modification in get_list_double. ABORT'
          END IF
        ELSE
          WRITE (*,90000) 'Please try again.'
          GO TO 10
        END IF
!     ..
!     .. Format Statements ..
90000   FORMAT (1X,1A)
90100   FORMAT (1X,2A)
90200   FORMAT (1X,3A)
90300   FORMAT (1X,1A,1P,1E15.6,1A)
90400   FORMAT (1X,1I10,2X,1P,1E15.6)

        RETURN

      CONTAINS

        FUNCTION approx_equal(zx,zy)
! Returns .TRUE. if the arguments are approximately equal
! otherwise returns .FALSE.
! This version is a floating point comparison.
! Note that MIN_DENOM and EPSILON are related to the precision
! of the floating point KIND being compared, and should be changed
! if the KIND is changed.
! .. Function Return Value ..
          LOGICAL :: approx_equal
! ..
! .. Parameters ..
          REAL (dpkind), PARAMETER :: epsilon = 1.0E-14_dpkind
          REAL (dpkind), PARAMETER :: min_denom = 1.0E-100_dpkind
! ..
! .. Scalar Arguments ..
          REAL (dpkind) :: zx, zy
! ..
! .. Intrinsic Functions ..
          INTRINSIC ABS, MAX
! ..
          approx_equal = ABS(zx-zy)/MAX(ABS(zx)+ABS(zy),min_denom) < &
            epsilon

          RETURN

        END FUNCTION approx_equal

      END SUBROUTINE get_list_double

!*********************************************************************

! Description:
! The following routines are intended to help automate the printing of
! user messages to the current output unit.  It is intended to be used
! with the output of format.specs, a Perl script.  It should be included
! with  USE statement(s) visible to all the code which incorporates the
! output of format.specs.  Typically, this will only require the programmer
! to assign values to SUB_STRINGS in order to produce attractive output.
! NOTE: The charcter strings in  SUB_STRINGS are processed to double any
! single  quotation marks  found,  and   the substitution  is   adjusted
! accordingly.  This should  make everything run  smoothly provided that
! the  string to be  substituted appears  as  it should appear on output
! (i.e.  the  quotation marks are  not pre-doubled), and that the format
! strings use  single   quotes  to  enclose character   constant  format
! elements.  If these conditions are not  met, the output may have twice
! the single  quotation marks as desired, but  the format  string should
! not cause an I/O error.
! Module Parameters:
! Note: If any of these three parameters  are changed, the corresponding
! variable in format.specs should probably be changed as well.
! mx_format_length -- The  maximum character length of a  format string.
!                     INTEGER, PARAMETER :: mx_format_length
! mx_subs          -- The maximum number of substitutions in each format
!                     string processed.
!                     INTEGER, PARAMETER :: mx_subs
! mx_sub_length    -- The maximum character length of  each substitution
!                     string.
!                     INTEGER, PARAMETER :: mx_sub_length
! Module Arrays:
! sub_pos          -- The indices of the  current format string at which
!                     the  substitutions  are  to  be  made.   The  I'th
!                     substitution  is made to the sub-string indexed by
!                     ( sub_pos(2*I-1) : sub_pos(2*I) ).
!                     INTEGER, DIMENSION(2*mx_subs) :: sub_pos
! sub_strings      -- The sub-strings  to be substituted into the format
!                     string. The I'th substitution gets sub_strings(I).
!      CHARACTER(LEN=mx_sub_length), DIMENSION( mx_subs ) :: sub_strings
! Module Scalars:
! always_print     -- Logical indicator of  under what circumstances the
!                     current format string should be printed.
!                        .TRUE.  -- Print it regarless of print_level.
!                        .FALSE. -- Print it according to print_level.
!                     Help messages generally set the value to .FALSE.
!                     LOGICAL :: always_print
! message_format   -- The 'default'  format string  to be  printed.  The
!                     Perl  script  format.specs  assigns  unnamed  text
!                     blocks to this character string,  so most messages
!                     will  probably  use  it.   It  doesn't need  to be
!                     declared in code using this module,  and it can be
!                     printed with a call to print_message_format().
!                     CHARACTER(LEN=mx_format_length) :: message_format
! num_subs         -- The number of substitutions to be preformed on the
!                     current format string.
!                     INTEGER :: num_subs
! print_level      -- The (user specified)  level at which  messages are
!                     to be printed.
!                     Allowable values of print_level:
!                        1 -- always print help
!                        2 -- ask whether help should be printed
!                        3 -- never print help
!                     INTEGER, SAVE :: print_level
! format_printed   -- Logical  value   indicating  whether  the  message
!                     format  was  printed  or  not.   This  is  set  by
!                     print_format( format_string ),  and is intended to
!                     help the calling routine  control the scrolling of
!                     messages,  particularly  in the  case of  optional
!                     messages   (where  always_print   ==   .FALSE.  ).
!                     LOGICAL, SAVE :: format_printed
! clear_screen_before_print -- logical, intended to be set by programmer.
! Module Subroutines:
! print_message_format()
!     calls print_format( message_format ).
! print_format( format_string )
!   CHARACTER(LEN=*), INTENT(INOUT) :: format_string
!     Prints  out format_string   using   the  above variables.    Calls
!     edit_format(  format_string  ),  and  sets  format_printed.   This
!     subroutine uses the complib routines GTCUIO and QYNGT.
! edit_message_format()
!     calls edit_format( message_format ).
! edit_format( format_string )
!   CHARACTER(LEN=*), INTENT(INOUT) :: format_string
!     Prepares format_string for output  by performing the substitutions
!     indicated.  It is called by print_format.
! for CLEAR_SCREEN( ), windowing routines, etc.
!*********************************************************************

      SUBROUTINE edit_format(format_string)
!                     EDIT_FORMAT( FORMAT_STRING )
! Substitutes the   first   NUM_SUBS elements of    SUB_STRINGS  for the
! sub-strings in FORMAT_STRING defined by (SUB_POS(2*I-1):SUB_POS(2*I)).
! Intended  to prepare FORMAT_STRING for output,  i.e.  after  a call to
! EDIT_FORMAT fills in the desired values,
!                    WRITE( stdout, FORMAT_STRING)
! should produce the desired message.
! NOTE: This code  assumes   that FORMAT_STRING uses single  quotes   to
! delimit the character strings, and therefore doubles any single quotes
! found in the SUB_STRINGS, adjusting the positions of the substitutions
! and copying FORMAT_STRING appropriately.
! Argument:
! FORMAT_STRING <-> A character string intended to be the template for a
!                   message.  On return the positions bounded by SUB_POS
!                   are altered.
!                     CHARACTER(LEN=*), INTENT(INOUT) :: FORMAT_STRING
! Module symbols used (see above):
! NUM_SUBS, SUB_POS, SUB_STRINGS
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (INOUT) :: format_string
! ..
! .. Local Scalars ..
        INTEGER :: extra_fmt_chars, extra_sub_chars, i, jj
        CHARACTER (1) :: temp_char
        CHARACTER (mx_sub_length) :: temp_sub_string
! ..
        extra_fmt_chars = 0

SUBSTITUTE: DO i = 2, 2*num_subs, 2

! copy SUB_STRINGS( I/2 ) into TEMP_SUB_STRING,
! checking for single quotes along the way

          extra_sub_chars = 0

! loop over the significant characters in SUB_STRINGS( I/2 )

          DO jj = 1, sub_pos(i) - sub_pos(i-1) + 1

! copy this character

            temp_char = sub_strings(i/2) (jj:jj)

            temp_sub_string(jj+extra_sub_chars:jj+extra_sub_chars) &
              = temp_char

! if this character was a single quote, double it and keep track

            IF (temp_char=='''') THEN

! Advance TEMP_SUB_STRING index by one

              extra_sub_chars = extra_sub_chars + 1

! put in the duplicate single quote

              temp_sub_string(jj+extra_sub_chars:jj+extra_sub_chars) &
                = ''''

            END IF

          END DO

! Now TEMP_SUB_STRING has SUB_STRINGS( I/2 ) with dup. single quotes
! and EXTRA_SUB_CHARS is how many single quotes there were


! Shift remainder of FORMAT_STRING, if necessary
! NOTE: last EXTRA_SUB_CHARS characters will be lost
          IF (extra_sub_chars>0) format_string(sub_pos(i)+1+ &
            extra_fmt_chars+extra_sub_chars:) = format_string(sub_pos(i)+ &
            1+extra_fmt_chars:)

! Substitute

          format_string(sub_pos(i-1)+extra_fmt_chars: &
            sub_pos(i)+extra_fmt_chars+extra_sub_chars) = temp_sub_string

! Keep track of total extra characters added so far

          extra_fmt_chars = extra_fmt_chars + extra_sub_chars

        END DO SUBSTITUTE

        RETURN

      END SUBROUTINE edit_format

!*********************************************************************

      SUBROUTINE edit_message_format

        CALL edit_format(message_format)

        RETURN

      END SUBROUTINE edit_message_format

!*********************************************************************

      SUBROUTINE print_format(format_string,force,unit,unit_only)
! DMS Added optional arguments unit and unit_only to allow printing
!     in a report file.
!     unit --> Fortran unit associated with a report file.
!     unit_only --> if present, flags that printing will be
!                   performed only to the report file.
! .. Scalar Arguments ..
        INTEGER, OPTIONAL, INTENT (IN) :: unit
        LOGICAL, OPTIONAL, INTENT (IN) :: force, unit_only
        CHARACTER (*), INTENT (INOUT) :: format_string
! ..
! .. Local Scalars ..
        INTEGER :: output_unit
        LOGICAL :: force_print, qerr, qyes, want_help
! ..
! .. Parameters ..
        CHARACTER (*), PARAMETER :: print_question = '(1X,"Want &
          &next help message? (y/n)"/  1X,"If you don''t know &
          &what it is you probably want it.")'
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
! See whether or not to really print
        force_print = .FALSE.

        IF (PRESENT(force)) THEN
          force_print = force
        END IF

        IF (print_off .AND. ( .NOT. force_print)) RETURN

! Set output unit

        IF (PRESENT(unit)) THEN
          output_unit = unit
        ELSE
          output_unit = -1
        END IF

! Do substitutions if necessary

        IF (num_subs>0) THEN
          CALL edit_format(format_string)
        END IF

! See if this message is desired

        IF (( .NOT. always_print) .AND. (print_level==2)) THEN

          CALL get_yn(print_question,qyes,qerr)

          IF (qerr) THEN
            want_help = .TRUE.
          ELSE
            want_help = qyes
          END IF

        END IF

! Print message

        IF ((force_print .OR. always_print .OR. print_level==1) .OR. &
            (print_level==2 .AND. want_help)) THEN

          IF (clear_screen_before_print) CALL clear_screen

! Always display on standard output

          IF ( .NOT. PRESENT(unit_only)) THEN
            WRITE (stdout,format_string)
          END IF

! If session has a report file, write in that file also

          IF (output_unit>0) THEN
            WRITE (output_unit,format_string)
          END IF

          format_printed = .TRUE.

        ELSE
          format_printed = .FALSE.
        END IF

        RETURN

      END SUBROUTINE print_format

!*********************************************************************

      SUBROUTINE print_message_format(force,unit,unit_only)
! .. Scalar Arguments ..
        INTEGER, OPTIONAL :: unit
        LOGICAL, OPTIONAL :: force, unit_only
! ..
        CALL print_format(message_format,force,unit,unit_only)

        RETURN

      END SUBROUTINE print_message_format

!*********************************************************************

      SUBROUTINE write_array(x,astr)
! Writes the array elements in x, according to the format in astr.
! x --> array to be written
! astr --> the format to display the elements of array x
! Example:
! x(6)
! astr="4F11.6 3X 2F11.6"
! 4F11.6 - the first 4 elements will be written using the format F11.6
! 3X     - three spaces will be added
! 2F11.6 - the last two elements of x will be displayed
! Rules:
! The sum of the fields (4+2) MUST be equal to the length of x
! field_format: F11.6, E11.3, etc.
! separator: 3X
! Due to the lack of string functions in F90, we simplify a lot.
! Assume:
!   * separator:
!       first character is a digit (# of blanks to be written)
!       second character: X
!   * field_format:
!       first character: digit (number of elements to be written)
!       second character: F
!       last character: number of decimal positions
! .. Use Statements ..
        USE biomath_strings_mod
! ..
! .. Array Arguments ..
        REAL (dpkind), INTENT (IN) :: x(:)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: dval
        INTEGER :: i, icount, iend, indov, istart, ival, lcval, ndec, &
          nrep, wdt
        LOGICAL :: qstart
        CHARACTER (2) :: curtyp
        CHARACTER (6) :: cval
        CHARACTER (10) :: fmt, i_fmt
        CHARACTER (79) :: out_string
! ..
! .. Scalar Arguments ..
        CHARACTER (80), INTENT (IN) :: astr
! ..
! .. Intrinsic Functions ..
        INTRINSIC ABS, TRIM
! ..
        qstart = .TRUE.
        istart = 1
        iend = 0
        icount = 0
        out_string(1:79) = ' '

        DO WHILE (qlex(astr,qstart,curtyp,cval,lcval,ival,dval,indov))

! Suppose the FORMAT is '4F10.6', that is:
! REPETITION factor: nrep=4 (between 1..9)
! Field WIDTH:       wdt=10 (>= 10)
! DECIMAL places:    ndec=6 (between 3 and 6)

          SELECT CASE (cval(2:2))
          CASE ('F')

! Extract REPETITION factor
            WRITE (curtyp,'(A1)') cval(1:1)
            READ (curtyp,'(I1)') nrep

! Extract field WIDTH
            WRITE (curtyp,'(A2)') cval(3:4)
            READ (curtyp,'(I2)') wdt

! Extract number of DECIMAL places
            WRITE (curtyp,'(A1)') cval(6:6)
            READ (curtyp,'(I1)') ndec

            fmt = TRIM('('//cval(2:6))
            fmt(7:7) = ')'

! Write the fields into the string

            DO i = 1, nrep
              icount = icount + 1
              dval = x(icount)

! The format depends on the value to be displayed

              i_fmt = fmt

              IF (ABS(dval)<1.0E-3 .OR. ABS(dval)>1.0E3) THEN
                i_fmt(2:2) = 'E'

! Number of decimal places is decresed by 3

                IF (ndec>4) ndec = ndec - 3

                WRITE (i_fmt(6:6),'(I1)') ndec

              ELSE
                i_fmt(2:2) = 'F'
              END IF

              iend = istart + wdt - 1

              WRITE (out_string(istart:iend),i_fmt) dval

              istart = istart + wdt
            END DO

          CASE ('X')
            WRITE (curtyp,'(A1)') cval(1:1)
            READ (curtyp,'(I1)') wdt
            istart = istart + wdt

          CASE DEFAULT
            PRINT *, 'INTERNAL ERROR IN WRITE_ARRAY!'
            STOP 'INTERNAL ERROR IN WRITE_ARRAY!'

          END SELECT

        END DO

        WRITE (stdout,*) out_string

! If a report file was opened, write to it

        IF (report_unit>0) THEN
          WRITE (report_unit,*) out_string
        END IF

      END SUBROUTINE write_array

!*********************************************************************

      SUBROUTINE write_message(msg)
! Displays the message contained in the string "msg" on standard output.
! If the application has an associated report file,
! writes the message in that file also.
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (IN) :: msg
! ..
! .. Intrinsic Functions ..
        INTRINSIC TRIM
! ..
        WRITE (stdout,'(A)') TRIM(msg)

        IF (report_unit>0) THEN
          WRITE (report_unit,'(A)') TRIM(msg)
        END IF

        RETURN

      END SUBROUTINE write_message

!*********************************************************************

      SUBROUTINE write_error(error_message)
! Displays the message contained in the string "error_message"
! on standard output.
! If the application has an associated report file,
! writes the message in that file also.
! Then, STOPs the application.
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (IN) :: error_message
! ..
! .. Intrinsic Functions ..
        INTRINSIC TRIM
! ..
        WRITE (stdout,'(A)') TRIM(error_message)

        IF (report_unit>0) THEN
          WRITE (report_unit,'(A)') TRIM(error_message)
        END IF

        STOP 'ERROR in stattab'

      END SUBROUTINE write_error

!*********************************************************************

    END MODULE biomath_interface_mod
