    MODULE biomath_sort_mod
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: sort_list
! ..
! .. Generic Interface Blocks ..
      INTERFACE sort_list
        MODULE PROCEDURE sort_list_char, sort_list_dbl, &
          sort_list_integer, sort_list_single
      END INTERFACE
! ..
    CONTAINS

!*********************************************************************

      SUBROUTINE sort_list_char(a,ncol,a_gt_b)
!     *** QUICKSORT
! .. Parameters ..
        INTEGER, PARAMETER :: chgsrt = 10, mxlen = 256, stckmx = 50
! ..
! .. Scalar Arguments ..
        INTEGER :: ncol
! ..
! .. Array Arguments ..
        CHARACTER (*) :: a(ncol)
! ..
! .. Local Scalars ..
        INTEGER :: i, icol, j, l, part, r, stckct, swap1, swap2
        LOGICAL :: q, qdone
        CHARACTER (mxlen) :: temp
! ..
! .. Local Arrays ..
        INTEGER :: stck(2,stckmx)
! ..
! .. Intrinsic Functions ..
        INTRINSIC MAX, MIN
! ..
! .. Non-Generic Interface Blocks ..
        INTERFACE
          FUNCTION a_gt_b(a,b)
! .. Function Return Value ..
            LOGICAL :: a_gt_b
! ..
! .. Scalar Arguments ..
            CHARACTER (*) :: a, b
! ..
          END FUNCTION a_gt_b
        END INTERFACE
! ..
        OPTIONAL :: a_gt_b

        stckct = 0
        l = 1
        r = ncol
        qdone = (r-l) < chgsrt

!     *** 'MEDIAN-OF-THREE' MODIFICATION.  SEE 'OVERVIEW' ABOVE

        DO WHILE ( .NOT. (qdone))
          swap1 = (l+r)/2
          swap2 = l + 1
          CALL swap_elements

          IF (qgr(a(l+1),a(r))) THEN
            swap1 = l + 1
            swap2 = r
            CALL swap_elements
          END IF

          IF (qgr(a(l),a(r))) THEN
            swap1 = l
            swap2 = r
            CALL swap_elements
          END IF

          IF (qgr(a(l+1),a(l))) THEN
            swap1 = l + 1
            swap2 = l
            CALL swap_elements
          END IF

!     *** PARTITIONING OF SUBFILE

          i = l + 1
          j = r
          part = l
          DO

            DO
              i = i + 1
              IF ( .NOT. (qgr(a(part),a(i)))) EXIT
            END DO

            DO
              j = j - 1
              IF ( .NOT. (qgr(a(j),a(part)))) EXIT
            END DO

            IF (j<i) THEN
              EXIT
            ELSE
              swap1 = i
              swap2 = j

              CALL swap_elements
            END IF

          END DO

          swap1 = l
          swap2 = j

          CALL swap_elements


!     *** RECURSION STEP

!     *** DEFER FURTHER SORTING OF BOTH SUBFILES TO INSERTION SORT

          IF (MAX((j-l),(r-i+1))<=chgsrt) THEN
            IF (stckct==0) THEN
              qdone = .TRUE.
!     *** POP STACK AND CONTINUE QUICKSORT
            ELSE
              l = stck(1,stckct)
              r = stck(2,stckct)
              stckct = stckct - 1
            END IF

!     *** CONTINUE QUICKSORT ON AT LEAST ONE SUBFILE

          ELSE
!     *** DEFER SMALL TO INSERTION SORT, CONTINUE QUICKSORT ON LARGE

            IF (MIN((j-l),(r-i+1))<=chgsrt) THEN
!     *** LEFT SUBFILE IS LARGE ONE
!     L = L
              IF ((j-l)>=(r-i+1)) THEN
                r = j - 1
!     *** RIGHT SUBFILE IS LARGE ONE
              ELSE
                l = i
!     R = R
              END IF

!     *** CONTINUE QUICKSORT ON BOTH SUBFILES

            ELSE
!     *** STACK IS FULL
              IF (stckct>=stckmx) THEN
                STOP ' STACK OVERFLOW IN SORT_LIST_CHAR'
              END IF

!     *** PUSH LARGE SUBFILE ONTO STACK, CONTINUE QUICKSORT WITH SMALL
              stckct = stckct + 1
!     *** LEFT SUBFILE IS LARGE ONE

              IF ((j-l)>=(r-i+1)) THEN
                stck(1,stckct) = l
                stck(2,stckct) = j - 1
                l = i
!     R = R

!     *** RIGHT SUBFILE IS LARGE ONE

              ELSE
                stck(1,stckct) = i
                stck(2,stckct) = r
!     L = L
                r = j - 1

              END IF

            END IF

          END IF


        END DO


!     *** INSERTION SORT

        DO i = (ncol-1), 1, -1
!     *** FOR EACH POSITION PRECEDING THE LAST POSITION ...

!     ***
          IF (qgr(a(i),a(i+1))) THEN
            j = i + 1

!     *** A(I) Needs to be moved further out

            DO
              j = j + 1

              IF (j>ncol) THEN
                q = .FALSE.

              ELSE
                q = qgr(a(i),a(j))

              END IF

              IF ( .NOT. (q)) EXIT

            END DO


            temp = a(i)
            DO icol = i, j - 2
              a(icol) = a(icol+1)
            END DO
            a(j-1) = temp

          END IF

        END DO

        RETURN

      CONTAINS

!.....................................................................

        FUNCTION qgr(a,b)
! .. Function Return Value ..
          LOGICAL :: qgr
! ..
! .. Scalar Arguments ..
          CHARACTER (*), INTENT (IN) :: a, b
! ..
! .. Intrinsic Functions ..
          INTRINSIC PRESENT
! ..
          IF (PRESENT(a_gt_b)) THEN
            qgr = a_gt_b(a,b)
          ELSE
            qgr = a >= b
          END IF

          RETURN

        END FUNCTION qgr

!.....................................................................

        SUBROUTINE swap_elements
! .. Local Scalars ..
          CHARACTER (mxlen) :: temp
! ..
          temp = a(swap1)
          a(swap1) = a(swap2)
          a(swap2) = temp

        END SUBROUTINE swap_elements

      END SUBROUTINE sort_list_char

!*********************************************************************

      SUBROUTINE sort_list_dbl(a,ncol,a_gt_b)
!     *** QUICKSORT
! .. Parameters ..
        INTEGER, PARAMETER :: chgsrt = 10, stckmx = 50
! ..
! .. Scalar Arguments ..
        INTEGER :: ncol
! ..
! .. Array Arguments ..
        REAL (dpkind) :: a(ncol)
! ..
! .. Local Scalars ..
        REAL (dpkind) :: temp
        INTEGER :: i, icol, j, l, part, r, stckct, swap1, swap2
        LOGICAL :: q, qdone
! ..
! .. Local Arrays ..
        INTEGER :: stck(2,stckmx)
! ..
! .. Intrinsic Functions ..
        INTRINSIC MAX, MIN
! ..
! .. Non-Generic Interface Blocks ..
        INTERFACE
          FUNCTION a_gt_b(a,b)
! .. Use Statements ..
            USE biomath_constants_mod
! ..
! .. Function Return Value ..
            LOGICAL :: a_gt_b
! ..
! .. Scalar Arguments ..
            REAL (dpkind) :: a, b
! ..
          END FUNCTION a_gt_b
        END INTERFACE
! ..
        OPTIONAL :: a_gt_b

        stckct = 0
        l = 1
        r = ncol
        qdone = (r-l) < chgsrt

!     *** 'MEDIAN-OF-THREE' MODIFICATION.  SEE 'OVERVIEW' ABOVE

        DO WHILE ( .NOT. (qdone))
          swap1 = (l+r)/2
          swap2 = l + 1

          CALL swap_elements

          IF (qgr(a(l+1),a(r))) THEN
            swap1 = l + 1
            swap2 = r

            CALL swap_elements
          END IF

          IF (qgr(a(l),a(r))) THEN
            swap1 = l
            swap2 = r
            CALL swap_elements
          END IF

          IF (qgr(a(l+1),a(l))) THEN
            swap1 = l + 1
            swap2 = l

            CALL swap_elements
          END IF

!     *** PARTITIONING OF SUBFILE

          i = l + 1
          j = r
          part = l

          DO

            DO
              i = i + 1
              IF ( .NOT. (qgr(a(part),a(i)))) EXIT
            END DO

            DO
              j = j - 1
              IF ( .NOT. (qgr(a(j),a(part)))) EXIT
            END DO

            IF (j<i) THEN
              EXIT
            ELSE
              swap1 = i
              swap2 = j

              CALL swap_elements
            END IF

          END DO

          swap1 = l
          swap2 = j

          CALL swap_elements

!     *** RECURSION STEP

!     *** DEFER FURTHER SORTING OF BOTH SUBFILES TO INSERTION SORT

          IF (MAX((j-l),(r-i+1))<=chgsrt) THEN
            IF (stckct==0) THEN
              qdone = .TRUE.
!     *** POP STACK AND CONTINUE QUICKSORT
            ELSE
              l = stck(1,stckct)
              r = stck(2,stckct)
              stckct = stckct - 1
            END IF

!     *** CONTINUE QUICKSORT ON AT LEAST ONE SUBFILE

          ELSE
!     *** DEFER SMALL TO INSERTION SORT, CONTINUE QUICKSORT ON LARGE

            IF (MIN((j-l),(r-i+1))<=chgsrt) THEN
!     *** LEFT SUBFILE IS LARGE ONE
!     L = L
              IF ((j-l)>=(r-i+1)) THEN
                r = j - 1
!     *** RIGHT SUBFILE IS LARGE ONE
              ELSE
                l = i
!     R = R
              END IF

!     *** CONTINUE QUICKSORT ON BOTH SUBFILES

            ELSE
!     *** STACK IS FULL
              IF (stckct>=stckmx) THEN
                STOP ' STACK OVERFLOW IN SORT_LIST_DBL'
              END IF

!     *** PUSH LARGE SUBFILE ONTO STACK, CONTINUE QUICKSORT WITH SMALL
              stckct = stckct + 1
!     *** LEFT SUBFILE IS LARGE ONE

              IF ((j-l)>=(r-i+1)) THEN
                stck(1,stckct) = l
                stck(2,stckct) = j - 1
                l = i
!     R = R

!     *** RIGHT SUBFILE IS LARGE ONE

              ELSE
                stck(1,stckct) = i
                stck(2,stckct) = r
!     L = L
                r = j - 1

              END IF

            END IF

          END IF


        END DO

!     *** INSERTION SORT

        DO i = (ncol-1), 1, -1
!     *** FOR EACH POSITION PRECEDING THE LAST POSITION ...

!     ***
          IF (qgr(a(i),a(i+1))) THEN
            j = i + 1

!     *** A(I) Needs to be moved further out

            DO
              j = j + 1

              IF (j>ncol) THEN
                q = .FALSE.

              ELSE
                q = qgr(a(i),a(j))

              END IF

              IF ( .NOT. (q)) EXIT

            END DO


            temp = a(i)
            DO icol = i, j - 2
              a(icol) = a(icol+1)
            END DO
            a(j-1) = temp

          END IF

        END DO
        RETURN

      CONTAINS

!.....................................................................

        FUNCTION qgr(a,b)
! .. Function Return Value ..
          LOGICAL :: qgr
! ..
! .. Scalar Arguments ..
          REAL (dpkind), INTENT (IN) :: a, b
! ..
! .. Intrinsic Functions ..
          INTRINSIC PRESENT
! ..
          IF (PRESENT(a_gt_b)) THEN
            qgr = a_gt_b(a,b)
          ELSE
            qgr = a >= b
          END IF

          RETURN

        END FUNCTION qgr

!.....................................................................

        SUBROUTINE swap_elements
! .. Local Scalars ..
          REAL (dpkind) :: temp
! ..
          temp = a(swap1)
          a(swap1) = a(swap2)
          a(swap2) = temp

        END SUBROUTINE swap_elements

!.....................................................................

      END SUBROUTINE sort_list_dbl

!*********************************************************************

      SUBROUTINE sort_list_integer(a,ncol,a_gt_b)
!     *** QUICKSORT
! .. Parameters ..
        INTEGER, PARAMETER :: chgsrt = 10, stckmx = 50
! ..
! .. Scalar Arguments ..
        INTEGER :: ncol
! ..
! .. Array Arguments ..
        INTEGER :: a(ncol)
! ..
! .. Local Scalars ..
        INTEGER :: i, icol, j, l, part, r, stckct, swap1, swap2, temp
        LOGICAL :: q, qdone
! ..
! .. Local Arrays ..
        INTEGER :: stck(2,stckmx)
! ..
! .. Intrinsic Functions ..
        INTRINSIC MAX, MIN
! ..
! .. Non-Generic Interface Blocks ..
        INTERFACE
          FUNCTION a_gt_b(a,b)
! .. Function Return Value ..
            LOGICAL :: a_gt_b
! ..
! .. Scalar Arguments ..
            INTEGER :: a, b
! ..
          END FUNCTION a_gt_b
        END INTERFACE
! ..
        OPTIONAL :: a_gt_b

        stckct = 0
        l = 1
        r = ncol
        qdone = (r-l) < chgsrt

!     *** 'MEDIAN-OF-THREE' MODIFICATION.  SEE 'OVERVIEW' ABOVE

        DO WHILE ( .NOT. (qdone))
          swap1 = (l+r)/2
          swap2 = l + 1

          CALL swap_elements

          IF (qgr(a(l+1),a(r))) THEN
            swap1 = l + 1
            swap2 = r

            CALL swap_elements
          END IF

          IF (qgr(a(l),a(r))) THEN
            swap1 = l
            swap2 = r

            CALL swap_elements
          END IF

          IF (qgr(a(l+1),a(l))) THEN
            swap1 = l + 1
            swap2 = l

            CALL swap_elements
          END IF


!     *** PARTITIONING OF SUBFILE

          i = l + 1
          j = r
          part = l
          DO

            DO
              i = i + 1
              IF ( .NOT. (qgr(a(part),a(i)))) EXIT
            END DO

            DO
              j = j - 1
              IF ( .NOT. (qgr(a(j),a(part)))) EXIT
            END DO

            IF (j<i) THEN
              EXIT
            ELSE
              swap1 = i
              swap2 = j

              CALL swap_elements
            END IF

          END DO

          swap1 = l
          swap2 = j

          CALL swap_elements

!     *** RECURSION STEP

!     *** DEFER FURTHER SORTING OF BOTH SUBFILES TO INSERTION SORT

          IF (MAX((j-l),(r-i+1))<=chgsrt) THEN
            IF (stckct==0) THEN
              qdone = .TRUE.
!     *** POP STACK AND CONTINUE QUICKSORT
            ELSE
              l = stck(1,stckct)
              r = stck(2,stckct)
              stckct = stckct - 1
            END IF

!     *** CONTINUE QUICKSORT ON AT LEAST ONE SUBFILE

          ELSE
!     *** DEFER SMALL TO INSERTION SORT, CONTINUE QUICKSORT ON LARGE

            IF (MIN((j-l),(r-i+1))<=chgsrt) THEN
!     *** LEFT SUBFILE IS LARGE ONE
!     L = L
              IF ((j-l)>=(r-i+1)) THEN
                r = j - 1
!     *** RIGHT SUBFILE IS LARGE ONE
              ELSE
                l = i
!     R = R
              END IF

!     *** CONTINUE QUICKSORT ON BOTH SUBFILES

            ELSE
!     *** STACK IS FULL
              IF (stckct>=stckmx) THEN
                STOP ' STACK OVERFLOW IN SORT_LIST_INTEGER'
              END IF

!     *** PUSH LARGE SUBFILE ONTO STACK, CONTINUE QUICKSORT WITH SMALL
              stckct = stckct + 1
!     *** LEFT SUBFILE IS LARGE ONE

              IF ((j-l)>=(r-i+1)) THEN
                stck(1,stckct) = l
                stck(2,stckct) = j - 1
                l = i
!     R = R

!     *** RIGHT SUBFILE IS LARGE ONE

              ELSE
                stck(1,stckct) = i
                stck(2,stckct) = r
!     L = L
                r = j - 1

              END IF

            END IF

          END IF


        END DO


!     *** INSERTION SORT

        DO i = (ncol-1), 1, -1
!     *** FOR EACH POSITION PRECEDING THE LAST POSITION ...

!     ***
          IF (qgr(a(i),a(i+1))) THEN
            j = i + 1

!     *** A(I) Needs to be moved further out

            DO
              j = j + 1

              IF (j>ncol) THEN
                q = .FALSE.
              ELSE
                q = qgr(a(i),a(j))
              END IF

              IF ( .NOT. (q)) EXIT

            END DO

            temp = a(i)

            DO icol = i, j - 2
              a(icol) = a(icol+1)
            END DO

            a(j-1) = temp

          END IF

        END DO

        RETURN
      CONTAINS

!.....................................................................

        FUNCTION qgr(a,b)
! .. Function Return Value ..
          LOGICAL :: qgr
! ..
! .. Scalar Arguments ..
          INTEGER, INTENT (IN) :: a, b
! ..
! .. Intrinsic Functions ..
          INTRINSIC PRESENT
! ..
          IF (PRESENT(a_gt_b)) THEN
            qgr = a_gt_b(a,b)
          ELSE
            qgr = a >= b
          END IF

          RETURN

        END FUNCTION qgr

!.....................................................................

        SUBROUTINE swap_elements
! .. Local Scalars ..
          INTEGER :: temp
! ..
          temp = a(swap1)
          a(swap1) = a(swap2)
          a(swap2) = temp

        END SUBROUTINE swap_elements

!.....................................................................

      END SUBROUTINE sort_list_integer

!*********************************************************************

      SUBROUTINE sort_list_single(a,ncol,a_gt_b)
!     *** QUICKSORT
! .. Parameters ..
        INTEGER, PARAMETER :: chgsrt = 10, stckmx = 50
! ..
! .. Scalar Arguments ..
        INTEGER :: ncol
! ..
! .. Array Arguments ..
        REAL (spkind) :: a(ncol)
! ..
! .. Local Scalars ..
        REAL (spkind) :: temp
        INTEGER :: i, icol, j, l, part, r, stckct, swap1, swap2
        LOGICAL :: q, qdone
! ..
! .. Local Arrays ..
        INTEGER :: stck(2,stckmx)
! ..
! .. Intrinsic Functions ..
        INTRINSIC MAX, MIN
! ..
! .. Non-Generic Interface Blocks ..
        INTERFACE
          FUNCTION a_gt_b(a,b)
! .. Use Statements ..
            USE biomath_constants_mod
! ..
! .. Function Return Value ..
            LOGICAL :: a_gt_b
! ..
! .. Scalar Arguments ..
            REAL (spkind) :: a, b
! ..
          END FUNCTION a_gt_b
        END INTERFACE
! ..
        OPTIONAL :: a_gt_b

        stckct = 0
        l = 1
        r = ncol
        qdone = (r-l) < chgsrt

!     *** 'MEDIAN-OF-THREE' MODIFICATION.  SEE 'OVERVIEW' ABOVE

        DO WHILE ( .NOT. (qdone))
          swap1 = (l+r)/2
          swap2 = l + 1

          CALL swap_elements

          IF (qgr(a(l+1),a(r))) THEN
            swap1 = l + 1
            swap2 = r

            CALL swap_elements
          END IF

          IF (qgr(a(l),a(r))) THEN
            swap1 = l
            swap2 = r
            CALL swap_elements
          END IF

          IF (qgr(a(l+1),a(l))) THEN
            swap1 = l + 1
            swap2 = l
            CALL swap_elements
          END IF

!     *** PARTITIONING OF SUBFILE

          i = l + 1
          j = r
          part = l
          DO

            DO
              i = i + 1
              IF ( .NOT. (qgr(a(part),a(i)))) EXIT
            END DO

            DO
              j = j - 1
              IF ( .NOT. (qgr(a(j),a(part)))) EXIT
            END DO

            IF (j<i) THEN
              EXIT
            ELSE
              swap1 = i
              swap2 = j
              CALL swap_elements
            END IF

          END DO

          swap1 = l
          swap2 = j

          CALL swap_elements

!     *** RECURSION STEP

!     *** DEFER FURTHER SORTING OF BOTH SUBFILES TO INSERTION SORT

          IF (MAX((j-l),(r-i+1))<=chgsrt) THEN
            IF (stckct==0) THEN
              qdone = .TRUE.
!     *** POP STACK AND CONTINUE QUICKSORT
            ELSE
              l = stck(1,stckct)
              r = stck(2,stckct)
              stckct = stckct - 1
            END IF

!     *** CONTINUE QUICKSORT ON AT LEAST ONE SUBFILE

          ELSE
!     *** DEFER SMALL TO INSERTION SORT, CONTINUE QUICKSORT ON LARGE

            IF (MIN((j-l),(r-i+1))<=chgsrt) THEN
!     *** LEFT SUBFILE IS LARGE ONE
!     L = L
              IF ((j-l)>=(r-i+1)) THEN
                r = j - 1
!     *** RIGHT SUBFILE IS LARGE ONE
              ELSE
                l = i
!     R = R
              END IF

!     *** CONTINUE QUICKSORT ON BOTH SUBFILES

            ELSE
!     *** STACK IS FULL
              IF (stckct>=stckmx) THEN
                STOP ' STACK OVERFLOW IN SORT_LIST_SINGLE'
              END IF

!     *** PUSH LARGE SUBFILE ONTO STACK, CONTINUE QUICKSORT WITH SMALL
              stckct = stckct + 1
!     *** LEFT SUBFILE IS LARGE ONE

              IF ((j-l)>=(r-i+1)) THEN
                stck(1,stckct) = l
                stck(2,stckct) = j - 1
                l = i
!     R = R

!     *** RIGHT SUBFILE IS LARGE ONE

              ELSE
                stck(1,stckct) = i
                stck(2,stckct) = r
!     L = L
                r = j - 1

              END IF

            END IF

          END IF

        END DO

!     *** INSERTION SORT

        DO i = (ncol-1), 1, -1
!     *** FOR EACH POSITION PRECEDING THE LAST POSITION ...

!     ***
          IF (qgr(a(i),a(i+1))) THEN
            j = i + 1

!     *** A(I) Needs to be moved further out

            DO
              j = j + 1

              IF (j>ncol) THEN
                q = .FALSE.

              ELSE
                q = qgr(a(i),a(j))

              END IF

              IF ( .NOT. (q)) EXIT

            END DO

            temp = a(i)

            DO icol = i, j - 2
              a(icol) = a(icol+1)
            END DO

            a(j-1) = temp

          END IF

        END DO

        RETURN

      CONTAINS

!.....................................................................

        FUNCTION qgr(a,b)
! .. Function Return Value ..
          LOGICAL :: qgr
! ..
! .. Scalar Arguments ..
          REAL (spkind), INTENT (IN) :: a, b
! ..
! .. Intrinsic Functions ..
          INTRINSIC PRESENT
! ..
          IF (PRESENT(a_gt_b)) THEN
            qgr = a_gt_b(a,b)
          ELSE
            qgr = a >= b
          END IF

          RETURN

        END FUNCTION qgr

!.....................................................................

        SUBROUTINE swap_elements
! .. Local Scalars ..
          REAL (spkind) :: temp
! ..
          temp = a(swap1)
          a(swap1) = a(swap2)
          a(swap2) = temp

        END SUBROUTINE swap_elements

!.....................................................................

      END SUBROUTINE sort_list_single

!*********************************************************************

    END MODULE biomath_sort_mod
