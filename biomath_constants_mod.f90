    MODULE biomath_constants_mod
! .. Default Accessibility ..
      PUBLIC
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: dpkind = KIND(0.0D0)
      INTEGER, PARAMETER :: spkind = KIND(0.0)
      INTEGER, PARAMETER :: stdin = 5, stdout = 6
      REAL (dpkind), PARAMETER :: eight = 8.0E0_dpkind
      REAL (dpkind), PARAMETER :: five = 5.0E0_dpkind
      REAL (dpkind), PARAMETER :: four = 4.0E0_dpkind
      REAL (dpkind), PARAMETER :: hundred = 100.0E0_dpkind
      REAL (dpkind), PARAMETER :: nine = 9.0E0_dpkind
      REAL (dpkind), PARAMETER :: one = 1.0E0_dpkind
      REAL (dpkind), PARAMETER :: seven = 7.0E0_dpkind
      REAL (dpkind), PARAMETER :: six = 6.0E0_dpkind
      REAL (dpkind), PARAMETER :: sixth = one/six
      REAL (dpkind), PARAMETER :: ten = 10.0E0_dpkind
      REAL (dpkind), PARAMETER :: tenth = one/ten
      REAL (dpkind), PARAMETER :: thousand = 1000.0E0_dpkind
      REAL (dpkind), PARAMETER :: thousandth = one/thousand
      REAL (dpkind), PARAMETER :: three = 3.0E0_dpkind
      REAL (dpkind), PARAMETER :: twelve = 12.0E0_dpkind
      REAL (dpkind), PARAMETER :: two = 2.0E0_dpkind
      REAL (dpkind), PARAMETER :: zero = 0.0E0_dpkind
      REAL (dpkind), PARAMETER :: eighth = one/eight
      REAL (dpkind), PARAMETER :: fifth = one/five
      REAL (dpkind), PARAMETER :: fourth = one/four
      REAL (dpkind), PARAMETER :: half = one/two
      REAL (dpkind), PARAMETER :: hundredth = one/hundred
      REAL (dpkind), PARAMETER :: third = one/three
! ..
! .. Intrinsic Functions ..
      INTRINSIC KIND
! ..
    END MODULE biomath_constants_mod
