    MODULE cdf_f_mod
! ----------------------------------------------------------------------

!                                cdf_f_mod
!                                *=*=*=*=*

!  -  SUBROUTINE CDF_F( WHICH, CUM, CCUM, F, DFN, DFD, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CUM_F(F,  DFN, DFD, STATUS, CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION CCUM_F(F, DFN, DFD, STATUS,  CHECK_INPUT)
!  -  REAL (dpkind) FUNCTION INV_F( CUM, CCUM, DFN, DFD, STATUS, CHECK_INPUT)

!                             The Distribution
!                             ================

!  F is the distribution of the ratio of two independent random
!  variables. The numerator random variable is distributed as chi-squared
!  with DF degrees of freedom divided by DF. The denominator random
!  variable is distributed as chi-squared with DFD degrees of freedom
!  divided by DFD.The density of the f distribution is defined on x in
!  [0,+Infinity] and is proportional to:

!                                (DFN-2)/2)        
!                               x                  
!                        ------------------------- 
!                                      (DFN+DFD)/2 
!                        [1+(DFN/DFD)x]            

!                                Arguments
!                                =========

!  - INTEGER, INTENT(IN) :: WHICH. Integer indicating which of the next
!    four arguments is to be calculated.
!  Input Range: [ 1:2 ]
!     1.  CUM and CCUM 
!     2.  F 
!  NOTE: DFN and DFD will not be computed because CUM is not monotone in
!    either argument.
!  - REAL (dpkind), OPTIONAL :: CUM. The CDF of the f distribution.
!  Range: [ 0:1-10^-10 ]
!  - REAL (dpkind), OPTIONAL :: CCUM. One minus the CDF of the f
!    distribution.
!  Range: [ 10^-10:1 ]
!  - REAL (dpkind) :: F. The upper limit of integration of the f density.
!    The lower limit is 0.
!  Input Range: [ 0:10^100 ]
!  - REAL (dpkind) :: DFN. The numerator degrees of freedom.
!  Range: [ 10^-3:10^10 ]
!  - REAL (dpkind) :: DFD. The denominator degrees of freedom.
!  Range: [ 10^-3:10^10 ]
!  - INTEGER, OPTIONAL, INTENT(OUT) :: STATUS. Return code. 
!     -1 WHICH outside input range 
!     -2 CUM outside range 
!     -3 CCUM outside range 
!     -4 F outside range 
!     -5 DFN outside range 
!     -6 DFD outside range 
!      3 CUM + CCUM is not nearly one 
!     10 can not solve cdf_beta calling from local_cum_f
!    -50 Answer (if any) is BELOW the LOWER search bound 
!     50 Answer (if any) is ABOVE the UPPER search bound 
!  - LOGICAL, INTENT(IN), OPTIONAL :: CHECK_INPUT. If PRESENT and
!       .TRUE. input argument values are not checked for validity.

!   NOTE: CUM and CCUM must add to (nearly) one.

!   NOTE: The value of the CDF of the f distribution is NOT necessarily
!         monotone in either degree of freedom argument. There may thus
!         be two values that provide a given DCF value. 
!         This routine assumes monotonicity and will find an arbitrary one
!         of the two values.
! ----------------------------------------------------------------------
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: ccum_f, cdf_f, cum_f, inv_f
! ..
! .. Parameters ..
      REAL (dpkind), PARAMETER :: tiny = 1.0E-100_dpkind
! ..
    CONTAINS

!*********************************************************************

      FUNCTION ccum_f(f,dfn,dfd,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: ccum_f
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: dfd, dfn, f
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_f(which=1,ccum=ccum_f,f=f,dfn=dfn,dfd=dfd,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION ccum_f

!*********************************************************************

      SUBROUTINE cdf_f(which,cum,ccum,f,dfn,dfd,status,check_input)
! .. Use Statements ..
        USE cdf_aux_mod
        USE zero_finder
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind) :: dfd, dfn, f
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
        INTEGER :: zf_status
        LOGICAL :: has_status, local_check_input, match_cum
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        has_status = PRESENT(status)

! status = 0 means NO error

        IF (has_status) THEN
          status = 0
        END IF

        CALL check_complements(cum,ccum,the_f%name,'cum','ccum', &
          local_cum,local_ccum,set_values=(which/=1),bad_status=3, &
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

        IF (PRESENT(check_input)) THEN
          local_check_input = check_input
        ELSE
          local_check_input = .TRUE.
        END IF

        IF (local_check_input) THEN
          CALL validate_parameters(the_f,which,params,status)

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

          CALL local_cum_f(f,dfn,dfd,local_cum,local_ccum,status)

          IF (PRESENT(cum)) cum = local_cum
          IF (PRESENT(ccum)) ccum = local_ccum

        CASE (2)
          f = five
          zf_status = 0

          CALL cdf_set_zero_finder(the_f,3,local)

          DO
            CALL rc_interval_zf(zf_status,f,fx,local)

            IF (zf_status/=1) EXIT

            CALL local_cum_f(f,dfn,dfd,try_cum,try_ccum,status)

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

! Can NOT solve for Degrees of Freedom

        END SELECT

        IF (has_status) THEN
! Set the status of the zero finder
          CALL cdf_finalize_status(local,status)
        END IF

        RETURN

      END SUBROUTINE cdf_f

!*********************************************************************

      FUNCTION cum_f(f,dfn,dfd,status,check_input)
! .. Function Return Value ..
        REAL (dpkind) :: cum_f
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (IN) :: dfd, dfn, f
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_f(which=1,cum=cum_f,f=f,dfn=dfn,dfd=dfd,status=status, &
          check_input=check_input)

        RETURN

      END FUNCTION cum_f

!*********************************************************************

      FUNCTION inv_f(cum,ccum,dfn,dfd,status,check_input)
! Solves the INVERSE F-distribution problem:
! Given cum (or/and ccum), dfn, and dfd, computes f.
! .. Function Return Value ..
        REAL (dpkind) :: inv_f
! ..
! .. Scalar Arguments ..
        REAL (dpkind), OPTIONAL :: ccum, cum
        REAL (dpkind), INTENT (IN) :: dfd, dfn
        INTEGER, OPTIONAL, INTENT (OUT) :: status
        LOGICAL, OPTIONAL, INTENT (IN) :: check_input
! ..
        CALL cdf_f(which=2,cum=cum,ccum=ccum,f=inv_f,dfn=dfn,dfd=dfd, &
          status=status,check_input=check_input)

        RETURN

      END FUNCTION inv_f

!*********************************************************************

      SUBROUTINE local_cum_f(f,dfn,dfd,cum,ccum,status)
!----------------------------------------------------------------------

!                              Function

!     Computes the integral from 0 to F of the f-density.
!     f    --> Upper limit of integration of the f-density.
!     dfn  --> Degrees of freedom of the numerator
!     dfd  --> Degrees of freedom of the denominator
!     cum  <-- Cumulative f-distribution.
!     ccum <-- Complement of Cumulative f-distribution.
!     status <-- tells the caller whether the problem was successfully solved.
!                Set to 10 if cdf_beta has no answer and status was used in
!                the initial call to cdf_f (otherwise, it breaks down)

!                              Method

!     Formula   26.6.2   of   Abramowitz   and   Stegun,  Handbook  of
!     Mathematical  Functions (1966) is used to reduce the computation
!     of the  cumulative  distribution function for the  F  variate to
!     that of an incomplete beta.

!                              Note
!     If F is less than or equal to 0, 0 is returned.
!------------------------------------------------------------------
! .. Use Statements ..
        USE cdf_beta_mod
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: ccum, cum
        REAL (dpkind), INTENT (IN) :: dfd, dfn, f
        INTEGER, OPTIONAL, INTENT (OUT) :: status
! ..
! .. Local Scalars ..
        REAL (dpkind) :: dsum, prod, xx, yy
        INTEGER :: beta_status
! ..
! .. Intrinsic Functions ..
        INTRINSIC PRESENT
! ..
        IF (f<=zero) THEN
          cum = zero
          ccum = one
          RETURN
        END IF

        beta_status = 0

        prod = dfn*f

!     XX is such that the incomplete beta with parameters
!     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM

!     YY is 1 - XX

!     Calculate the smaller of XX and YY accurately

        dsum = dfd + prod

        IF (xx>half) THEN
          yy = prod/dsum
          xx = one - yy
        ELSE
          xx = dfd/dsum
          yy = one - xx
        END IF

        CALL cdf_beta(1,ccum,cum,xx,yy,half*dfd,half*dfn, &
          status=beta_status)

        IF (beta_status/=0) THEN
! cdf_beta has NO answer.

          IF (PRESENT(status)) THEN
            status = 10
          ELSE
            WRITE (*,*) 'Error in local_cum_f call to cdf_beta'
            WRITE (*,*) 'Status: ', beta_status

            STOP 'Error in local_cum_f call to cdf_beta'
          END IF

        END IF

        RETURN

      END SUBROUTINE local_cum_f

!*********************************************************************

    END MODULE cdf_f_mod
