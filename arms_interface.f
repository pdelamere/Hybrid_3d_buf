      module arms_interface
      use iso_c_binding, only: C_DOUBLE, c_ptr
      implicit none
      interface
          subroutine arms_vs_setup(ninit, b, state)
     x               bind(C,name="arms_vs_setup_")
          use iso_c_binding, only: C_INT, C_DOUBLE, c_ptr
          integer(kind=C_INT), intent(in) :: ninit
          real(kind=C_DOUBLE), intent(in) :: b
          type(c_ptr), intent(out) :: state
          end subroutine
          subroutine arms_sample(xsamp, nsamp, state)
     x               bind(C,name="arms_sample_")
          use iso_c_binding, only: C_INT, C_DOUBLE, c_ptr
          integer(kind=C_INT), intent(in)  :: nsamp
          real(kind=C_DOUBLE), intent(out) :: xsamp(nsamp)
          type(c_ptr), intent(in) :: state
          end subroutine
      end interface
      contains
      subroutine arms_sample_sp(xsamp, nsamp, state)
      ! A single precision (sp) version of arms_sample
      integer, intent(in) :: nsamp
      real, intent(out) :: xsamp(nsamp)
      type(c_ptr), intent(in) :: state

      real(kind=C_DOUBLE) :: dxsamp(nsamp)

      call arms_sample(dxsamp, nsamp, state)
      xsamp = real(dxsamp)
      end subroutine
      end module
