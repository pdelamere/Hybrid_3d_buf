      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 181, ny = 79, nz = 21)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 8000000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
