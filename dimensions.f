      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 121, ny = 69, nz = 21)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 8000000)
      PARAMETER (Ni_max_buf = Ni_max/5)

      END MODULE DIMENSIONS
