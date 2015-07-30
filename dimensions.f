      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 119, ny = 59, nz = 13)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 4000000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
