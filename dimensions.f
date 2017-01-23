      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 200, ny = 100, nz = 12)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 300000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
