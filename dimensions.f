      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 75, ny = 75, nz = 10)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 1500000)
      PARAMETER (Ni_max_buf = Ni_max/3)

      END MODULE DIMENSIONS
