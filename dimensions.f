      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 220, ny = 159, nz = 13)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 6000000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
