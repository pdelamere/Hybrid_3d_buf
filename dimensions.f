      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 121/2, ny = 69/2, nz = 21/2)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 6000000)
      PARAMETER (Ni_max_buf = Ni_max/6)

      END MODULE DIMENSIONS
