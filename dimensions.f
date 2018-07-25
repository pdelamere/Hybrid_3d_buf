      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 150, ny = 150, nz = 3)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 25000000)
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
