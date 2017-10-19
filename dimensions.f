      MODULE DIMENSIONS

c simulation domain dimensions

      PARAMETER (nx = 75, ny = 75, nz = 6)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 4000000 )
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
