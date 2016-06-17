      MODULE DIMENSIONS

c simulation domain dimensions

<<<<<<< HEAD
      PARAMETER (nx = 109, ny = 109, nz = 11)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 7000000)
=======
      PARAMETER (nx = 150, ny = 100, nz = 10)
c particle array dimensions

      integer*4 Ni_max, Ni_max_buf
      PARAMETER (Ni_max = 12000000)
>>>>>>> origin/master
      PARAMETER (Ni_max_buf = Ni_max/10)

      END MODULE DIMENSIONS
