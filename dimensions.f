      MODULE DIMENSIONS

c simulation domain dimensions
      integer nx, ny, nz
      integer num_cells
      integer num_buf_cells

      PARAMETER (nx = 25, ny = 25, nz = 6)
      PARAMETER (num_cells = nx*ny*nz)
      PARAMETER (num_buf_cells = ny*nz)
c particle array dimensions

      integer Ni_max
      PARAMETER (Ni_max = 20000000)

      END MODULE DIMENSIONS
