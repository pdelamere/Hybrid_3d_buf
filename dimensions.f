      MODULE DIMENSIONS

c simulation domain dimensions
      integer nx, ny, nz
      integer num_cells
      integer num_buf_cells

      PARAMETER (nx = 50, ny = 50, nz = 20)
      PARAMETER (num_cells = nx*ny*nz)
      PARAMETER (num_buf_cells = ny*nz)
c particle array dimensions

      integer Ni_max
      PARAMETER (Ni_max = 2000000)

      END MODULE DIMENSIONS
