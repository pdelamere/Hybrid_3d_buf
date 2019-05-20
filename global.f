      MODULE GLOBAL


      USE inputs
      USE dimensions
      save

c contains the common data structures to be used by the various
c simulation subroutines 

c raw grid coordinate data
      real gz(nz)
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction 
      real dx_grid(nx),dx_cell(nx)
      real dy_grid(ny),dy_cell(ny)
      real xrat(nx), yrat(ny), zrat(nz)
      
c Total number of ions produced at a given time during the simulation
      integer Ni_tot, Ni_tot_0
      integer Ni_tot_buf
      integer Ni_tot_buf_0
      integer Ni_thermal_H, Ni_thermal_He, Ni_shell_H, Ni_shell_He
      integer Ni_thermal_H_buf, Ni_thermal_He_buf
      integer Ni_shell_H_buf, Ni_shell_He_buf
      integer Ni_thermal_H_buf_0, Ni_thermal_He_buf_0
      integer Ni_shell_H_buf_0, Ni_shell_He_buf_0

      ! Ni_s(1) === number of thermal H macroparticles
      ! Ni_s(2) === number of thermal He macroparticles
      ! Ni_s(3) === number of shell H macroparticles
      ! Ni_s(4) === number of shell He macroparticles
      integer, dimension(4) :: Ni_s

c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical in_bounds(Ni_max)

c Random seed
      integer*4 seed

c Indices for boundary cells
      integer ip,im,jp,jm,kp,km

c Total input energy (ions) 
      real input_E, bndry_Eflux
      real input_chex, input_bill

c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)

c Dummy variables for doing vector operations and particle flux
c calculation....etc.
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)

c Variable time step for fluid velocity update used in both
c particle update and fluid update

      real dtf,ntf,dtsub

c Weight variables for trilinear interpolation
      real wght(Ni_max,8)!, wquad(Ni_max,3)

c variable for anomlous resistivity
      real nuei

c variable for particle scale
      real beta, beta_p(Ni_max)

c mass array for multi-ion species
      real mrat(Ni_max)

c mixing array
      real mixed(nx,ny,nz)

c parallel processor info
      integer procnum, my_rank, up_proc, down_proc, cart_rank, cartcomm

c Tags for particle origin
      real tags(Ni_max)

      real vsw_us(ny,nz,3)
      real b0_us(ny,nz,3)

c Buffer variables
      real, dimension(:,:), allocatable :: xp_buf
      real, dimension(:,:), allocatable :: vp_buf
      real, dimension(:), allocatable :: mrat_buf
      real, dimension(:), allocatable :: beta_p_buf
      real, dimension(:), allocatable :: tags_buf

      real, dimension(:,:), allocatable :: Ep_buf
      real, dimension(:,:), allocatable :: vplus_buf
      real, dimension(:,:), allocatable :: vminus_buf

      real, dimension(:,:), allocatable :: xp_out_buf
      real, dimension(:,:), allocatable :: vp_out_buf
      real, dimension(:,:), allocatable :: E_out_buf
      real, dimension(:,:), allocatable :: B_out_buf
      real, dimension(:), allocatable :: mrat_out_buf
      logical, dimension(:), allocatable :: in_bounds_buf

c Computed constants for the buffer
      real dx_buf

      integer Ni_max_buf

c beta_p values for the various subpopulations
      real b_sw_thermal_H, b_sw_thermal_He, b_sw_shell_H, b_sw_shell_He

      end MODULE GLOBAL
