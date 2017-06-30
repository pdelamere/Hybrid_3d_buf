      MODULE inter_proc_comm

      USE global
      USE misc
      USE mpi
      USE boundary
      USE grid_interp

      contains

c----------------------------------------------------------------------
      SUBROUTINE exchange_part(var, Ni_out, out_part, in_part, dest, 
     x     source, dim)
c exchanges ions between processors on parallel spatial domain decomposition
c----------------------------------------------------------------------

      real var(Ni_max,dim)

      integer Ni_out, Ni_in
      integer source, dest
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)
      integer stat(MPI_STATUS_SIZE)

      real, dimension(:,:), allocatable :: in_part
      real, dimension(:,:), allocatable :: out_part

      allocate(out_part(Ni_out,dim))

      call MPI_SEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
     x     cartcomm, ierr)
      call MPI_RECV(Ni_in, 1, MPI_INTEGER, source, tag,
     x     cartcomm, stat, ierr)
      
      allocate(in_part(Ni_in,3))

      do m = 1,dim
         out_part(1:Ni_out,m) = 
     x        pack(var(1:Ni_tot,m), .not.in_bounds(1:Ni_tot))
c         xp(1:Ni_tot_in,m) = pack(xp(1:Ni_tot,m), in_bounds(1:Ni_tot))
      enddo

      call pack_pd(xp, in_bounds, dim)

      call MPI_ISEND(out_part, dim*Ni_out, realtype, dest, tag, 
     x     cartcomm, reqs(1), ierr)
      call MPI_IRECV(in_part, dim*Ni_in, realtype, source, tag,
     x     cartcomm, reqs(2), ierr)
      
      call MPI_WAITALL(2, reqs, stats, ierr)

                 
      var(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


      deallocate(out_part)
      deallocate(in_part)

      end SUBROUTINE exchange_part


      end MODULE inter_proc_comm
