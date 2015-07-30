      program main
      implicit none
c      integer :: i, count, seed
      real, dimension(:,:), allocatable :: r
      real, dimension(:), allocatable :: a
      logical, dimension(:), allocatable :: msk
      integer*4 i,nx,nmax,cnt,b(2)
      
      parameter(nx = 9000000)
      parameter(nmax = 100000)

      allocate(r(nx,3))
      allocate(msk(nx))

c      b = shape(r)
c      write(*,*) 'size r....',b(2)
c      stop

      do i = 1,nx
         r(i,:) = i
      enddo

      msk(:) = .true.
      msk(1:nmax) = .false.

      allocate(a(count(msk)))
      cnt = 1
      do i = 1,nx
         if (msk(i) .eq. .true.) then
            r(cnt,2) = r(i,2)
            cnt = cnt + 1
         endif
      enddo
 

 
c      a(:) = pack(r(:,2), msk(:))

c      r(1:count(msk),2) = a(:) 

    
      write(*,*) r(1,2),r(count(msk),2)
c      write(*,*) 'msk....',pack(r(:,2), msk(:))

      


      deallocate(r)
c      deallocate(a) 
      deallocate(msk)

      end Program main
