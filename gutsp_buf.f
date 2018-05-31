      MODULE gutsp_buf

      USE dimensions
      USE global
      USE misc
      USE part_init

      contains


c----------------------------------------------------------------------
      SUBROUTINE part_setup_buf(xp_buf,vp_buf)
c----------------------------------------------------------------------

      real vp_buf(Ni_max_buf,3)
      real xp_buf(Ni_max_buf,3)
      real rnd,f,v
      real vx,vy,vz
      real vol_buf

      integer flg
      integer Ni_tot_buf_1

      vol_buf = (qy(ny-1)-qy(1))*(qz(nz-1)-qz(1))*dx_buf

      Ni_tot_buf = nint(nf_init*vol_buf*beta)

      

c initialize protons

      do 10 l = 1,Ni_tot_buf

         xp_buf(l,1) = qx(nx)+(1.0-pad_ranf())*dx_buf
         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         
         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         call maxwl_init(vth,vx,vy,vz)

         vp_buf(l,1) = -vsw + vx 
         vp_buf(l,2) = vy 
         vp_buf(l,3) = vz 

 10   continue
         
      
      mrat_buf(1:Ni_tot_buf) = 1.0
      mrat_buf(Ni_tot_buf+1:) = 1.0/m_pu       
      beta_p_buf(1:Ni_tot_buf) = 1.0
      beta_p_buf(Ni_tot_buf+1:) = 1.0
      tags_buf(1:Ni_tot_buf) = 1

c initialize He++ (m/q =2) 

      Ni_tot_buf_1 = Ni_tot_buf
      Ni_tot_buf = Ni_tot_buf_1 + f_mq_2*Ni_tot_buf_1

      do 20 l = Ni_tot_buf_1+1,Ni_tot_buf

         xp_buf(l,1) = qx(nx)+(1.0-pad_ranf())*dx_buf
         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         
         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         call maxwl_init(vth,vx,vy,vz)


         vp_buf(l,1) = -vsw + vx 
         vp_buf(l,2) = vy 
         vp_buf(l,3) = vz 

         mrat_buf(l) = 1./2.
         beta_p_buf(l) = b_mq_2
         tags_buf(l) = 1

 20   continue
         
c add shell distribution

c      Ni_tot_buf_1 = Ni_tot_buf 
c      
c      Ni_tot_buf = Ni_tot_buf_1 + f_shl*Ni_tot_buf_1
c      
c      do 69 l = Ni_tot_buf_1+1,Ni_tot_buf
c         
c         xp_buf(l,1) = qx(nx)+(1.0-pad_ranf())*dx_buf
c         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
c         xp_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
c         
c         mrat_buf(l) = 1.0
c         beta_p_buf(l) = b_shl
c         tags_buf(l) = 1
c         
c         vz = pad_ranf()*2 - 1
c         tmp = sqrt(1-vz**2)
c         phi = 2*PI*pad_ranf()
c         vx = tmp*cos(phi)
c         vy = tmp*sin(phi)
c         
c         vp_buf(l,1) = -vsw+vsw*vx !+dvx
c         vp_buf(l,2) = vsw*vy !+dvz 
c         vp_buf(l,3) = vsw*vz
c         
c 69   enddo


      return
      end SUBROUTINE part_setup_buf
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE get_Ep_buf(Ep_buf,b0,xp_buf,up)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real Ep_buf(Ni_max_buf,3)
      real b0(nx,ny,nz,3)
      real xp_buf(Ni_max_buf,3)
      real up(nx,ny,nz,3)

      real aa(3), bb(3), cc(3)    !dummy vars for doing cross product
      real ntot                   !total plasma density
      real fnf,fnp                !fraction nf,np of total n
      real up3(3),  !vector field values at Ba position
     x      btc3(3)

c      real eoverm
c      parameter (eoverm = q/mO)

      do 10 l=1,Ni_tot_buf

cc         i = floor(xp_buf(l,1)/dx)
c         j = floor(xp_buf(l,2)/dy)
c         k = floor(xp_buf(l,3)/delz)
         
c         i=0
c 31      continue
c         i = i + 1
c         if (xp(l,1) .gt. qx(i)) go to 31 !find i on non-uniform 
c         i = i-1
c         ijkp(l,1)= i

         j=0
 33      continue
         j = j + 1
         if (xp_buf(l,2) .gt. qy(j)) go to 33 !find i on non-uniform 
         j = j-1
         ijkp(l,2)= j
        
         k=0
 31      continue
         k = k + 1
         if (xp_buf(l,3) .gt. qz(k)) go to 31 !find k on non-uniform 
         k = k-1
         ijkp(l,3)= k
         

c         up3(1) = -vsw
c         up3(2) = 0.0
c         up3(3) = 0.0

         up3(1) = up(nx,j,k,1)
         up3(2) = up(nx,j,k,2)
         up3(3) = up(nx,j,k,3)

         btc3(1) = b0(1,1,1,1)!0.0
         btc3(2) = b0(1,1,1,2)!b0_init*eoverm
         btc3(3) = b0(1,1,1,3)!0.0


         do 20 m=1,3
            aa(m) = - up3(m)
            bb(m) = btc3(m)                   
 20         continue

         cc(1) = aa(2)*bb(3) - aa(3)*bb(2)    !do cross product
         cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
         cc(3) = aa(1)*bb(2) - aa(2)*bb(1)

         do 30 m=1,3
            Ep_buf(l,m) = cc(m) 
            Ep_buf(l,m) = Ep_buf(l,m)*mrat_buf(l) !O_to_Ba
 30         continue


 10      continue

      return
      end SUBROUTINE get_Ep_buf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vplus_vminus_buf(Ep_buf,vp_buf,vplus_buf,
     x                                vminus_buf,b0)
c----------------------------------------------------------------------

c      include 'incurv.h'

      real Ep_buf(Ni_max_buf,3),
     x     vp_buf(Ni_max_buf,3),      !particle velocities at t level n-1/2
     x     vplus_buf(Ni_max_buf,3),
     x     vminus_buf(Ni_max_buf,3),
     x     b0(nx,ny,nz,3)

      real a1, a2, a3,       !coefficients for calculating vplus
     x     a_d,              !denominator for a1,a2,a3
     x     B2,dt2,           !B*B,dt*dt
     x     Bx,By,Bz,         !B for cross product call
     x     vminus_buf_x_B(3),    !v- x B
     x     vminus_buf_dot_B     !v- . B

      real btc3(3)
      
      do 10 m=1,3
         do 10 l=1,Ni_tot_buf 
            vminus_buf(l,m) = vp_buf(l,m) + 0.5*dt*Ep_buf(l,m)
 10         continue

   
      do 30 l=1,Ni_tot_buf

         btc3(1) = b0(1,1,1,1)!0.0
         btc3(2) = b0(1,1,1,2)!b0_init*eoverm
         btc3(3) = b0(1,1,1,3)!0.0

         vminus_buf_x_B(1) = vminus_buf(l,2)*btc3(3)*mrat_buf(l) - !O_to_Ba - 
     x                     vminus_buf(l,3)*btc3(2)*mrat_buf(l)   !O_to_Ba
         vminus_buf_x_B(2) = vminus_buf(l,3)*btc3(1)*mrat_buf(l) - !O_to_Ba - 
     x                     vminus_buf(l,1)*btc3(3)*mrat_buf(l)   !O_to_Ba
         vminus_buf_x_B(3) = vminus_buf(l,1)*btc3(2)*mrat_buf(l) - !O_to_Ba -
     x                     vminus_buf(l,2)*btc3(1)*mrat_buf(l)   !O_to_Ba

         vminus_buf_dot_B = vminus_buf(l,1)*btc3(1)*mrat_buf(l) + !O_to_Ba +
     x                     vminus_buf(l,2)*btc3(2)*mrat_buf(l) + !O_to_Ba +
     x                     vminus_buf(l,3)*btc3(3)*mrat_buf(l)   !O_to_Ba

         Bx = btc3(1)*mrat_buf(l) !O_to_Ba
         By = btc3(2)*mrat_buf(l) !O_to_Ba
         Bz = btc3(3)*mrat_buf(l) !O_to_Ba
      
         B2 = Bx**2 + By**2 + Bz**2
         dt2 = dt**2

         a_d = 1 + (B2*dt2/4.0)
         a1 = (1 - (B2*dt2/4.0)) / a_d
         a2 = dt / a_d
         a3 = 0.5*dt2 / a_d

         do 40 m=1,3
            vplus_buf(l,m) = a1*vminus_buf(l,m) + a2*vminus_buf_x_B(m) + 
     x           a3*vminus_buf_dot_B*btc3(m)*mrat_buf(l) !O_to_Ba
 40      continue

 30   continue

      return
      end SUBROUTINE get_vplus_vminus_buf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_vp_buf_final(Ep_buf,vp_buf,vplus_buf)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real Ep_buf(Ni_max_buf,3),
     x     vp_buf(Ni_max_buf,3),    !particle velocities at t level n+1/2
     x     vplus_buf(Ni_max_buf,3)

      do 10 m=1,3
         do 10 l = 1,Ni_tot_buf
            vp_buf(l,m) = vplus_buf(l,m) + 0.5*dt*Ep_buf(l,m)
 10         continue

      return
      end SUBROUTINE get_vp_buf_final
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE exchange_ion_in_buf(xp_buf,vp_buf,xp,vp,vp1)
c----------------------------------------------------------------------
c Exchange ions from the inflow buffer into the main domain.

      real xp_buf(Ni_max_buf,3),
     x     vp_buf(Ni_max_buf,3),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)
      integer Ni_tot_in, Ni_out
     
      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
      real, dimension(:), allocatable :: out_mrat
      real, dimension(:), allocatable :: out_beta_p
      real, dimension(:), allocatable :: out_tags

      in_bounds_buf(1:Ni_tot_buf) = .true.
      in_bounds_buf(Ni_tot_buf+1:) = .false.

      where(xp_buf(1:Ni_tot_buf,1) .lt. qx(nx))
     X     in_bounds_buf(1:Ni_tot_buf) = .false.

      Ni_tot_in = count(in_bounds_buf)
      Ni_out = count(.not.in_bounds_buf(1:Ni_tot_buf))

      allocate(out_xp(Ni_out,3))
      allocate(out_vp(Ni_out,3))
      allocate(out_mrat(Ni_out))
      allocate(out_beta_p(Ni_out))
      allocate(out_tags(Ni_out))
      
      do m = 1,3 

        out_xp(1:Ni_out,m) = pack(xp_buf(1:Ni_tot_buf,m), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))
        out_vp(1:Ni_out,m) = pack(vp_buf(1:Ni_tot_buf,m), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))
        xp_buf(1:Ni_tot_in,m) = pack(xp_buf(1:Ni_tot_buf,m), 
     x                               in_bounds_buf(1:Ni_tot_buf))
        vp_buf(1:Ni_tot_in,m) = pack(vp_buf(1:Ni_tot_buf,m), 
     x                               in_bounds_buf(1:Ni_tot_buf))

      enddo

      out_mrat(1:Ni_out) = pack(mrat_buf(1:Ni_tot_buf), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))

      out_beta_p(1:Ni_out) = pack(beta_p_buf(1:Ni_tot_buf), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))
      out_tags(1:Ni_out) = pack(tags_buf(1:Ni_tot_buf), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))

      mrat_buf(1:Ni_tot_in) = pack(mrat_buf(1:Ni_tot_buf), 
     x                             in_bounds_buf(1:Ni_tot_buf))
      
      beta_p_buf(1:Ni_tot_in) = pack(beta_p_buf(1:Ni_tot_buf), 
     x                             in_bounds_buf(1:Ni_tot_buf))
      tags_buf(1:Ni_tot_in) = pack(tags_buf(1:Ni_tot_buf), 
     x                             in_bounds_buf(1:Ni_tot_buf))
      



      do m = 1,3
         xp(Ni_tot+1:Ni_tot+Ni_out,m) = out_xp(:,m)
         vp(Ni_tot+1:Ni_tot+Ni_out,m) = out_vp(:,m)
         vp1(Ni_tot+1:Ni_tot+Ni_out,m) = out_vp(:,m)
      enddo


      do l = Ni_tot+1,Ni_tot+Ni_out 

         ii=0
 51      continue
         ii = ii + 1
         if (xp(l,1) .gt. qx(ii)) go to 51 !find i on non-uniform 
         ii = ii-1
         ijkp(l,1)= ii

         jj=0
 53      continue
         jj = jj + 1
         if (xp(l,2) .gt. qy(jj)) go to 53 !find i on non-uniform 
         jj = jj-1
         ijkp(l,2)= jj

         kk=0
 50      continue
         kk = kk + 1
         if (xp(l,3) .gt. qz(kk)) go to 50 !find k on non-uniform 
         kk = kk-1
         ijkp(l,3)= kk

      enddo

      mrat(Ni_tot+1:Ni_tot+Ni_out) = out_mrat(:)
      beta_p(Ni_tot+1:Ni_tot+Ni_out) = out_beta_p(:)
      tags(Ni_tot+1:Ni_tot+Ni_out) = out_tags(:)
      
      do l = Ni_tot+1,Ni_tot+Ni_out 
         do m=1,3
            input_E = input_E + 
     x         0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /(beta*beta_p(l))
         enddo
      enddo

      Ni_tot = Ni_tot + Ni_out

      Ni_tot_buf = count(in_bounds_buf)

      deallocate(out_xp)
      deallocate(out_vp)
      deallocate(out_mrat)
      deallocate(out_beta_p)
      deallocate(out_tags)


      return
      end SUBROUTINE exchange_ion_in_buf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
c----------------------------------------------------------------------

      real xp_buf(Ni_max_buf,3),
     x     vp_buf(Ni_max_buf,3),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)
      integer zindx(1)
      integer Ni_tot_in, Ni_out
     
      real dth                 !half time step

      real v2
      real vx,vy,vz
      
      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
      real, dimension(:), allocatable :: out_mrat

      dth = dt/2


      do 10 l=1,Ni_tot_buf                   !make 1/2 time step advance

         xp_buf(l,1) = xp_buf(l,1) + dth*vp_buf(l,1)
         xp_buf(l,2) = xp_buf(l,2) + dth*vp_buf(l,2)
         xp_buf(l,3) = xp_buf(l,3) + dth*vp_buf(l,3)


 10      continue

      where (xp_buf(:,2) .ge. qy(ny-1))
         xp_buf(:,2) = qy(1) + ( xp_buf(:,2) - qy(ny-1) )
      endwhere

      where (xp_buf(:,2) .le. qy(1)) 
         xp_buf(:,2) = qy(ny-1) - (qy(1) - xp_buf(:,2))
      endwhere


      where (xp_buf(:,3) .ge. qz(nz))
         xp_buf(:,3) = qz(2) + ( xp_buf(:,3) - qz(nz) )
      endwhere

      where (xp_buf(:,3) .le. qz(2)) 
         xp_buf(:,3) = qz(nz) - (qz(2) - xp_buf(:,3))
      endwhere

      return
      end SUBROUTINE move_ion_half_buf
c----------------------------------------------------------------------
      end MODULE gutsp_buf

