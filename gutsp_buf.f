      MODULE gutsp_buf

      USE dimensions
      USE global
      USE misc
      USE part_init

      contains


c----------------------------------------------------------------------
      SUBROUTINE part_setup_buf(xp_buf,vp_buf)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vp_buf(Ni_max_buf,3)
      real xp_buf(Ni_max_buf,3)
      real rnd,f,v
      real vx,vy,vz
      real vol_buf

      integer flg
      integer Ni_tot_buf_1

      vol_buf = (qy(ny-1)-qy(1))*(qz(nz-1)-qz(1))*dx_buf

      Ni_tot_buf = nint(nf_init*vol_buf*beta)

c      write(*,*) 'Ni_tot_buf....',Ni_tot_buf,Ni_max_buf,Ni_tot,my_rank
      

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
         
      
c      m_arr_buf(1:Ni_tot_buf) = mproton
c      m_arr_buf(Ni_tot_buf+1:) = m_pu*mproton 
      mrat_buf(1:Ni_tot_buf) = 1.0
      mrat_buf(Ni_tot_buf+1:) = 1.0/m_pu       
      beta_p_buf(1:Ni_tot_buf) = 1.0
      beta_p_buf(Ni_tot_buf+1:) = 1.0

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

c         m_arr_buf(l) = 2*mproton
         mrat_buf(l) = 1./2.
         beta_p_buf(l) = b_mq_2

 20   continue
         
c add shell distribution

      Ni_tot_buf_1 = Ni_tot_buf 
      
      Ni_tot_buf = Ni_tot_buf_1 + f_shl*Ni_tot_buf_1
      
      do 69 l = Ni_tot_buf_1+1,Ni_tot_buf
         
         xp_buf(l,1) = qx(nx)+(1.0-pad_ranf())*dx_buf
         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         
c         m_arr_buf(l) = mproton
         mrat_buf(l) = 1.0
         beta_p_buf(l) = b_shl
         
         theta = pad_ranf()*PI
         phi = pad_ranf()*2*PI
         
         vp_buf(l,1) = -vsw+vsw*cos(phi)*sin(theta) !+dvx
         vp_buf(l,2) = vsw*sin(phi)*sin(theta) !+dvz 
         vp_buf(l,3) = vsw*cos(theta)
         
 69   enddo


      return
      end SUBROUTINE part_setup_buf
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE part_setup_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
     x            B_out_buf,mrat_out_buf,b0)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real vp_out_buf(Ni_max_buf,3)
      real xp_out_buf(Ni_max_buf,3)
      real E_out_buf(Ni_max_buf,3)
      real B_out_buf(Ni_max_buf,3)
      real mrat_out_buf(Ni_max_buf)
c      real m_arr_out_buf(Ni_max_buf)
      real b0(nx,ny,nz,3)
      real rnd,f,v
      real vx,vy,vz
      real vol_buf

      integer flg
      
      vol_buf = (qy(ny-1)-qy(1))*(qz(nz-1)-qz(1))*dx_buf

      Ni_tot_out_buf = nint(nf_init*vol_buf*beta)

      do 10 l = 1,Ni_tot_out_buf

         xp_out_buf(l,1) = qx(1)-(1.0-pad_ranf())*dx_buf
         xp_out_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp_out_buf(l,3) = qz(2)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         
         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
                flg = 1
               vz = v
            endif
 44      continue

         vp_out_buf(l,1) = -vsw + vx 
         vp_out_buf(l,2) = vy 
         vp_out_buf(l,3) = vz 

         E_out_buf(l,1) = 0.0
         E_out_buf(l,2) = 0.0!-vsw*b0(1,1,1,3)
         E_out_buf(l,3) = vsw*b0(1,1,1,2)!b0_init*eoverm

         B_out_buf(l,1) = b0(1,1,1,1)
         B_out_buf(l,2) = b0(1,1,1,2)!b0_init*eoverm
         B_out_buf(l,3) = b0(1,1,1,3)!0.0


 10   continue
         
      
c      m_arr_out_buf(1:Ni_tot_out_buf) = mproton
c      m_arr_out_buf(Ni_tot_out_buf+1:) = m_pu*mproton 
      mrat_out_buf(1:Ni_tot_out_buf) = 1.0
      mrat_out_buf(Ni_tot_out_buf+1:) = 1.0/m_pu       

      return
      end SUBROUTINE part_setup_out_buf
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
      SUBROUTINE exchange_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp_buf(Ni_max_buf,3),
     x     vp_buf(Ni_max_buf,3),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)
      integer Ni_tot_in, Ni_out
     
      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
c      real, dimension(:), allocatable :: out_m_arr
      real, dimension(:), allocatable :: out_mrat
      real, dimension(:), allocatable :: out_beta_p

      in_bounds_buf(1:Ni_tot_buf) = .true.
      in_bounds_buf(Ni_tot_buf+1:) = .false.

      where(xp_buf(1:Ni_tot_buf,1) .lt. qx(nx))
     X     in_bounds_buf(1:Ni_tot_buf) = .false.

      Ni_tot_in = count(in_bounds_buf)
      Ni_out = count(.not.in_bounds_buf(1:Ni_tot_buf))

      allocate(out_xp(Ni_out,3))
      allocate(out_vp(Ni_out,3))
c      allocate(out_m_arr(Ni_out))
      allocate(out_mrat(Ni_out))
      allocate(out_beta_p(Ni_out))
      
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

c      out_m_arr(1:Ni_out) = pack(m_arr_buf(1:Ni_tot_buf), 
c     x                            .not.in_bounds_buf(1:Ni_tot_buf))
      out_mrat(1:Ni_out) = pack(mrat_buf(1:Ni_tot_buf), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))

      out_beta_p(1:Ni_out) = pack(beta_p_buf(1:Ni_tot_buf), 
     x                            .not.in_bounds_buf(1:Ni_tot_buf))

c      m_arr_buf(1:Ni_tot_in) = pack(m_arr_buf(1:Ni_tot_buf), 
c     x                              in_bounds_buf(1:Ni_tot_buf))
      mrat_buf(1:Ni_tot_in) = pack(mrat_buf(1:Ni_tot_buf), 
     x                             in_bounds_buf(1:Ni_tot_buf))
      
      beta_p_buf(1:Ni_tot_in) = pack(beta_p_buf(1:Ni_tot_buf), 
     x                             in_bounds_buf(1:Ni_tot_buf))
      



      do m = 1,3
         xp(Ni_tot+1:Ni_tot+Ni_out,m) = out_xp(:,m)
         vp(Ni_tot+1:Ni_tot+Ni_out,m) = out_vp(:,m)
         vp1(Ni_tot+1:Ni_tot+Ni_out,m) = out_vp(:,m)
      enddo


c      ijkp(Ni_tot+1:Ni_tot+Ni_out,1) = 
c     x          floor(xp(Ni_tot+1:Ni_tot+Ni_out,1)/dx)
cc      wquad(Ni_tot+1:Ni_tot+Ni_out,1) = -1.0
c      ijkp(Ni_tot+1:Ni_tot+Ni_out,2) = 
c     x          floor(xp(Ni_tot+1:Ni_tot+Ni_out,2)/dy)

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


c         k=1
c         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c            ijkp(l,3) = k       !grid
c            k=k+1
c 50      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif

c         beta_p(l) = 1.0    !add as solar wind ions

      enddo

c      m_arr(Ni_tot+1:Ni_tot+Ni_out) = out_m_arr(:)
      mrat(Ni_tot+1:Ni_tot+Ni_out) = out_mrat(:)
      beta_p(Ni_tot+1:Ni_tot+Ni_out) = out_beta_p(:)
      
      do l = Ni_tot+1,Ni_tot+Ni_out 
         do m=1,3
            input_E = input_E + 
     x         0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /(beta*beta_p(l))
c            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
         enddo
      enddo

      Ni_tot = Ni_tot + Ni_out

      Ni_tot_buf = count(in_bounds_buf)

      deallocate(out_xp)
      deallocate(out_vp)
c      deallocate(out_m_arr)
      deallocate(out_mrat)
      deallocate(out_beta_p)


      return
      end SUBROUTINE exchange_ion_half_buf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE move_ion_half_buf(xp_buf,vp_buf,xp,vp,vp1)
c----------------------------------------------------------------------
c      include 'incurv.h'

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
c      real, dimension(:), allocatable :: out_m_arr
      real, dimension(:), allocatable :: out_mrat

      dth = dt/2

c      write(*,*) 'Ni_tot_buf before move...',Ni_tot_buf,dNi_sw
c      stop

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



c add protons

c      do l = Ni_tot_buf+1,Ni_tot_buf+dNi_sw

c            xp_buf(l,1) = (qx(nx)+dx_buf) - 10*pad_ranf()*vsw*dt/2 
c            xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
c            xp_buf(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(2))

c            vth = 0.5*(vth_top + vth_bottom) + 
c     x      0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

c            call maxwl_init(vth,vx,vy,vz)
            
c            vp_buf(l,1) = -vsw + vx
c            vp_buf(l,2) = vy
c            vp_buf(l,3) = vz
            
cc            m_arr_buf(l) = mproton
c            mrat_buf(l) = 1.0
c            beta_p_buf(l) = 1.0
c      enddo

c      Ni_tot_buf = Ni_tot_buf + dNi_sw


cc add He ++

c      Ni_tot_buf_1 = Ni_tot_buf
c      Ni_tot_buf = Ni_tot_buf_1 + f_mq_2*dNi_sw
c      
c      do l = Ni_tot_buf_1+1,Ni_tot_buf

c            xp_buf(l,1) = (qx(nx)+dx_buf) - 10*pad_ranf()*vsw*dt/2 
c            xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
c            xp_buf(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(2))

c            vth = 0.5*(vth_top + vth_bottom) + 
c     x      0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)
            
c            call maxwl_init(vth,vx,vy,vz)

c            vp_buf(l,1) = -vsw + vx
c            vp_buf(l,2) = vy
c            vp_buf(l,3) = vz
c            
cc            m_arr_buf(l) = 2.0*mproton
c            mrat_buf(l) = 1.0/2.0
c            beta_p_buf(l) = b_mp_2

c      enddo

cc add shell distribution

c      Ni_tot_buf_1 = Ni_tot_buf 
      
c      Ni_tot_buf = Ni_tot_buf_1 + f_shl*dNi_sw
      
c      do 69 l = Ni_tot_buf_1+1,Ni_tot_buf
         
c         xp_buf(l,1) = (qx(nx)+dx_buf) - 10*pad_ranf()*vsw*dt/2 
c         xp_buf(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
c         xp_buf(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(2))
         
cc         m_arr_buf(l) = mproton
c         mrat_buf(l) = 1.0
c         beta_p_buf(l) = b_shl
         
c         theta = pad_ranf()*PI
c         phi = pad_ranf()*2*PI
         
c         vp_buf(l,1) = -vsw+vsw*cos(phi)*sin(theta) !+dvx
c         vp_buf(l,2) = vsw*sin(phi)*sin(theta) !+dvz 
c         vp_buf(l,3) = vsw*cos(theta)
         
c 69   enddo



      return
      end SUBROUTINE move_ion_half_buf
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE exchange_ion_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf,xp,vp,vp1)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp_out_buf(Ni_max_buf,3),
     x     vp_out_buf(Ni_max_buf,3),
     x     E_out_buf(Ni_max_buf,3),
     x     B_out_buf(Ni_max_buf,3),
     x     mrat_out_buf(Ni_max_buf),
c     x     m_arr_out_buf(Ni_max_buf),
     x     xp(Ni_max,3),
     x     vp(Ni_max,3),
     x     vp1(Ni_max,3)

      real vpls_out_buf(3)
      integer Ni_tot_in, Ni_out
     
      real, dimension(:,:), allocatable :: out_xp
      real, dimension(:,:), allocatable :: out_vp
c      real, dimension(:), allocatable :: out_m_arr
      real, dimension(:), allocatable :: out_mrat


c move back into main domain

      in_bounds_buf(1:Ni_tot_out_buf) = .true.
      in_bounds_buf(Ni_tot_out_buf+1:) = .false.

      where(xp_out_buf(1:Ni_tot_out_buf,1) .gt. qx(1))
     X     in_bounds_buf(1:Ni_tot_out_buf) = .false.

      Ni_tot_in = count(in_bounds_buf)
      Ni_out = count(.not.in_bounds_buf(1:Ni_tot_out_buf))

      allocate(out_xp(Ni_out,3))
      allocate(out_vp(Ni_out,3))
c      allocate(out_m_arr(Ni_out))
      allocate(out_mrat(Ni_out))
      
      do m = 1,3 

        out_xp(1:Ni_out,m) = pack(xp_out_buf(1:Ni_tot_out_buf,m), 
     x                            .not.in_bounds_buf(1:Ni_tot_out_buf))
        out_vp(1:Ni_out,m) = pack(vp_out_buf(1:Ni_tot_out_buf,m), 
     x                            .not.in_bounds_buf(1:Ni_tot_out_buf))
        xp_out_buf(1:Ni_tot_in,m) = pack(xp_out_buf(1:Ni_tot_out_buf,m), 
     x                               in_bounds_buf(1:Ni_tot_out_buf))
        vp_out_buf(1:Ni_tot_in,m) = pack(vp_out_buf(1:Ni_tot_out_buf,m), 
     x                               in_bounds_buf(1:Ni_tot_out_buf))

      enddo


c      out_m_arr(1:Ni_out) = pack(m_arr_out_buf(1:Ni_tot_out_buf), 
c     x                            .not.in_bounds_buf(1:Ni_tot_out_buf))
      out_mrat(1:Ni_out) = pack(mrat_out_buf(1:Ni_tot_out_buf), 
     x                            .not.in_bounds_buf(1:Ni_tot_out_buf))

c      m_arr_out_buf(1:Ni_tot_in) = pack(m_arr_buf(1:Ni_tot_out_buf), 
c     x                              in_bounds_buf(1:Ni_tot_out_buf))
      mrat_out_buf(1:Ni_tot_in) = pack(mrat_out_buf(1:Ni_tot_out_buf), 
     x                             in_bounds_buf(1:Ni_tot_out_buf))


      do m = 1,3
         xp(Ni_tot+1:Ni_tot+Ni_out,m) = out_xp(:,m)
         vp(Ni_tot+1:Ni_tot+Ni_out,m) = out_vp(:,m)
         vp1(Ni_tot+1:Ni_tot+Ni_out,m) = out_vp(:,m)
      enddo

c      ijkp(Ni_tot+1:Ni_tot+Ni_out,1) = 
c     x          floor(xp(Ni_tot+1:Ni_tot+Ni_out,1)/dx)
cc      wquad(Ni_tot+1:Ni_tot+Ni_out,1) = -1.0
c      ijkp(Ni_tot+1:Ni_tot+Ni_out,2) = 
c     x          floor(xp(Ni_tot+1:Ni_tot+Ni_out,2)/dy)

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


c         k=1
c         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
c            ijkp(l,3) = k       !grid
c            k=k+1
c 50      continue
c         k=ijkp(l,3)
c         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
c            ijkp(l,3) = k+1
c         endif
         
         write(*,*) 'k in...',ijkp(l,3),xp(l,3),Ni_out
      enddo



c      m_arr(Ni_tot+1:Ni_tot+Ni_out) = out_m_arr(:)
      mrat(Ni_tot+1:Ni_tot+Ni_out) = out_mrat(:)
      beta_p(Ni_tot+1:Ni_tot+Ni_out) = 1.0 !need particle history still

      do l = Ni_tot+1,Ni_tot+Ni_out 
         do m=1,3
            input_E = input_E + 
     x         0.5*(mion/mrat(l))*(vp(l,m)*km_to_m)**2 /(beta*beta_p(l))
c            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
         enddo
      enddo


      Ni_tot = Ni_tot + Ni_out

      Ni_tot_out_buf = count(in_bounds_buf)

      deallocate(out_xp)
      deallocate(out_vp)
c      deallocate(out_m_arr)
      deallocate(out_mrat)


c remove completely from out buf after 10 dx.

      in_bounds_buf(1:Ni_tot_out_buf) = .true.
      in_bounds_buf(Ni_tot_out_buf+1:) = .false.

      where(xp_out_buf(1:Ni_tot_out_buf,1) .lt. qx(1)-dx_buf)
     x     in_bounds_buf(1:Ni_tot_out_buf) = .false.

      Ni_tot_in = count(in_bounds_buf)
      Ni_out = count(.not.in_bounds_buf(1:Ni_tot_out_buf))

c      write(*,*) 'Ni out removed past 10 dx...',Ni_out
         
      do m = 1,3
            
        xp_out_buf(1:Ni_tot_in,m) = pack(xp_out_buf(1:Ni_tot_out_buf,m), 
     x           in_bounds_buf(1:Ni_tot_out_buf))
        vp_out_buf(1:Ni_tot_in,m) = pack(vp_out_buf(1:Ni_tot_out_buf,m), 
     x           in_bounds_buf(1:Ni_tot_out_buf))
        E_out_buf(1:Ni_tot_in,m) = pack(E_out_buf(1:Ni_tot_out_buf,m), 
     x           in_bounds_buf(1:Ni_tot_out_buf))
        B_out_buf(1:Ni_tot_in,m) = pack(B_out_buf(1:Ni_tot_out_buf,m), 
     x           in_bounds_buf(1:Ni_tot_out_buf))

      enddo

        mrat_out_buf(1:Ni_tot_in) = 
     x           pack(mrat_out_buf(1:Ni_tot_out_buf), 
     x           in_bounds_buf(1:Ni_tot_out_buf))
c        m_arr_out_buf(1:Ni_tot_in) = 
c     x           pack(m_arr_out_buf(1:Ni_tot_out_buf), 
c     x           in_bounds_buf(1:Ni_tot_out_buf))



      Ni_tot_out_buf = count(in_bounds_buf)


      return
      end SUBROUTINE exchange_ion_out_buf
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE move_ion_out_buf(xp_out_buf,vp_out_buf,E_out_buf,
     x        B_out_buf,mrat_out_buf)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real xp_out_buf(Ni_max_buf,3),
     x     vp_out_buf(Ni_max_buf,3),
     x     E_out_buf(Ni_max_buf,3),
     x     B_out_buf(Ni_max_buf,3),
     x     mrat_out_buf(Ni_max_buf)
c     x     m_arr_out_buf(Ni_max_buf)
c     x     xp(Ni_max,3),
c     x     vp(Ni_max,3),
c     x     vp1(Ni_max,3)
      integer source, dest
      integer reqs(2)
      integer stats(MPI_STATUS_SIZE,2)
      integer stat(MPI_STATUS_SIZE)

      real vpls_out_buf(3)
      integer zindx(1)
      integer Ni_tot_in, Ni_out, Ni_in
     
      real dth                 !half time step

      real v2

      real, dimension(:,:), allocatable :: out_part
      real, dimension(:,:), allocatable :: in_part
      real, dimension(:), allocatable :: out_mass
      real, dimension(:), allocatable :: in_mass

      logical ib_out_buf(Ni_max_buf) 

      dth = dt/2      

c      write(*,*) 'Ni_tot_buf before move...',Ni_tot_buf,dNi_sw
c      stop

      do 10 l=1,Ni_tot_out_buf                   !make 1/2 time step advance

         call push_part_test(E_out_buf(l,:),B_out_buf(l,:),
     x                       xp_out_buf(l,:),vp_out_buf(l,:),
     x                       mrat_out_buf(l),vpls_out_buf)

         vp_out_buf(l,:) = vpls_out_buf(:)

         xp_out_buf(l,1) = xp_out_buf(l,1) + dth*vpls_out_buf(1)
         xp_out_buf(l,2) = xp_out_buf(l,2) + dth*vpls_out_buf(2)
         xp_out_buf(l,3) = xp_out_buf(l,3) + dth*vpls_out_buf(3)

 10      continue

c         write(*,*) 'xp, vp...',xp_out_buf(1,:),vp_out_buf(1,:)

      where (xp_out_buf(:,2) .ge. qy(ny-1))
         xp_out_buf(:,2) = qy(1) + ( xp_out_buf(:,2) - qy(ny-1) )
      endwhere

      where (xp_out_buf(:,2) .le. qy(1)) 
         xp_out_buf(:,2) = qy(ny-1) - (qy(1) - xp_out_buf(:,2))
      endwhere

c -------------------z exchange, up-----------------------------

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      ib_out_buf(1:Ni_tot_out_buf) = .true.
c      ib_out_buf(Ni_tot_out_buf+1:) = .false.


c      where (xp_out_buf(1:Ni_tot_out_buf,3) .gt. qz(nz))
c         ib_out_buf(1:Ni_tot_out_buf)= .false.
c         xp_out_buf(1:Ni_tot_out_buf,3) = qz(2)+
c     x        (xp_out_buf(1:Ni_tot_out_buf,3)-qz(nz))
c      endwhere

c      Ni_tot_in = count(ib_out_buf(1:Ni_tot_out_buf))

c      Ni_out = count(.not.ib_out_buf(1:Ni_tot_out_buf))

cc      write(*,*) 'out buf exchange up....',Ni_out

c      allocate(out_part(Ni_out,3))
c      allocate(out_mass(Ni_out))

c      dest = nbrs(n_up)
c      source = nbrs(n_down)

c      call MPI_SEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(Ni_in, 1, MPI_INTEGER, source, tag,
c     x     cartcomm, stat, ierr)
      
c      allocate(in_part(Ni_in,3))
c      allocate(in_mass(Ni_in))

c      do m = 1,3
c        out_part(1:Ni_out,m) = 
c     x        pack(xp_out_buf(1:Ni_tot_out_buf,m),  
c     x        .not.ib_out_buf(1:Ni_tot_out_buf))
c        xp_out_buf(1:Ni_tot_in,m) = pack(xp_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo




c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      xp_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x        pack(vp_out_buf(1:Ni_tot_out_buf,m), .not.
c     x        ib_out_buf(1:Ni_tot_out_buf))
c         vp_out_buf(1:Ni_tot_in,m) = 
c     x        pack(vp_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      vp_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x        pack(E_out_buf(1:Ni_tot_out_buf,m), .not.
c     x        ib_out_buf(1:Ni_tot_out_buf))
c         E_out_buf(1:Ni_tot_in,m) = pack(E_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      E_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x        pack(B_out_buf(1:Ni_tot_out_buf,m), .not.
c     x        ib_out_buf(1:Ni_tot_out_buf))
c         B_out_buf(1:Ni_tot_in,m) = pack(B_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      B_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


c      out_mass(1:Ni_out) = 
c     x     pack(m_arr_out_buf(1:Ni_tot_out_buf), .not.
c     x     ib_out_buf(1:Ni_tot_out_buf))
c      m_arr_out_buf(1:Ni_tot_in) = 
c     x     pack(m_arr_out_buf(1:Ni_tot_out_buf), 
c     x     ib_out_buf(1:Ni_tot_out_buf))


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      m_arr_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


c      out_mass(1:Ni_out) = 
c     x     pack(mrat_out_buf(1:Ni_tot_out_buf), .not.
c     x     ib_out_buf(1:Ni_tot_out_buf))
c      mrat_out_buf(1:Ni_tot_in) = 
c     x     pack(mrat_out_buf(1:Ni_tot_out_buf), 
c     x     ib_out_buf(1:Ni_tot_out_buf))


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      mrat_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


cc      write(*,*) 'Ni_tot_out_buf...',Ni_tot_out_buf,Ni_out,Ni_in,qz(2)

c      Ni_tot_out_buf = Ni_tot_out_buf - Ni_out + Ni_in

c      deallocate(out_part)
c      deallocate(out_mass)

c      deallocate(in_part)
c      deallocate(in_mass)





c -------------------z exchange, down-----------------------------

c      call MPI_Barrier(MPI_COMM_WORLD,ierr)

c      ib_out_buf(1:Ni_tot_out_buf) = .true.
c      ib_out_buf(Ni_tot_out_buf+1:) = .false.



c      where (xp_out_buf(1:Ni_tot_out_buf,3) .le. qz(2))
c         ib_out_buf(1:Ni_tot_out_buf)= .false.
c         xp_out_buf(1:Ni_tot_out_buf,3) = qz(nz)-
c     x        (qz(2)-xp_out_buf(1:Ni_tot_out_buf,3))
c      endwhere

c      Ni_tot_in = count(ib_out_buf(1:Ni_tot_out_buf))

c      Ni_out = count(.not.ib_out_buf(1:Ni_tot_out_buf))


c      allocate(out_part(Ni_out,3))
c      allocate(out_mass(Ni_out))

c      dest = nbrs(n_down)
c      source = nbrs(n_up)

c      call MPI_SEND(Ni_out, 1, MPI_INTEGER, dest, tag, 
c     x     cartcomm, ierr)
c      call MPI_RECV(Ni_in, 1, MPI_INTEGER, source, tag,
c     x     cartcomm, stat, ierr)
      
c      allocate(in_part(Ni_in,3))
c      allocate(in_mass(Ni_in))

c      do m = 1,3
c        out_part(1:Ni_out,m) = 
c     x        pack(xp_out_buf(1:Ni_tot_out_buf,m), .not. 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c        xp_out_buf(1:Ni_tot_in,m) = pack(xp_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo

cc      write(*,*) 'out_part...',out_part(1:Ni_out,1)

c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      xp_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)

c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x        pack(vp_out_buf(1:Ni_tot_out_buf,m), .not.
c     x        ib_out_buf(1:Ni_tot_out_buf))
c         vp_out_buf(1:Ni_tot_in,m) = 
c     x        pack(vp_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      vp_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x        pack(E_out_buf(1:Ni_tot_out_buf,m), .not.
c     x        ib_out_buf(1:Ni_tot_out_buf))
c         E_out_buf(1:Ni_tot_in,m) = pack(E_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      E_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


c      do m = 1,3
c         out_part(1:Ni_out,m) = 
c     x        pack(B_out_buf(1:Ni_tot_out_buf,m), .not.
c     x        ib_out_buf(1:Ni_tot_out_buf))
c         B_out_buf(1:Ni_tot_in,m) = pack(B_out_buf(1:Ni_tot_out_buf,m), 
c     x        ib_out_buf(1:Ni_tot_out_buf))
c      enddo


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      B_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in,:) = in_part(:,:)


c      out_mass(1:Ni_out) = 
c     x     pack(m_arr_out_buf(1:Ni_tot_out_buf), .not.
c     x     ib_out_buf(1:Ni_tot_out_buf))
c      m_arr_out_buf(1:Ni_tot_in) = 
c     x     pack(m_arr_out_buf(1:Ni_tot_out_buf), 
c     x     ib_out_buf(1:Ni_tot_out_buf))


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      m_arr_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


c      out_mass(1:Ni_out) = 
c     x     pack(mrat_out_buf(1:Ni_tot_out_buf), .not.
c     x     ib_out_buf(1:Ni_tot_out_buf))
c      mrat_out_buf(1:Ni_tot_in) = 
c     x     pack(mrat_out_buf(1:Ni_tot_out_buf), 
c     x     ib_out_buf(1:Ni_tot_out_buf))


c      call MPI_ISEND(out_part, 3*Ni_out, MPI_REAL, dest, tag, 
c     x     cartcomm, reqs(1), ierr)
c      call MPI_IRECV(in_part, 3*Ni_in, MPI_REAL, source, tag,
c     x     cartcomm, reqs(2), ierr)
      
c      call MPI_WAITALL(2, reqs, stats, ierr)

c      mrat_out_buf(Ni_tot_in+1:Ni_tot_in+Ni_in) = in_mass(:)


c      Ni_tot_out_buf = Ni_tot_out_buf - Ni_out + Ni_in

c      deallocate(out_part)
c      deallocate(out_mass)

c      deallocate(in_part)
c      deallocate(in_mass)

      where (xp_out_buf(:,3) .ge. qz(nz))
         xp_out_buf(:,3) = qz(2) + ( xp_out_buf(:,3) - qz(nz) )
      endwhere

      where (xp_out_buf(:,3) .le. qz(2)) 
         xp_out_buf(:,3) = qz(nz) - (qz(2) - xp_out_buf(:,3))
      endwhere


      return
      end SUBROUTINE move_ion_out_buf
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE push_part_test(E,B,x,v,mr,vpls_out_buf)
c----------------------------------------------------------------------
c      include 'incurv.h'

      real E(3)
      real B(3)
      real x(3)
      real v(3)
      real mr
      real vpls_out_buf(3)
      real vmin(3),vpls(3)
      real vmin_x_B(3),vmin_dot_B
      real aa(3),bb(3),cc(3)

      real a1,a2,a3,a_d
      real Bx, By, Bz
      real B2
      real dt2,dth

      dth = dt/2

c      write(*,*) 'part push 1...',v(:)

c      do m=1,3
c         aa(m) = -v(m)
c         bb(m) = B(m)                   
c      enddo
      
c      cc(1) = aa(2)*bb(3) - aa(3)*bb(2) !do cross product
c      cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
c      cc(3) = aa(1)*bb(2) - aa(2)*bb(1)
      
c         E(:)= E(:)*mr 
      
c         write(*,*) 'part push...',E(3),B(2),
c     x               x(1),v(1),mr


      do m=1,3
         vmin(m) = v(m) + 0.5*dth*E(m)*mr
      enddo

      vmin_x_B(1) = vmin(2)*B(3)*mr- !O_to_Ba - 
     x     vmin(3)*B(2)*mr      !O_to_Ba
      vmin_x_B(2) = vmin(3)*B(1)*mr - !O_to_Ba - 
     x     vmin(1)*B(3)*mr      !O_to_Ba
      vmin_x_B(3) = vmin(1)*B(2)*mr - !O_to_Ba -
     x     vmin(2)*B(1)*mr      !O_to_Ba
      
      vmin_dot_B = vmin(1)*B(1)*mr + !O_to_Ba +
     x     vmin(2)*B(2)*mr +    !O_to_Ba +
     x     vmin(3)*B(3)*mr      !O_to_Ba
      
      Bx = B(1)*mr              !O_to_Ba
      By = B(2)*mr              !O_to_Ba
      Bz = B(3)*mr              !O_to_Ba
      
      B2 = Bx**2 + By**2 + Bz**2
      dt2 = dth**2
      
      a_d = 1 + (B2*dt2/4.0)
      a1 = (1 - (B2*dt2/4.0)) / a_d
      a2 = dth / a_d
      a3 = 0.5*dt2 / a_d

      do m=1,3
         vpls(m) = a1*vmin(m) + a2*vmin_x_B(m) + 
     x        a3*vmin_dot_B*B(m)*mr !O_to_Ba
      enddo

      vpls_out_buf(:) = vpls(:) + 0.5*dth*E(:)*mr

c      write(*,*) 'part push 2...',vpls_out_buf(:)

      return
      end SUBROUTINE push_part_test
c----------------------------------------------------------------------

      end MODULE gutsp_buf

