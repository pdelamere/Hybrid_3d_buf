
      MODULE gutsf

      USE global
      USE boundary
      USE grid_interp
      USE iso_fortran_env, only: error_unit
      
      contains


c----------------------------------------------------------------------
      SUBROUTINE f_update_tlev(b1,b12,b1p2,bt,b0)
c----------------------------------------------------------------------
      
      real b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b0(nx,ny,nz,3)

      b12 = b1
      b1 = b1p2
      bt = b0 + b1
      return
      end SUBROUTINE f_update_tlev
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE crossf2(aa, btc, cc)
c The cross product is formed at the main cell center.  aa and btc must
c be given already extrapolated to the main cell center.
c----------------------------------------------------------------------
      !include 'incurv.h'

      real aa(nx,ny,nz,3)        !main cell contravarient vector 
      real btc(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)

      real cus(ny,nz,3)

      real ax,ay,az,bx,by,bz    !dummy vars

      integer i,j,k

      call periodic(aa)
      call periodic(btc)

      ct(:,:,:,1) = aa(:,:,:,2)*btc(:,:,:,3) - aa(:,:,:,3)*btc(:,:,:,2)
      ct(:,:,:,2) = aa(:,:,:,3)*btc(:,:,:,1) - aa(:,:,:,1)*btc(:,:,:,3)
      ct(:,:,:,3) = aa(:,:,:,1)*btc(:,:,:,2) - aa(:,:,:,2)*btc(:,:,:,1)


c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
      
      do 60 k=2,nz-1
         do 60 j=2,ny-1
            do 60 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

      cus = ct(nx,:,:,:)
      call periodic(cc)

      return
      end SUBROUTINE crossf2
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE curlB(b1,np,aj)
c Calculates curl B / n*alpha.  The resulting "current" is called aj
c which is used in several other places in the code.  This curl is 
c performed on the main cell where B is covarient.  The resulting
c current is main cell contravarient.  Note that dz_cell is used for
c the cell dimensions since dz_grid is not equal to dz_cell on non-
c uniform grid.
c----------------------------------------------------------------------
CVD$R VECTOR
      !include 'incurv.h'

      real b1(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)

      real curl_B(nx,ny,nz,3)      !dummy for holding curl vector
      real ntot(nx,ny,nz,3)        !total density, np + nf
      integer i,j,k

c      call periodic_scalar(np)
c      call periodic_scalar(nf)
      call periodic(b1)
cc     call fix_normal_b(b1)

      do 10 k=2,nz-1   
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               ip = i+1
               jp = j+1
               kp = k+1
               ntot(i,j,k,1) = 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(i,j,k,2) = 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(i,j,k,3) = 0.5*(np(i,j,k)+np(i,j,kp))
 10    continue
      call periodic(ntot)

      do 20 k=2,nz-1   
         do 20 j=2,ny-1
            do 20 i=2,nx-1
               curl_B(i,j,k,1) = (b1(i,j,k,3) - 
     x              b1(i,j-1,k,3))/dy_cell(j) +
     x              (b1(i,j,k-1,2) - b1(i,j,k,2))/dz_cell(k)
               curl_B(i,j,k,2) = (b1(i,j,k,1) - 
     x              b1(i,j,k-1,1))/dz_cell(k) +
     x              (b1(i-1,j,k,3) - b1(i,j,k,3))/dx_cell(i)
               curl_B(i,j,k,3) = (b1(i,j,k,2) - 
     x              b1(i-1,j,k,2))/dx_cell(i) + 
     x              (b1(i,j-1,k,1) - b1(i,j,k,1))/dy_cell(j)
 20            continue
      call periodic(curl_B)

      if(any(ntot .le. nf_init/50.0)) then
          write(error_unit,*) "Density capped in curlB"
      endif
      where(ntot .gt. nf_init/50.0)
      aj = curl_B/(ntot*alpha)
      elsewhere
      aj = curl_B/((nf_init/50.0)*alpha)
      endwhere
      

      return
      end SUBROUTINE curlB
c----------------------------------------------------------------------




c----------------------------------------------------------------------
      SUBROUTINE curlE2(E,curl_E)
c E is dual cell covarient, and curl_E will be returned as main
c cell covarient...as all magnetic fields are.  All i,j,k exclude
c boundaries.  Boundaries are taken care of in main fluid code.
c----------------------------------------------------------------------
c      include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges
      integer i,j,k

c      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               lx = qx(i+1) - qx(i)
               ly = qy(j+1) - qy(j)
               lz = qz(k+1) - qz(k)

               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/ly) - (E(i,j,k,3)/ly)
     x                       + (E(i,j,k,2)/lz) - (E(i,j,k+1,2)/lz)
               curl_E(i,j,k,2) =  (E(i,j,k,3)/lx) - (E(i+1,j,k,3)/lx)
     x                       + (E(i,j,k+1,1)/lz) - (E(i,j,k,1)/lz)
               curl_E(i,j,k,3) =  (E(i,j,k,1)/ly) - (E(i,j+1,k,1)/ly)
     x                       + (E(i+1,j,k,2)/lx) - (E(i,j,k,2)/lx)

 10          continue

c      call periodic(curl_E)

      return
      end SUBROUTINE curlE2
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE curlE(E,curl_E)
c E is dual cell covarient, and curl_E will be returned as main
c cell covarient...as all magnetic fields are.  All i,j,k exclude
c boundaries.  Boundaries are taken care of in main fluid code.
c----------------------------------------------------------------------
      !include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges
      integer i,j,k

c      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

c               lx = qx(i+1) - qx(i)
c               ly = qy(j+1) - qy(j)
c               lz = qz(k+1) - qz(k)

c               lx = dx_grid(i)
c               ly = dy_grid(j)
c               lz = dz_grid(k)

c               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/dy_grid(j)) - 
c     x              (E(i,j,k,3)/dy_grid(j))
c     x              + (E(i,j,k,2)/dz_grid(k)) - 
c     x              (E(i,j,k+1,2)/dz_grid(k))
c               curl_E(i,j,k,2) =  (E(i,j,k,3)/dx_grid(i)) - 
c     x              (E(i+1,j,k,3)/dx_grid(i))
c     x              + (E(i,j,k+1,1)/dz_grid(k)) - 
c     x              (E(i,j,k,1)/dz_grid(k))
c               curl_E(i,j,k,3) =  (E(i,j,k,1)/dy_grid(j)) - 
c     x              (E(i,j+1,k,1)/dy_grid(j))
c     x              + (E(i+1,j,k,2)/dx_grid(i)) - 
c     x              (E(i,j,k,2)/dx_grid(i))


               curl_E(i,j,k,1) =  (E(i,j+1,k,3)- 
     x              E(i,j,k,3))/dy_grid(j)
     x              + (E(i,j,k,2)- 
     x              E(i,j,k+1,2))/dz_grid(k)
               curl_E(i,j,k,2) =  (E(i,j,k,3) - 
     x              E(i+1,j,k,3))/dx_grid(i)
     x              + (E(i,j,k+1,1) - 
     x              E(i,j,k,1))/dz_grid(k)
               curl_E(i,j,k,3) =  (E(i,j,k,1) - 
     x              E(i,j+1,k,1))/dy_grid(j)
     x              + (E(i+1,j,k,2) - 
     x              E(i,j,k,2))/dx_grid(i)



 10          continue

c      call periodic(curl_E)

      return
      end SUBROUTINE curlE
c----------------------------------------------------------------------




cc----------------------------------------------------------------------
c      SUBROUTINE get_um_dot_BB(u,b,cc)
cc uf and btmf are gathered at main cell center and uf.B*B 
cc calculated.  Result returned to main cell contravarient
cc postion.
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
c     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
c     x     cc(nx,ny,nz,3)   !(uf.B)*B

c      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
c      real temp                !used to vectorize loop
cc      real ct(nx,ny,nz,3)      !result are main cell center
c      real udotb               !u dot b
c      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

c! first gather everything at center

c      call periodic(u)
c      call periodic(b)
cc      call fix_normal_b(b)

c      do 5 k=1,nz
c         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
c 5       continue

c      do 10 k=2,nz-1      
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1

c               im = i-1
c               jm = j-1     
c               km = k-1

c               ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
c               bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))

c               uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
c               by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))

c               uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
c               bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)

c               udotb = ux*bx + uy*by + uz*bz

c               ct(i,j,k,1) = udotb*bx
c               ct(i,j,k,2) = udotb*by
c               ct(i,j,k,3) = udotb*bz

c 10            continue

c      call periodic(ct)


cc extrapolate back to main cell contravarient positions.
cc ...just average across cells.

c      do 60 k=2,nz-1
c         do 60 j=2,ny-1
c            do 60 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
c               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
c               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

c 60            continue

c      call periodic(cc)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_um_dot_BB_old(u,b,cc)
cc uf and btmf are gathered at main cell center and uf.B*B 
cc calculated.  Result returned to main cell contravarient
cc postion.
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
c     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
c     x     cc(nx,ny,nz,3)   !(uf.B)*B

c      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
c      real temp                !used to vectorize loop
c      real ct(nx,ny,nz,3)      !result are main cell center
c      real udotb               !u dot b
c      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

c! first gather everything at center

c      call periodic(u)
c      call periodic(b)
cc      call fix_normal_b(b)

c      do 5 k=1,nz
c         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
c 5       continue

c      do 10 i=2,nx-1      
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1

c               im = i-1
c               jm = j-1     
c               km = k-1

c               if (i .eq. 2) then 
c                  ux = 2.0*u(2,j,k,1) - 
c     x                 0.5*(u(3,j,k,1) + u(2,j,k,1))
c                  bx = 2.0*b(2,j,k,1) - 
c     x                 0.5*(b(3,j,k,1) + b(2,j,k,1))
c               else
c                  ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
c                  bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))
c               endif

c               if (j .eq. 2) then 
c                  uy = 2.0*u(i,2,k,2) - 
c     x                 0.5*(u(i,3,k,2) + u(i,2,k,2))
c                  by = 2.0*b(i,2,k,2) - 
c     x                 0.5*(b(i,3,k,2) + b(i,2,k,2))
c               else
c                  uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
c                  by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))
c               endif

c               if (k .eq. 2) then
c                  uz = 2.0*u(i,j,2,3) - 
c     x                 zfrc(k)*(u(i,j,3,3) - u(i,j,2,3)) + u(i,j,2,3) 
cc                  uz = 2.0*u(i,j,2,3) - 
cc     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
c                  bz = b(i,j,2,3)
cc                  bz = 2.0*b(i,j,2,3) - 
cc     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
c               else
c                  uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
c                  bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)
cc                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
cc                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))            
c               endif

c               udotb = ux*bx + uy*by + uz*bz

c               ct(i,j,k,1) = udotb*bx
c               ct(i,j,k,2) = udotb*by
c               ct(i,j,k,3) = udotb*bz

c 10            continue

c      call periodic(ct)

cc extrapolate back to main cell contravarient positions.
cc ...just average across cells.

c      do 60 i=2,nx
c         do 60 j=2,ny
c            do 60 k=2,nz

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               if (i .eq. nx-1) then 
c                  cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
c     x                           0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
c               else
c                  cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
c               endif

c               if (j .eq. ny-1) then 
c                  cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
c     x                           0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
c               else
c                  cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
c               endif
                  
c               if (k .eq. nz-1) then
cc                  temp = 2.0*ct(i,j,nz,3) - 
cc     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
c                   temp = 0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)) +
c     x                   (2.0*dz_cell(nz)/dz_grid(nz))*(ct(i,j,nz,3) -
c     x                    0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)))
c               else
c                  temp = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
c               endif             
  
c               cc(i,j,k,3) = temp

c 60            continue

c      call periodic(cc)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_ugradu_Lax(uf,ugradu,delta_t)
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real uf(nx,ny,nz,3),
c     x     ugradu(nx,ny,nz,3)

c      real ufc(nx,ny,nz,3)
c      real ax1,ax2,ay1,ay2,az1,az2       !misc const
c      real u1,u2,u3                      !temp vars

c      parameter(ad = 0.001)                 !coefficient to add extra
                                         !diffusion
c      call periodic(uf)

c      call face_to_center(uf,ufc)

c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               ip=i+1
c               jp=j+1
c               kp=k+1
c               im=i-1
c               jm=j-1
c               km=k-1

cc xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c               ax1 = 0.5*ufc(i,j,k,1)/dx
c               ax2 = ad*abs(ufc(ip,j,k,1) - ufc(im,j,k,1))
c               ay1 = 0.5*ufc(i,j,k,2)/dy
c               ay2 = ad*abs(ufc(ip,j,k,2) - ufc(im,j,k,2))
c               u1 = ax1*(ufc(ip,j,k,1) - ufc(im,j,k,1)) - 
c     x              ax2*(ufc(im,j,k,1) - 2.0*ufc(i,j,k,1) +
c     x                      ufc(ip,j,k,1))
c               u2 = ay1*(ufc(i,jp,k,1) - ufc(i,jm,k,1)) - 
c     x              ay2*(ufc(i,jm,k,1) - 2.0*ufc(i,j,k,1) +
c     x                      ufc(i,jp,k,1)) 
c               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
c               az2 = ad*abs(ufc(ip,j,k,3) - ufc(im,j,k,3))
c               u3 = az1*(ufc(i,j,kp,1)-ufc(i,j,km,1)) -
c     x              az2*(ufc(i,j,km,1) - 2.0*ufc(i,j,k,1) +
c     x                   ufc(i,j,kp,1))
c               ct(i,j,k,1) = u1 + u2 + u3

cc yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

c               ax1 = 0.5*ufc(i,j,k,1)/dx
c               ax2 = ad*abs(ufc(i,jp,k,1) - ufc(i,jm,k,1))
c               ay1 = 0.5*ufc(i,j,k,2)/dy
c               ay2 = ad*abs(ufc(i,jp,k,2) - ufc(i,jm,k,2))
c               u1 = ax1*(ufc(ip,j,k,2) - ufc(im,j,k,2)) - 
c     x              ax2*(ufc(im,j,k,2) - 2.0*ufc(i,j,k,2) +
c     x                      ufc(ip,j,k,2))
c               u2 = ay1*(ufc(i,jp,k,2) - ufc(i,jm,k,2)) - 
c     x              ay2*(ufc(i,jm,k,2) - 2.0*ufc(i,j,k,2) +
c     x                      ufc(i,jp,k,2)) 
c               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
c               az2 = ad*abs(ufc(i,jp,k,3) - ufc(i,jm,k,3))
c               u3 = az1*(ufc(i,j,kp,2)-ufc(i,j,km,2)) -
c     x              az2*(ufc(i,j,km,2) - 2.0*ufc(i,j,k,2) +
c     x                   ufc(i,j,kp,2))
c               ct(i,j,k,2) = u1 + u2 + u3

cc zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
c               ax1 = 0.5*ufc(i,j,k,1)/dx
c               ax2 = ad*abs(ufc(i,j,kp,1) - ufc(i,j,km,1))
c               ay1 = 0.5*ufc(i,j,k,2)/dy
c               ay2 = ad*abs(ufc(i,j,kp,2) - ufc(i,j,km,2))
c               u1 = ax1*(ufc(ip,j,k,3) - ufc(im,j,k,3)) - 
c     x              ax2*(ufc(im,j,k,3) - 2.0*ufc(i,j,k,3) +
c     x                      ufc(ip,j,k,3))
c               u2 = ay1*(ufc(i,jp,k,3) - ufc(i,jm,k,3)) - 
c     x              ay2*(ufc(i,jm,k,3) - 2.0*ufc(i,j,k,3) +
c     x                      ufc(i,jp,k,3)) 
c               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
c               az2 = ad*abs(ufc(i,j,kp,3) - ufc(i,j,km,3))
c               u3 = az1*(ufc(i,j,kp,3)-ufc(i,j,km,3)) -
c     x              az2*(ufc(i,j,km,3) - 2.0*ufc(i,j,k,3) +
c     x                   ufc(i,j,kp,3))
c               ct(i,j,k,3) = u1 + u2 + u3

c 10            continue

c      call periodic(ct)

cc interpolate back to contravarient positions.

c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1
c               ip = i+1
c               jp = j+1
c               kp = k+1
c               ugradu(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(ip,j,k,1))
c               ugradu(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,jp,k,2))
c               ugradu(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,kp,3))
c 20         continue

c       call periodic(ugradu)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_Ef(Ef,aj,np,nf,up,uf,btmf,nu,ugradu,delta_t,
c     x                  gradPf)
cc Need to treat boundaries separately!!
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     aj(nx,ny,nz,3),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
c     x     nu(nx,ny,nz),
c     x     ugradu(nx,ny,nz,3),
c     x     gradPf(nx,ny,nz,3)
 
c      real ntot(3)                !total plasma density
c      real fnp(3)                 !fraction, np/n
c      real aac(3),bbc(3),ccc(3)
c      real cc(nx,ny,nz,3)

c      call periodic_scalar(np)
c      call periodic_scalar(nf)

c      do 10 k=2,nz-1 
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

cc               if (ip .gt. nx) then ip = nx
cc               if (jp .gt. ny) then jp = ny
cc               if (kp .gt. nz) then kp = nz

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
c     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
c     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
c     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

c               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
c               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
c               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)
               
c               do 10 m=1,3
c                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m)
c 10               continue

c      call crossf(a,btmf,c)

c      call get_ugradu_Lax(uf,ugradu,delta_t)

c      do 20 k=2,nz-1 
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

cc               if (ip .gt. nx) then ip = nx
cc               if (jp .gt. ny) then jp = ny
cc               if (kp .gt. nz) then kp = nz

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
c     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
c     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
c     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

c               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
c               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
c               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)

c               do 20 m=1,3
c                  Ef(i,j,k,m) = c(i,j,k,m) - ugradu(i,j,k,m) 
c     x                          + nu(i,j,k)*fnp(m)*up(i,j,k,m)
c     x                          + nuei*aj(i,j,k,m) - gradPf(i,j,k,m)
cc     x                          + etar(i,j,k,m)*aj(i,j,k,m)
cc                  Ef(i,j,k,m) = c(i,j,k,m) + 
cc     x                          nu(i,j,k)*fnp(m)*up(i,j,k,m)
c 20            continue

c      call periodic(Ef)
cc      call fix_tangential_E(Ef)

c      Ef(nx-1:nx,:,:,3) = 0.0
c      Ef(nx-1:nx,:,:,2) = 0.0
c      Ef(nx-1:nx,:,:,1) = 0.0

c      return
c      end
cc----------------------------------------------------------------------



cc----------------------------------------------------------------------
c      SUBROUTINE get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
c     x                            delta_t)
cc This is the heart of the fluid velocity update.  It solves eqn. 18
cc (Dan's paper) for uf+
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3), !using particle velocities at t level n-1/2
c     x     nu(nx,ny,nz),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     uplus(nx,ny,nz,3),
c     x     uminus(nx,ny,nz,3)

c      real a1,b1,b2,b3      !4 coefficients (see notebook solution)
c      real PP,QQ            !intermediate variables for a1,b1,b2,b3
c      real eta              !intermediate variable for PP,QQ
cc      real B(3)             !B for cross product call
cc      real Bsqrd                     !B*B
c      real um_x_B(nx,ny,nz,3)        !uf- X B
c      real um_dot_BB(nx,ny,nz,3)     !uf- . B
c      real ntot(3)                   !total density np + nf
c      real npave(3)
c      real btc(nx,ny,nz,3)
c      real bsqrd(nx,ny,nz),bsq(3)

c      do 10 i=2,nx-1    
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1
c               do 10 m=1,3
c                  uminus(i,j,k,m) = uf(i,j,k,m) + 
c     x                              0.5*delta_t*Ef(i,j,k,m)
c 10            continue

c      call crossf(uminus, btmf, um_x_B)
c      call get_um_dot_BB(uminus , btmf, um_dot_BB)

c      call face_to_center(btmf,btc)

c      do 15 k=2,nz-1 
c         do 15 j=2,ny-1
c            do 15 i=2,nx-1
c               bsqrd(i,j,k) =  btc(i,j,k,1)**2 + btc(i,j,k,2)**2 + 
c     x               btc(i,j,k,3)**2
c 15            continue

c      call periodic_scalar(bsqrd)
              
c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

cc               if (ip .gt. nx) then ip = nx
cc               if (jp .gt. ny) then jp = ny
cc               if (kp .gt. nz) then kp = nz

c               bsq(1) = 0.5*(bsqrd(i,j,k) + bsqrd(ip,j,k))
c               bsq(2) = 0.5*(bsqrd(i,j,k) + bsqrd(i,jp,k))
c               bsq(3) = 0.5*(bsqrd(i,j,k) + bsqrd(i,j,kp))

c               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
c               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
c               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + npave(1)
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + npave(2)
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + npave(3)

c               do 20 m=1,3

c                  eta = npave(m)*delta_t/(2.0*ntot(m))
c                  QQ = eta/(1.0+nu(i,j,k)*eta)
c                  PP = (1.0-nu(i,j,k)*eta)/(1.0+nu(i,j,k)*eta)
c                  a1 = 1.0/(1.0 + QQ*QQ*bsq(m))
c                  b1 = PP - (QQ*QQ*bsq(m))
c                  b2 = (QQ*PP) + QQ
c                  b3 = (QQ*QQ*PP) + (QQ*QQ)

c                  uplus(i,j,k,m) = a1*(b1*uminus(i,j,k,m) + 
c     x                             b2*um_x_B(i,j,k,m) + 
c     x                             b3*um_dot_BB(i,j,k,m))

c 20            continue

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf,uplus, 
c     x                      uminus,ugradu,up,gradP,nuin,bdp,pf1)
cc Calculate the fluid velocity, uf,  at the new time step and replace
cc uf1 with the new value, uf, in preparation for the next time step.
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     b0(ny),
c     x     b1(nx,ny,nz,3),
c     x     b12(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
c     x     ufp2(nx,ny,nz,3),
c     x     nu(nx,ny,nz),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     uplus(nx,ny,nz,3), 
c     x     uminus(nx,ny,nz,3),
c     x     ugradu(nx,ny,nz,3),
c     x     up(nx,ny,nz,3),
cc     x     gradP(nx,ny,nz,3),
c     x     nuin(nx,ny,nz),
c     x     bdp(nx,ny,nz,3),
c     x     pf1(nx,ny,nz)

c      real b1h(nx,ny,nz,3)
c      real bth(nx,ny,nz,3)
c      real btmfh(nx,ny,nz,3)
c      real ajh(nx,ny,nz,3)
c      real gradPf(nx,ny,nz,3)

c      real m_den
c      real delta_t

c      delta_t = 2.0*dtsub

c      do k=2,nz-1
c         do j=2,ny-1
c            do i=2,nx-1
c               m_den = mO*0.5*(nf(i+1,j,k)+nf(i,j,k))
c               gradPf(i,j,k,1) = (pf1(i+1,j,k)-pf1(i,j,k))/(dx*m_den)
c               m_den = mO*0.5*(nf(i,j+1,k)+nf(i,j,k))
c               gradPf(i,j,k,2) = (pf1(i,j+1,k)-pf1(i,j,k))/(dy*m_den)
c               m_den = mO*0.5*(nf(i,j,k+1)+nf(i,j,k))
c               gradPf(i,j,k,3) = (pf1(i,j,k+1)-
c     x                             pf1(i,j,k))/(dz_grid(k)*m_den)
c            enddo
c         enddo
c      enddo


c      do 10 k=1,nz
c         do 10 j=1,ny
c            do 10 i=1,nx
c               b1h(i,j,k,1) = 0.5*(b1(i,j,k,1) + b12(i,j,k,1))
c               b1h(i,j,k,2) = 0.5*(b1(i,j,k,2) + b12(i,j,k,2))
c               b1h(i,j,k,3) = 0.5*(b1(i,j,k,3) + b12(i,j,k,3))
c               bth(i,j,k,1) = b1h(i,j,k,1) + bdp(i,j,k,1)
c               bth(i,j,k,2) = b1h(i,j,k,2) + b0(j) + bdp(i,j,k,2)
c               bth(i,j,k,3) = b1h(i,j,k,3) + bdp(i,j,k,3)
c 10            continue

c      call cov_to_contra(bth,btmfh)
c      call curlB(b1h,nf,np,ajh)

c      call get_Ef(Ef,ajh,np,nf,up,uf,btmfh,nu,ugradu,delta_t,gradPf)
c      call get_uplus_uminus(Ef,btmfh,uf2,nu,np,nf,uplus,uminus,
c     x                      delta_t)

c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1
c               do 20 m=1,3
c                  ufp2(i,j,k,m) = uplus(i,j,k,m) + 
c     x                            0.5*delta_t*Ef(i,j,k,m) !-
cc     x                        0.5*delta_t*nuin(i,j,k)*uplus(i,j,k,m)
c 20            continue

cc      ufp2(nx-1:nx,:,:,1) = -vsw

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
c     x                      ugradu,aj,up,ufp1,gradP,nuin,pf)
cc Calculate the fluid velocity, uf,  at the new time step and replace
cc uf1 with the new value, uf, in preparation for the next time step.
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
c     x     ufp2(nx,ny,nz,3),
c     x     nu(nx,ny,nz),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     uplus(nx,ny,nz,3), 
c     x     uminus(nx,ny,nz,3),
c     x     ugradu(nx,ny,nz,3),
c     x     aj(nx,ny,nz,3),
c     x     up(nx,ny,nz,3),
c     x     ufp1(nx,ny,nz,3),
cc     x     gradP(nx,ny,nz,3),
c     x     nuin(nx,ny,nz),
c     x     pf(nx,ny,nz)
 
c      real m_den
c      real delta_t
c      real gradPf(nx,ny,nz,3)

c      delta_t = dtsub

c      do k=2,nz-1
c         do j=2,ny-1
c            do i=2,nx-1
c               m_den = mO*0.5*(nf(i+1,j,k)+nf(i,j,k))
c               gradPf(i,j,k,1) = (pf(i+1,j,k)-pf(i,j,k))/(dx*m_den)
c               m_den = mO*0.5*(nf(i,j+1,k)+nf(i,j,k))
c               gradPf(i,j,k,2) = (pf(i,j+1,k)-pf(i,j,k))/(dy*m_den)
c               m_den = mO*0.5*(nf(i,j,k+1)+nf(i,j,k))
c               gradPf(i,j,k,3) = (pf(i,j,k+1)-
c     x                             pf(i,j,k))/(dz_grid(k)*m_den)
c            enddo
c         enddo
c      enddo

      
c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               do 10 m=1,3
c                  ufp1(i,j,k,m) = 0.5*(uf(i,j,k,m) + ufp2(i,j,k,m))
c 10              continue

c      call get_Ef(Ef,aj,np,nf,up,ufp1,btmf,nu,ugradu,delta_t,gradPf)
c      call get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
c     x                      delta_t)

c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1
c               do 20 m=1,3
c                  uf2(i,j,k,m) = uf(i,j,k,m)
c                  uf(i,j,k,m) = uplus(i,j,k,m) + 0.5*dtsub*Ef(i,j,k,m)
cc     x                          - 0.5*dtsub*nuin(i,j,k)*uplus(i,j,k,m)
c 20            continue

cc      uf(nx-1:nx,:,:,1) = -vsw

c      return
c      end
cc----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_E(E,b0,bt,aj,up,np,nu)
      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nu(nx,ny,nz)

      real btc(nx,ny,nz,3)
      real aus(ny,nz,3)
      real bus(ny,nz,3)



      real aa(nx,ny,nz,3)
      real us(ny,nz,3)

      call face_to_center(aj,aa)

      a = aa - up

      call edge_to_center(bt,btc)
      call crossf2(a,btc,c)

      E = c + spread(nu, 4, 3)*aj

      return
      end SUBROUTINE get_E
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu)
c Predictor step in magnetic field update.
c----------------------------------------------------------------------
CVD$R VECTOR
c      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)   !curl of E
      real bus(ny,nz,3)
      integer i,j,k,m

c      call cov_to_contra(bt,btmf) 
c      call edge_to_center(bt,btc)
c      call get_E(E,b0,bt,btmf,aj,up,np,nu)  !E at time level m 
      call get_E(E,b0,bt,aj,up,np,nu)  !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
c      call periodic(E)
c      call fix_tangential_E(E)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3

c                  b1p2(i,j,k,m)=lww1*(b12(i+1,j,k,m)+
c     x                 b12(i-1,j,k,m)+
c     x                 b12(i,j+1,k,m)+b12(i,j-1,k,m)+
c     x                 b12(i,j,k+1,m)+b12(i,j,k-1,m))+
c     x                 lww2*b12(i,j,k,m) -
c     x                 2.0*dtsub*curl_E(i,j,k,m)

                  b1p2(i,j,k,m) = b12(i,j,k,m) - 
     x                            2.0*dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)
c      call obstacle_boundary_B(b0,b1p2)
c      call fix_normal_b(b1p2)

      return
      end SUBROUTINE predict_B
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)
c The main feature here is that E must be calculated at time level
c m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
c calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
c m + 1/2, so they are used as is. 
c----------------------------------------------------------------------

      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nu(nx,ny,nz)

      real b1p1(nx,ny,nz,3)   !b1 at time level m + 1/2
      real btp1(nx,ny,nz,3)   !bt at time level m + 1/2
      real btp1mf(nx,ny,nz,3) !btp1 at contravarient position
      real btc(nx,ny,nz,3) 
      real aa(nx,ny,nz,3) 
      real aus(ny,nz,3)
      real bus(ny,nz,3)
    

      real ntot(3)            !total density np + nf
      real fnp(3),fnf(3)      !fraction np and nf of n
      real npave(3)
      integer i,j,k

      do 5 k=1,nz
         do 5 j=1,ny
            do 5 i=1,nx
               btp1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1)) +
     x              b0(i,j,k,1)
               b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
               btp1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2)) + 
     x              b0(i,j,k,2)
               b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
               btp1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3)) + 
     x              b0(i,j,k,3)
               b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
 5             continue

      call curlB(b1p1,np,aj)
      call face_to_center(aj,aa)
      a = aa - up

      call edge_to_center(btp1,btc)

      call crossf2(a,btc,c)
       
      E = c + spread(nu,4,3)*aj

      return
      end SUBROUTINE get_Ep1
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE correct_B(b0,b1,b1p2,E,aj,up,np,nu)
c Corrector step in magnetic field update.
c----------------------------------------------------------------------
CVD$R VECTOR
c      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3),
c     x     bdp(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)            !curl of E
      real bus(ny,nz,3)
      integer i,j,k,m

      call get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)  
                                                   !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
c      call periodic(E)
c      call fix_tangential_E(E)

c      write(*,*) 'E cb...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3

c                  b1p2(i,j,k,m)=lww1*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+
c     x                 b1(i,j+1,k,m)+b1(i,j-1,k,m)+
c     x                 b1(i,j,k+1,m)+b1(i,j,k-1,m))+
c     x                 lww2*b1(i,j,k,m) -
c     x                 dtsub*curl_E(i,j,k,m)

                  b1p2(i,j,k,m) = b1(i,j,k,m) - 
     x                            dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)
c      call obstacle_boundary_B(b0,b1p2)
c      call fix_normal_b(b1p2)

      return
      end SUBROUTINE correct_B
c----------------------------------------------------------------------

c----------------------------------------------------------------      
      SUBROUTINE check_time_step(b0,b1,bt,np,step,error_file)
c----------------------------------------------------------------      
      
      real b0(nx,ny,nz,3)
      real b1(nx,ny,nz,3)
      real bt(nx,ny,nz,3)
      real np(nx,ny,nz)
      real ak, btot, a1, a2, womega, phi, deltat
      integer error_file
      integer step
      integer i,j,k

      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               ak = 2./dx
               btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + 
     x              bt(i,j,k,3)**2)
               a1 = ak**2*Btot/(alpha*(np(i,j,k)))
               a2 = (ak*Btot)**2/(alpha*(np(i,j,k)))
               womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
               phi = womega/ak
               deltat = dx/phi
               if(deltat .le. 2.0*dtsub) then 
                   ! rank plus 1 is proc number
                  write(*,*) 'time stepping error...', my_rank+1
                  write(error_file) step
                  write(error_file) my_rank+1
                  write(error_file) i,j,k
                  write(error_file) qx(i),qy(j),qz(k),gz(k)
                  write(error_file) np
                  write(error_file) b0
                  write(error_file) b1
                  write(error_file) bt

                  dtsub = dtsub/2.0
                  ntf = ntf*2.0
               endif

            enddo
         enddo
      enddo

      if (ntf .gt. 100) then 
          write(error_unit,*) 'Aborting due to excesive subcycling'
          stop
      endif
      
      return
      end subroutine check_time_step
c----------------------------------------------------------------      

      


      end MODULE gutsf













