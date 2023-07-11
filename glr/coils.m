# $Id: coils.m,v 1.5 2003/05/22 16:14:36 bulmer Exp $
 
      subroutine bvcal(bs,zs,ass,als,z,n,zcc,bcc,znorm,bv,bvp,bvpp,b4,b5,z4,z5,nc)
      implicit none

      real bs(*),zs(*),ass(*),als(*),z(*),bv(*),bvp(*),bvpp(*)
      real zcc,bcc,znorm,b4,b5,z4,z5
      integer n,nc,i
      real bcorr,bcorrp,bcorrpp,tanhyp

      call bcccal(bs,zs,ass,als,z,n,bv,bvp,bvpp,b4,b5,z4,z5,nc)

      do i = 1,n
          bcorr = bcc
          bcorrp = 0
          bcorrpp = 0
          tanhyp = tanh((z(i) - zcc)/znorm)
          bcorrp = -bcc*0.5*(1 - tanhyp**2)/znorm
          bcorrpp = bcc*tanhyp*(1 - tanhyp**2)/znorm**2
          bcorr = bcc*0.5*(1 - tanhyp)
          bv(i) = bv(i) + bcorr
          bvp(i) = bvp(i) + bcorrp
          bvpp(i) = bvpp(i) + bcorrpp
          if (z(i) == z4) b4 = b4 + bcorr
          if (z(i) == z5) b5 = b5 + bcorr
      enddo

      return
      end # bvcal

      subroutine bcccal(bs,zs,ass,als,z,n,bv,bvp,bvpp,b4,b5,z4,z5,nc)
      implicit none

      integer n,nc,i,j,i4,i5,is(3)
      real z(n),bv(n),bvp(n),bvpp(n)
      real bs(*),zs(*),ass(*),als(*)
      real b4,b5,z4,z5,det
      real alpha(3,3),ak(3)
      real bfun

      # compute matrix elements
      do j = 1,nc
          do i = 1,nc
              alpha(i,j) = bfun(zs(i),zs(j),als(j),ass(j),0)
          enddo
      enddo
 
      # determinant
      if (nc == 2) alpha(3,3) = 1
      det = alpha(1,1)*(alpha(2,2)*alpha(3,3) - alpha(3,2)*alpha(2,3)) 
     *     - alpha(1,2)*(alpha(2,1)*alpha(3,3) - alpha(3,1)*alpha(2,3)) 
     *     + alpha(1,3)*(alpha(2,1)*alpha(3,2) - alpha(3,1)*alpha(2,2))
 
      # solution
      ak(1) = (bs(1)*(alpha(2,2)*alpha(3,3) - alpha(3,2)*alpha(2,3)) 
     *       -  bs(2)*(alpha(1,2)*alpha(3,3) - alpha(3,2)*alpha(1,3)) 
     *       +  bs(3)*(alpha(1,2)*alpha(2,3) - alpha(2,2)*alpha(1,3))) 
     *       /det
      ak(2) = (-bs(1)*(alpha(2,1)*alpha(3,3) - alpha(3,1)*alpha(2,3)) 
     *       +   bs(2)*(alpha(1,1)*alpha(3,3) - alpha(3,1)*alpha(1,3)) 
     *       -   bs(3)*(alpha(1,1)*alpha(2,3) - alpha(2,1)*alpha(1,3))) 
     *       /det
      ak(3) = (bs(1)*(alpha(1,2)*alpha(2,3) - alpha(1,3)*alpha(2,2)) 
     *       -  bs(2)*(alpha(1,1)*alpha(2,3) - alpha(1,3)*alpha(2,1)) 
     *       +  bs(3)*(alpha(1,1)*alpha(2,2) - alpha(1,2)*alpha(2,1))) 
     *       /det
 
      # fields and derivatives
      do i = 1,n
          bv(i) = 0
          bvp(i) = 0
          bvpp(i) = 0
          do j = 1,nc
              bv(i) = bv(i) + ak(j)*bfun(z(i),zs(j),als(j),ass(j),0)
              bvp(i) = bvp(i) + ak(j)*bfun(z(i),zs(j),als(j),ass(j),1)
              bvpp(i) = bvpp(i) + ak(j)*bfun(z(i),zs(j),als(j),ass(j),2)
          enddo
      enddo
 
      # minima and their positions
      do j = 1,nc
          do i = 2,n
              if (zs(j) < z(i)) break
          enddo
          is(j) = i - 1
      enddo

      if (n == 1) return

      b4 = bv(is(1))
      do i = is(1),is(2)
          if (bv(i) < b4) then
              b4 = bv(i)
              i4 = i
          endif
      enddo
      if (i4 > 0) then
          z4 = z(i4)
      endif

      if (nc == 2) return

      b5 = bv(is(2))
      do i = is(2),is(3)
          if (bv(i) < b5) then
              b5 = bv(i)
              i5 = i
          endif
      enddo
      if (i5 > 0) then
          z5 = z(i5)
      endif

      return
      end # bcccal

      function bfun(z,zx,al,a,ind)
      implicit none

      real bfun,gfun,gfun1,gfun2
      real z,zx,al,a,up,um,t1,t2
      integer ind,ind1

      ind1 = ind + 1
      up = zx + 0.5*al
      um = zx - 0.5*al
      if (ind1 == 1) then
          t1 = gfun(z,up,a)
          t2 = gfun(z,um,a)
      elseif (ind1 == 2) then
          t1 = gfun1(z,up,a)
          t2 = gfun1(z,um,a)
      elseif (ind1 == 3) then
          t1 = gfun2(z,up,a)
          t2 = gfun2(z,um,a)
      endif
      bfun = (t1 - t2)/a**2

      return
      end # bfun

      function gfun(z,u,a)
      implicit none

      real gfun
      real z,u,a,x1,x2

      x1 = u - z
      x2 = sqrt(x1**2 + a**2)
      gfun = x1/x2

      return
      end # gfun

      function gfun1(z,u,a)
      implicit none

      real gfun1
      real z,u,a,x1,x2

      x1 = u - z
      x2 = (x1**2 + a**2)**1.5
      gfun1 = -a**2/x2

      return
      end # gfun1

      function gfun2(z,u,a)
      implicit none

      real gfun2
      real z,u,a,x1,x2

      x1 = u - z
      x2 = (x1**2 + a**2)**2.5
      gfun2 = -3.0*x1*a**2/x2

      return
      end # gfun2
