# $Id: diagno.m,v 1.5 2003/06/04 16:50:00 bulmer Exp $ 

      subroutine diagno
      implicit none

c.....sets up arrays for spatial plots

      Use(ArraySizes)
      Use(Const_1a)
      Use(Const_5)
      Use(Fstor)

      integer np,i,kp1,kp2,kpp,ll,ixpl,j
 
      np=np+1
      do i=2,ix-1
          kp1=i-1+(jx/2-1)*(ix-2)
          kp2=jx/2-1+(jx-2)*(i-2)
          kpp=.5*(1+isw)*kp1+.5*(1-isw)*kp2
          xrspz(i,np)=xro(kpp)
          xflute(i,np)=xro(kpp)*r(i,jx/2)*b(i,jx/2)
      enddo
cjb  xrspz(i,  amd xflute(i,     use only xro, not xio

      do ll=1,2
          ixpl=ix/2**ll
          do j=2,jx-1
              kp1=ixpl-1+(j-2)*(ix-2)
              kp2=j-1+(jx-2)*(ixpl-2)
              kpp=.5*(1+isw)*kp1+.5*(1-isw)*kp2
              xrsppsi(j,np+(ll-1)*nps)=xro(kpp)
              xrsppsi(j,np+(ll-1)*nps)=xro(kpp)
          enddo
      enddo
cjb xrsppsi(j,   also only uses xro, not xio

      npm=np
cjb  npm is in vdf    like counter for num calls diagno
      return
      end # diagno
