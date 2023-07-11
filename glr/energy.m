# $Id: energy.m,v 1.7 2003/06/04 16:51:42 bulmer Exp $

      subroutine energy
      implicit none

c...  energy calculates the total first order energy (tenergy) as
c...  an integral (using trapazoidal quadrature) over the volume,
c...  every nen'th time step, using two consecutive time steps to
c...  determine the value at n+.5 .

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_4)
      Use(Const_5)
      Use(Const_8)
      Use(GEnergy)
      Use(Fstor)

      integer lmove,i,j,k1,k2,k,l
      real tquadoo,tquado,tquad,tquad2oo,tquad2o,tquad2,tquad3oo
      real tquad3o,tquad3,tquad4oo,tquad4o,tquad4,psum,tquad1
       
      ltt=ltt+1
      if ((ltt == 2) .or. (nen == 1)) then
          lmove=jx*ix
          call mymove(drxro(1,1),drxr(1,1),lmove)
          call mymove(drxio(1,1),drxi(1,1),lmove)
 
          tquadoo=tquado
          tquado=tquad
          tquad2oo=tquad2o
          tquad2o=tquad2
          tquad3oo=tquad3o
          tquad3o=tquad3
          tquad4oo=tquad4o
          tquad4o=tquad4
      endif

      do i=2,ix-1
          do j=2,jx-1
              k1=i-1+(j-2)*(ix-2)
              k2=j-1+(jx-2)*(i-2)
              k=0.5*((1+isw)*k1+(1-isw)*k2)
              rxr(j)=r(i,j)*xro(k)
              rxi(j)=r(i,j)*xio(k)
          enddo
          call ddpsi(rxr,drxr,i,ix)
          call ddpsi(rxi,drxi,i,ix)
      enddo
 
c     evaluates z quadrature
 
      do j=2,jx-1
          do i=2,ix-1
              k1=i-1+(j-2)*(ix-2)
              k2=j-1+(jx-2)*(i-2)
              k=0.5*((1+isw)*k1+(1-isw)*k2)
              if ((ltt == 2) .or. (nen == 1))then
                  xrgo(i)=xroo(k)
                  xigo(i)=xioo(k)
              endif
              xrg(i)=xro(k)
              xig(i)=xio(k)
          enddo
          do i=2,ix-1
              dum2(i)=r(i,j)*b(i,j)*xrg(i)
          enddo
          call ddz(dum2,dxrz2)
          do i=2,ix-1
              dum2(i)=r(i,j)*b(i,j)*xig(i)
          enddo
          call ddz(dum2,dxiz2)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*qub(i,j)*dxrz2(i)/((r(i,j)*b(i,j))**2)
          enddo
          psuma(j,1)=psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*qub(i,j)*dxiz2(i)/((r(i,j)*b(i,j))**2)
          enddo
          psuma(j,2)=psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=b(i,j)*drxi(i,j)
          enddo
          call ddz(dum1,dum2)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*qub(i,j)*r(i,j)**2*dum2(i)
          enddo
          psuma(j,3)=psum(dum1,ix)/vpsi(j)/mm**2
          do i=2,ix-1
              dum1(i)=b(i,j)*drxr(i,j)
          enddo
          call ddz(dum1,dum2)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*qub(i,j)*r(i,j)**2*dum2(i)
          enddo
          psuma(j,4)=psum(dum1,ix)/vpsi(j)/mm**2
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*xrg(i)**2*r(i,j)*(qubv(i,j)*rzz(i,j) \
                     +r(i,j)*ering(i,j)/b(i,j))
          enddo
          psuma(j,5)=psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*xig(i)**2*r(i,j)*(qubv(i,j)*rzz(i,j) \
                     +r(i,j)*ering(i,j)/b(i,j))
          enddo
          psuma(j,6)=psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*xrg(i)**2*yyy(i,j)/b(i,j)
          enddo
          psuma(j,7)=mm**2*psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*xig(i)**2*yyy(i,j)/b(i,j)
          enddo
          psuma(j,8)=mm**2*psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*yyy(i,j)*xrg(i)*r(i,j)*drxr(i,j)
          enddo
          psuma(j,9)=-2.0*psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*yyy(i,j)*r(i,j)**2*b(i,j)*drxr(i,j)**2
          enddo
          psuma(j,10)=psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=-(1.0/uuz(i))*yyy(i,j)*xig(i)*r(i,j)*drxi(i,j)
          enddo
          psuma(j,11)=2*psum(dum1,ix)/vpsi(j)
          do i=2,ix-1
              dum1(i)=(1.0/uuz(i))*yyy(i,j)*r(i,j)**2*b(i,j)*drxi(i,j)**2
          enddo
          psuma(j,12)=psum(dum1,ix)/vpsi(j)
          if ((ltt == 2) .or. (nen == 1)) then
              do i=2,ix-1
                  dum1(i)=(1.0/uuz(i))*rho(i,j)*((xrg(i)-xrgo(i))/dt)**2/b(i,j)
              enddo
              psumb(j,1)=psum(dum1,ix)/vpsi(j)
              do i=2,ix-1
                  dum1(i)=(1.0/uuz(i))*rho(i,j)*((xig(i)-xigo(i))/dt)**2/b(i,j)
              enddo
              psumb(j,2)=psum(dum1,ix)/vpsi(j)
              do i=2,ix-1
                  dum1(i)=(1.0/uuz(i))*(r(i,j)*(drxr(i,j) \
                         -drxro(i,j))/dt)**2*b(i,j)*rho(i,j)
              enddo
              psumb(j,3)=psum(dum1,ix)/vpsi(j)/mm**2
              do i=2,ix-1
                  dum1(i)=(1.0/uuz(i))*(r(i,j)*(drxi(i,j) \
                         -drxio(i,j))/dt)**2*b(i,j)*rho(i,j)
              enddo
              psumb(j,4)=psum(dum1,ix)/vpsi(j)/mm**2
          endif
      enddo

      do l=1,12
          quad(l)=psum(psuma(1:jx,l),jx)
      enddo
      tquad=0.0
      do l=1,12
          tquad=tquad+quad(l)*dv*du
      enddo       
      tquad2=0.0
      do l=1,4
          tquad2=tquad2+quad(l)
      enddo
      tquad3=quad(5)+quad(6)
      tquad4=0.0
      do l=7,12
          tquad4=quad(l)+tquad4
      enddo
      if ((ltt == 2) .or. (nen == 1)) then
          tquad1=0.0
          lee=lee+1
          do l=1,4
              quad1(l)=psum(psumb(1:jx,l),jx)
              tquad1=tquad1+quad1(l)
          enddo
          enpot(lee)=(3.0*tquad+6.0*tquado-tquadoo)/8.0
          enkin(lee)=tquad1*dv*du
          tenergy(lee)=abs(enpot(lee)+enkin(lee))
          if (n >= 1 .and. n <= nmax) then
              timengy(lee)=time(n)
          endif
          enbend(lee)=(3.0*tquad2+6.0*tquad2o-tquad2oo)/8.0*du*dv
          encurve(lee)=abs(3.0*tquad3+6.0*tquad3o-tquad3oo)/8.0*du*dv
          tenrel(lee)=abs(tenergy(lee))/(enkin(lee)+abs(enpot(lee)))*2
          enflr(lee)=abs(3.0*tquad4+6.0*tquad4o-tquad4oo)/8.0*du*dv
          ltt=0
      endif
      if (mod(n,iscreen) == 0) then
          call baswline(STDOUT,'end of subr energy')
      endif

      return
      end # energy

      subroutine ddpsi(f1,f2,i,imx)
      implicit none

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_4)

      integer i,j,imx
      real f1(*),f2(imx,*)

      do j=4,jx-1
          f2(i,j-1)=(f1(j)-f1(j-2))/(2.0*dv)*vpsi(j-1)
      enddo
      f2(i,2)=(2.0*f1(2)*vpsih(1)+(f1(3)-f1(2))*vpsih(2))/(dv*2)
      f2(i,jx-1)=(-(1.0-fjrx)*f1(jx-1)*vpsih(jx-1)+(f1(jx-1) \
                -f1(jx-2))*vpsih(jx-2))/(dv*2)
      return
      end # ddpsi

      subroutine ddz(f1,f2)
      implicit none

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_4)

      real f1(*),f2(*)
      integer i

      do i=4,ix-1
          f2(i-1)=((f1(i)-f1(i-2))/(2.0*du)*uuz(i-1))**2
      enddo
      f2(2)=(((1.0-fi1)*f1(2)*uuzh(1)+(f1(3)-f1(2))*uuzh(2))/(du*2))**2
      f2(ix-1)=((-(1.0-fizx)*f1(ix-1)*uuzh(ix-1)+(f1(ix-1)-f1(ix-2)) \
              *uuzh(ix-2))/(du*2))**2
      return
      end # ddz

      function psum(f,k)
      implicit none

      real psum, f(*)
      integer k,l,nq

      nq=k-4
      psum=0.05*(f(2)+f(k-1))
      do l=3,3+nq-1
          psum=psum+f(l)
      enddo

      return
      end # psum
