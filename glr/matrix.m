# $Id: matrix.m,v 1.11 2004/06/29 18:47:45 byers Exp $

      subroutine f1to11
      implicit none

c.....calculates the f1 to f11 functions needed to generate the a and
c.....b matrices . uses the equiibrium quantities r, rho, b, etc.

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_4)
      Use(Const_7)
      Use(Const_8)
      Use(Fstor)

      integer m2,i,j
      real r2,bb,vp,r4,f2t,dppdpsi,dpedpsi

      m2=mm**2
      do i=1,ix
          do j=2,jx
              r2=r(i,j)**2
              bb=b(i,j)
              vp=vpsi(j)
              r4=r2**2
              f1(i,j)=rho(i,j)*bb*r4
              if (j == 1) then
                  f2t=(1.-m2)*rho(i,j)/bb+r2*vp*(rho(i,j+1)-rho(i,j))/dv
              elseif (j == jx) then
                  f2t=(1.-m2)*rho(i,j)/bb+r2*vp*(rho(i,j)-rho(i,j-1))/dv
              else
                  f2t=(1.-m2)*rho(i,j)/bb+r2*vp*(rho(i,j+1)-rho(i,j-1))/(2.*dv)
              endif
              f2(i,j)=f2t/vp
              f3(i,j)=mm*xxx(i,j)*r4*bb
              f4(i,j)=(1.-m2)*mm*xxx(i,j)/bb
              f5(i,j)=-m2*yyy(i,j)*r4*bb
              f7(i,j)=(1.-m2)*(-m2)*yyy(i,j)/(bb*vp)
              g4(i,j)=qub(i,j)*r(i,j)**2
              g3(i,j)=r(i,j)*b(i,j)
              g2(i,j)=qub(i,j)/(r(i,j)*b(i,j))**2
              dppdpsi=dp2dpsi(j)*(abp(i)*b(i,j)**4+bbp(i)*b(i,j)**2+cbp(i)) 
     *                +p2(j)*(4*abp(i)*b(i,j)**3+2*bbp(i)*b(i,j))*dbdpsi(i,j)
              dpedpsi=dp1dpsi(j)*(abf(i)*b(i,j)**4+bbf(i)*b(i,j)**2+cbf(i)) 
     *                +p1(j)*(4*abf(i)*b(i,j)**3+2*bbf(i)*b(i,j))*dbdpsi(i,j)
              ering(i,j)=dppdpsi*dpedpsi
              g1(i,j)=+(mm)**2*r(i,j)*(rzz(i,j)*qubv(i,j) 
     *                +dppdpsi*dpedpsi*r(i,j)/b(i,j))

cjb dppdpsi used only *dpedpsi  in both ering and g1
cjb              so if dpedpsi=0 dont need dppdpsi

              g1(i,j)=g1(i,j)*swg1/vpsi(j)
              g2(i,j)=g2(i,j)*swg2
              g3(i,j)=g3(i,j)*swg3
              g4(i,j)=g4(i,j)*swg4
          enddo
      enddo

c..... fill in edge values
      do i=1,ix
          f1(i,1)=-f1(i,2)
          f2(i,1)=f2(i,2)
          f3(i,1)=-f3(i,2)
          f4(i,1)=f4(i,2)
          f5(i,1)=-f5(i,2)
          f7(i,1)=f7(i,2)
          g4(i,1)=-g4(i,2)
      enddo

      return
      end # f1to11

      subroutine amat
      implicit none
 
c.....calculates the matrix coefficients for a1, a2, a3, b1, b2
c.....in the equation a1*x(n+1)=a2*x(n)+a3*x(n-1)+b1*y(n)+b2*y(n-1) .
c.....uses f1 to f11 from subroutine f1to11 and equilibrium quantities.
 
      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_2)
      Use(Const_4)
      Use(Fstor)
      Use(Matrix)
 
      integer i,j,k,k1,k2,k1i,k1j,k2i,k2j,m
      real gam3,du2,dt2,dv2,vp,bijmh,bijph
      real bip1jph,bip1jmh,bim1jph,bim1jmh,g4iphjph,g4iphjmh,g4imhjph
      real g4imhjmh,g2iphj,g2imhj,f1ijph,f1ijmh,f5ijph
      real f5ijmh,uzbar,f3ijmh,f3ijph,denom,sfi1,sfj1,sfjrx
      real sfizx,fac3,fac4,rbt,rbt1,down2,up2,gm2,up1,down1,gm1,up,down
      real fac5

      gam3=-gam2
      du2=du**2
      dt2=dt**2
      dv2=dv**2
      do i=2,ix-1
          do j=2,jx-1
              k1=i-1+(j-2)*(ix-2)
              k2=j-1+(jx-2)*(i-2)
              k=.5*(1+isw)*k1+.5*(1-isw)*k2
              vp=vpsi(j)
              bijmh=(b(i,j)+b(i,j-1))*.5
              bijph=(b(i,j)+b(i,j+1))*.5
              bip1jph=(b(i+1,j+1)+b(i+1,j))*.5
              bip1jmh=(b(i+1,j-1)+b(i+1,j))*.5
              bim1jph=(b(i-1,j+1)+b(i-1,j))*.5
              bim1jmh=(b(i-1,j-1)+b(i-1,j))*.5
              g4iphjph=(g4(i+1,j+1)+g4(i,j)+g4(i+1,j)+g4(i,j+1))*.25*uuzh(i)
              g4iphjmh=(g4(i+1,j-1)+g4(i,j)+g4(i+1,j)+g4(i,j-1))*.25*uuzh(i)
              g4imhjph=(g4(i-1,j+1)+g4(i,j)+g4(i-1,j)+g4(i,j+1))*.25*uuzh(i)
              g4imhjmh=(g4(i-1,j-1)+g4(i,j)+g4(i-1,j)+g4(i,j-1))*.25*uuzh(i)
              g2iphj=(g2(i+1,j)+g2(i,j))*.5*uuzh(i)
              g2imhj=(g2(i-1,j)+g2(i,j))*.5*uuzh(i-1)
 
              f1ijph=(f1(i,j)+f1(i,j+1))*.5*vpsih(j)
              f1ijmh=(f1(i,j)+f1(i,j-1))*.5*vpsih(j-1)
              f5ijph=(f5(i,j)+f5(i,j+1))*.5*vpsih(j)
              f5ijmh=(f5(i,j)+f5(i,j-1))*.5*vpsih(j-1)
 
              if (j <= 2) then
                  f1ijmh=0.
              endif
              uzbar=-uuz(i)*r(i,j)/(du2*dv2)
              a1(k,1)=-gam1*bim1jmh*g4imhjmh*bijmh*uzbar*vpsih(j-1) 
     *                *r(i-1,j-1)
cjb a1(k,1)>>> i-1,j-1

              a2(k,1)=-gam2*bim1jmh*g4imhjmh*bijmh*uzbar*vpsih(j-1) 
     *                *r(i-1,j-1)
              a3(k,1)=-gam3*bim1jmh*g4imhjmh*bijmh*uzbar*vpsih(j-1) 
     *                *r(i-1,j-1)
              
              a1(k,2)=-f1ijmh/((dt*dv)**2)+gam1*(f5ijmh/dv2+bijmh**2*uzbar 
     *                *vpsih(j-1)*r(i,j-1)*(g4imhjmh+g4iphjmh))
cjb a1(k,2)>>> i,j-1

              a2(k,2)=-f1ijmh/((dt*dv)**2)+gam2*(f5ijmh/dv2+bijmh**2*uzbar 
     *                *vpsih(j-1)*r(i,j-1)*(g4imhjmh+g4iphjmh))
              a3(k,2)=-f1ijmh/((dt*dv)**2)+gam3*(f5ijmh/dv2+bijmh**2*uzbar 
     *                *vpsih(j-1)*r(i,j-1)*(g4imhjmh+g4iphjmh))
              
              a1(k,3)=-gam1*bip1jmh*g4iphjmh*bijmh*uzbar*vpsih(j-1)*r(i+1,j-1)
cjb a1(k,3)>>> i+1,j-1

              a2(k,3)=-gam2*bip1jmh*g4iphjmh*bijmh*uzbar*vpsih(j-1)*r(i+1,j-1)
              a3(k,3)=-gam3*bip1jmh*g4iphjmh*bijmh*uzbar*vpsih(j-1)*r(i+1,j-1)

              a1(k,4)=gam1*((bim1jmh*g4imhjmh*vpsih(j-1)*bijmh+bim1jph*g4imhjph 
     *                *vpsih(j)*bijph)*uzbar*r(i-1,j)+mm**2*b(i,j)*uzbar 
     *                *g2imhj*g3(i-1,j)*dv2/vpsi(j))
cjb a1(k,4)>>> i-1,j

              a2(k,4)=gam2*((bim1jmh*g4imhjmh*vpsih(j-1)*bijmh+bim1jph*g4imhjph 
     *                *vpsih(j)*bijph)*uzbar*r(i-1,j)+mm**2*b(i,j)*uzbar 
     *                *g2imhj*g3(i-1,j)*dv2/vpsi(j))
              a3(k,4)=gam3*((bim1jmh*g4imhjmh*vpsih(j-1)*bijmh+bim1jph*g4imhjph 
     *                *vpsih(j)*bijph)*uzbar*r(i-1,j)+mm**2*b(i,j)*uzbar 
     *                *g2imhj*g3(i-1,j)*dv2/vpsi(j))
 
              a1(k,5)=((f1ijph+f1ijmh)/dv2-f2(i,j))/dt2+gam1*(-(f5ijph+f5ijmh) 
     *                /dv2+f7(i,j)+g1(i,j)+(-bijmh**2*(g4imhjmh+g4iphjmh)*vpsih(j-1) 
     *                -bijph**2*(g4imhjph+g4iphjph)*vpsih(j))*r(i,j)*uzbar 
     *                -mm**2*b(i,j)*uzbar*(g2imhj+g2iphj)*g3(i,j)*dv2/vpsi(j))
cjb a1(k,5)>>> i,j

              a2(k,5)=((f1ijph+f1ijmh)/dv2-f2(i,j))/dt2+gam2*(-(f5ijph+f5ijmh) 
     *                /dv2+f7(i,j)+g1(i,j)+(-bijmh**2*(g4imhjmh+g4iphjmh)*vpsih(j-1) 
     *                -bijph**2*(g4imhjph+g4iphjph)*vpsih(j))*r(i,j)*uzbar 
     *                -mm**2*b(i,j)*uzbar*(g2imhj+g2iphj)*g3(i,j)*dv2/vpsi(j))
              a3(k,5)=((f1ijph+f1ijmh)/dv2-f2(i,j))/dt2+gam3*(-(f5ijph+f5ijmh) 
     *                /dv2+f7(i,j)+g1(i,j)+(-bijmh**2*(g4imhjmh+g4iphjmh)*vpsih(j-1) 
     *                -bijph**2*(g4imhjph+g4iphjph)*vpsih(j))*r(i,j)*uzbar 
     *                -mm**2*b(i,j)*uzbar*(g2imhj+g2iphj)*g3(i,j)*dv2/vpsi(j))
 
              a1(k,6)=gam1*((bip1jmh*g4iphjmh*bijmh*vpsih(j-1)+bip1jph*g4iphjph 
     *                *bijph*vpsih(j))*uzbar*r(i+1,j)+mm**2*b(i,j)*uzbar 
     *                *g2iphj*g3(i+1,j)*dv2/vpsi(j))
cjb a1(k,6)>>> i+1,j

              a2(k,6)=gam2*((bip1jmh*g4iphjmh*bijmh*vpsih(j-1)+bip1jph*g4iphjph 
     *                *bijph*vpsih(j))*uzbar*r(i+1,j)+mm**2*b(i,j)*uzbar 
     *                *g2iphj*g3(i+1,j)*dv2/vpsi(j))
              a3(k,6)=gam3*((bip1jmh*g4iphjmh*bijmh*vpsih(j-1)+bip1jph*g4iphjph 
     *                *bijph*vpsih(j))*uzbar*r(i+1,j)+mm**2*b(i,j)*uzbar 
     *                *g2iphj*g3(i+1,j)*dv2/vpsi(j))
                  
              a1(k,7)=gam1*(-bim1jph*g4imhjph*bijph*vpsih(j)*uzbar*r(i-1,j+1))
cjb a1(k,7)>>> 1-1,j+1

              a2(k,7)=gam2*(-bim1jph*g4imhjph*bijph*vpsih(j)*uzbar*r(i-1,j+1))
              a3(k,7)=gam3*(-bim1jph*g4imhjph*bijph*vpsih(j)*uzbar*r(i-1,j+1))
              a1(k,8)=-f1ijph/(dt2*dv2)+gam1*(f5ijph/dv2+(bijph**2*(g4imhjph 
     *                +g4iphjph)*vpsih(j)*r(i,j+1)*uzbar))
cjb a1(k,8)>>> i,j+1

              a2(k,8)=-f1ijph/(dt2*dv2)+gam2*(f5ijph/dv2+(bijph**2*(g4imhjph 
     *                +g4iphjph)*vpsih(j)*r(i,j+1)*uzbar))
              a3(k,8)=-f1ijph/(dt2*dv2)+gam3*(f5ijph/dv2+(bijph**2*(g4imhjph 
     *                +g4iphjph)*vpsih(j)*r(i,j+1)*uzbar))
              a1(k,9)=gam1*(-bip1jph*g4iphjph*bijph*vpsih(j)*uzbar*r(i+1,j+1))
cjb a1(k,9)>>> i+1,j+1

              a2(k,9)=gam2*(-bip1jph*g4iphjph*bijph*vpsih(j)*uzbar*r(i+1,j+1))
              a3(k,9)=gam3*(-bip1jph*g4iphjph*bijph*vpsih(j)*uzbar*r(i+1,j+1))
c.....b1 array for rhs
              f3ijmh=(f3(i,j)+f3(i,j-1))*.5*vpsih(j-1)
              f3ijph=(f3(i,j)+f3(i,j+1))*.5*vpsih(j)
              denom=1./(dv2)
              b1(k,1)=f3ijmh*denom
              b1(k,3)=f3ijph*denom
              b1(k,2)=-(f3ijmh+f3ijph-f4(i,j)*dv2/vp)*denom
          enddo
      enddo

c.....correct coefficients on boundaries
      sfi1=sign(1.0,fi1)
      sfj1=sign(1.0,fj1)
      sfjrx=sign(1.0,fjrx)
      sfizx=sign(1.0,fizx)
c..... set corners to 0
      k1i=ix-2
      k1j=1+(ix-3)*(jx-2)
      k1=.5*(1+isw)*k1i+.5*(1-isw)*k1j
      fac1=-1.
      if (sfj1 == 1 .and. sfizx == 1) then
          fac1=r(ix-1,2)*b(ix-1,2)/(r(ix,2)*b(ix,2))
      endif
      a1(k1,5)=a1(k1,5)+fac1*a1(k1,3)
      a2(k1,5)=a2(k1,5)+fac1*a2(k1,3)
      a3(k1,5)=a3(k1,5)+fac1*a3(k1,3)
      a1(k1,3)=0.
      a2(k1,3)=0.
      a3(k1,3)=0.
      k2i=1+(ix-2)*(jx-3)
      k2j=jx-2
      k2=.5*(1+isw)*k2i+.5*(1-isw)*k2j
      fac3=-dvout/dvin
      if (sfjrx == 1 .and. sfi1 == 1) then
          fac3=r(2,jx-1)*b(2,jx-1)/(r(1,jx-1)*b(1,jx-1))
      endif
      a1(k2,5)=a1(k2,5)+fac3*a1(k2,7)
      a2(k2,5)=a2(k2,5)+fac3*a2(k2,7)
      a3(k2,5)=a3(k2,5)+fac3*a3(k2,7)
      a1(k2,7)=0.
      a2(k2,7)=0.
      a3(k2,7)=0.
      fac2=-1.
      if (sfj1 == 1 .and. sfi1 == 1) then
          fac2=r(2,2)*b(2,2)/(r(1,2)*b(1,2))
      endif
      a1(1,5)=a1(1,5)+fac2*a1(1,1)
      a2(1,5)=a2(1,5)+fac2*a2(1,1)
      a3(1,5)=a3(1,5)+fac2*a3(1,1)
      a1(1,1)=0.
      a2(1,1)=0.
      a3(1,1)=0.
      fac4=-dvout/dvin
c     if(sfjrx.eq.1.and.sfizx.eq.1)fac4=r(ix-1,jx-1)*b(ix-1,jx-1)/
c    c (r(ix,jx-1)*b(ix,jx-1))
      a1(kxx,5)=a1(kxx,5)+fac4*a1(kxx,9)
      a2(kxx,5)=a2(kxx,5)+fac4*a2(kxx,9)
      a3(kxx,5)=a3(kxx,5)+fac4*a3(kxx,9)
      a1(kxx,9)=0.
      a2(kxx,9)=0.
      a3(kxx,9)=0.
      i=2
      do j=2,jx-1
          if (sfi1 == 1) then
              sfi1=r(2,j)*b(2,j)/(r(1,j)*b(1,j))
          endif
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=2,8,3
              a1(k,m)=a1(k,m)+sfi1*a1(k,m-1)
              a2(k,m)=a2(k,m)+sfi1*a2(k,m-1)
              a3(k,m)=a3(k,m)+sfi1*a3(k,m-1)
          enddo
      enddo
      i=ix-1
      do j=2,jx-1
c     bc(j)=0.
c     if(psi(j).lt.(psih-dpsih))bc(j)=1.
c     if(psi(j).ge.(psih-dpsih).and.psi(j).le.(psih+dpsih))
          bc(j)=.5*(tanh((-psi(j)+psih)/dpsih)+1)
          rbt=r(ix-1,j)*b(ix-1,j)
          rbt1=r(ix,j)*b(ix,j)
          down2=8.*bc(j)*uuzh(ix-1)/du+3.*(1.-bc(j))
          up2=1.-bc(j)
          gm2=up2/down2
          up1=-up2*rbt1*.75+bc(j)*rbt*uuzh(ix-1)/du
          down1=(bc(j)*uuzh(ix-1)/du+up2*3./8.)*rbt1
          gm1=up1/down1
          up=(bc(j)*r(ix-1,j)*b(ix-1,j)*uuzh(ix-1)/du+rbt*(bc(j)-1.))
          down=(bc(j)*r(ix,j)*b(ix,j)*uuzh(ix-1)/du+rbt*(-bc(j)+1.))
          sfizx=up/down
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=2,8,3
              a1(k,m)=a1(k,m)+gm1*a1(k,m+1)
              a1(k,m-1)=a1(k,m-1)+gm2*a1(k,m+1)
              a2(k,m)=a2(k,m)+gm1*a2(k,m+1)
              a2(k,m-1)=a2(k,m-1)+gm2*a2(k,m+1)
              a3(k,m)=a3(k,m)+gm1*a3(k,m+1)
              a3(k,m-1)=a3(k,m-1)+gm2*a3(k,m+1)
          enddo
      enddo

      i=2
      do j=2,jx-1
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=1,7,3
              a1(k,m)=0.
              a2(k,m)=0.
              a3(k,m)=0.
          enddo
      enddo

      i=ix-1
      do j=2,jx-1
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=3,9,3
              a1(k,m)=0.
              a2(k,m)=0.
              a3(k,m)=0.
          enddo
      enddo

      j=2
      do i=2,ix-1
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=4,6
              a1(k,m)=a1(k,m)+sfj1*a1(k,m-3)
              a2(k,m)=a2(k,m)+sfj1*a2(k,m-3)
              a3(k,m)=a3(k,m)+sfj1*a3(k,m-3)
          enddo
          b1(k,2)=b1(k,2)+sfj1*b1(k,1)
      enddo

      j=jx-1
      fac5=sfjrx
      if (sfjrx == -1) then
          fac5=fac5*dvout/dvin
      endif
      do i=2,ix-1
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=4,6
              a1(k,m)=a1(k,m)+fac5*a1(k,m+3)
              a2(k,m)=a2(k,m)+fac5*a2(k,m+3)
              a3(k,m)=a3(k,m)+fac5*a3(k,m+3)
          enddo
          b1(k,2)=b1(k,2)+fac5*b1(k,3)
      enddo

      j=2
      do i=2,ix-1
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=4,6
              a1(k,m-3)=0.
              a2(k,m-3)=0.
              a3(k,m-3)=0.
          enddo
          b1(k,1)=0.
      enddo

      j=jx-1
      do i=2,ix-1
          k1=i-1+(j-2)*(ix-2)
          k2=j-1+(jx-2)*(i-2)
          k=.5*(1+isw)*k1+.5*(1-isw)*k2
          do m=4,6
              a1(k,m+3)=0.
              a3(k,m+3)=0.
              a2(k,m+3)=0.
          enddo
          b1(k,3)=0.
      enddo

      return
      end # amat

      subroutine comat
      implicit none

c.....transforms the elements of the a1(k,m) array into into the
c.....elements of the compressed row matrix abar which will be
c.....operated upon by sgbco and sgbsl.

      Use(ArraySizes)
      Use(Const_1a)
      Use(Matrix)

      integer k1,k2,k3,ml,k,m,mbarj,mbari,mbar,kr,km
      real amaxabar,maxa1
      integer ibrx,jbrx
      character*80 msg

c.....choose bandwith in the smallest grid dimension
      if (isw == 1) then
          k1=ix-2
          k2=1
          k3=ix-2
          ml=ix-1
      else
          k1=1
          k2=jx-2
          k3=jx-2
          ml=jx-1
      endif

      call amax2(abar,lda,kxx,amaxabar,ibrx,jbrx)

      write(msg,'(a,e14.6)')'beginning of comat, amaxabar= ', amaxabar
      call baswline(STDOUT,msg)

      call amax2(a1,kxx,9,maxa1,ibrx,jbrx)

      write(msg,'(a,e14.6)')'begin of comat, maxa1= ', maxa1
      call baswline(STDOUT,msg)

      write(msg,'(a,2i10)') 'ibrx,jbrx=', ibrx,jbrx
      call baswline(STDOUT,msg)
      call baswline(STDOUT,'a1(ibrx-2:ibrx+1,jbrx)')
      write(msg,'(4e14.6)') a1(ibrx-2,jbrx),a1(ibrx-1,jbrx), 
     *           a1(ibrx,jbrx),a1(ibrx+1,jbrx)
      call baswline(STDOUT,msg)
      call baswline(STDOUT,'a1(ibrx,jbrx-2:jbrx+1)')
 
      write(msg,'(4e14.6)') a1(ibrx,jbrx-2),a1(ibrx,jbrx-1),a1(ibrx,jbrx), 
     *            a1(ibrx,jbrx+1)
      call baswline(STDOUT,msg)

      # ? Temporary fix...
      # In the following, subscripts kr and km can have illegal values, so
      # the IF test is required to only assign to abar with valid subscripts.
      do k=1,kxx
          do m=1,9
              mbarj=1+mod(m-1,3)*3+(m-1)/3
              mbari=m
              mbar=0.5*(isw+1)*mbari+0.5*(-isw+1)*mbarj
              kr=ml+3-mod(mbar-1,3)+(2-(mbar-1)/3)*k3
              km=k+((m-1)/3-1)*k1+(mod(m-1,3)-1)*k2
              if (kr >= 1 .and. kr <= lda .and. km >= 1 .and. km <= kxx) then
                  abar(kr,km)=a1(k,m)
              endif
          enddo
      enddo

      return
      end # comat

cjb  abar kr=50 m=9  km=k+49 cols 1,49 not set =lefttri
cjb   abar kr=51 m=6 km=k+48 cols 1,48 not set =lefttri 
cjb   abar kr=52 m=3 km=k+47 cols 1,48 not set  1more 0 than lefttri

cjb   abar kr=98  m=8 km=k+1  col 1 not set =lefttri
cjb   abar kr=100 m=2 km=k-1  col kxx not set =righttrip

cjb   abar kr=146 m=7 km=k-47  48 right cols not set  1 more 0 than righttri
cjb   abar kr=147 m=4 km=k-48  48 right cols of kr147 notset =righttri
cjb   abar kr=148 m=1 km=k-49  49 right cols of kr148  notset =righttri 
