# $Id: initial.m,v 1.6 2003/06/04 16:56:44 bulmer Exp $

      subroutine initial
      implicit none
 
c.....set up initial displacement vectors, xro and xio
c.....test case 1, cos(kz) in z, flat in psi
c.....set up 1/4/82 by r. freis

      Use(ArraySizes)
      Use(GConst)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_2)
      Use(Const_4)
      Use(Const_8)
      Use(Fstor)
      Use(Matrix)

      real ranf
      integer i,j,kzsp,k,k1,k2
      real rbf,bnorm,r1,r2,r3,r4,cosx

      rbf=1.
      bnorm=1.
c.... establish the random no. seed

cjb takeout for test      ra=ranset(0)
cjb  try putting back in as rnset  as seen in icepic
cjb  using rnset allows compile     but rnset then undefined symbol incode.err
cjbout need to repl with some newer seed method      ra=rnset(0)

      do j=2,jx-1
c..... kzs=1 for stiff perturbation (initially) in psi. Otherwise
c.....  initial perturbation is random in psi.
c..... If kzs=1, initial perturbation has 1/(bvac(i)) z dependence.
          if (kzs == 1) then
              if (psi(j) <= psi0) then
                  r1=1.
                  r2=1.
                  r3=1.
                  r4=1.
              else
                  r1=1.
                  r2=(psi(jx)-psi(j))/(psi(jx)-psi0)
                  r4=(psi(jx)-psi(j))/(psi(jx)-psi0)
                  r3=1.
              endif
          else
              r1=ranf(b1)
              r2=ranf(b1)
              r3=ranf(b1)
              r4=ranf(b1)
          endif
          do i=2,ix-1
              if (psih > psi(j)) then
                  rbf=1./(r(i,j)*b(i,j))
                  kzsp=0
              else
                  r1=ranf(b1)
                  r2=ranf(b1)
                  r3=ranf(b1)
                  r4=ranf(b1)
                  kzsp=1
              endif
              k1=i-1+(j-2)*(ix-2)
              k2=j-1+(jx-2)*(i-2)
              k=.5*(1+isw)*k1+.5*(1-isw)*k2
              if (kzs == 1) then
                  rbf=1.
                  bnorm=bvac(2)/bvac(i)
              endif
              cosx=cos(.5*Pi*(z(i)*kzsp)/zedge)*bnorm
              xro(k)=ex0*cosx*rbf*(r1+r2-1.)+ex1*cos(.5*Pi*(z(i))/zedge)
              xio(k)=ex0*cosx*rbf*(r3+r4-1.)+ex1*cos(.5*Pi*z(i)/zedge)
c     xio(k)=cos(theta0)*xro(k)


cjb coding above, at least for kzs=1, isw=-1 
cjb give const out to psi0 then linear dropoff to wall    and no zdep
cjb  for my tests with psi0 > psimax   this gave strong transient
cjb  and for psi0= .5psimax excites that  sudden jump in slope
cjb  that appears to be an allowed mode of the system at least for exb
cjb   idea is that d2dpsi2 is zero in the linear regions then a jump at onej

cjb  when imposing linear falloff thruout 
cjb  we both avoid strong transient  and  dont see that mode with jumpinslope
cjb  this linear falloff allpsi is close to J1(r) for ~const rho besselmodel
cjb  not sure we want it for other, flat better for rigidrot m=1

cjb      if(doflat) then
              if (doflat .and. doxrlin) then
                  xro(k)=ex0*(psi(jx)-psi(j))/psi(jx)
                  xio(k)=xro(k)
              elseif (doflat .and. doxrtanh) then
                  xro(k)=-ex0*tanh(.1*(j-jx))
                  xio(k)=xro(k)
              endif

cjb both the doxrlin and the doxrtanh branches avoid use of psi0
cjb  doxrlin=true could be sensible for m=2 expecting J1(r)
cjb  doxrtanh=true for m=1  flat low psi falls off smoothly to 0 j=jx
cjb  -tanh(-4.9) = .999889  value at jx-j = 49  j=1
cjb  -tanh(-2.5) = .986614  value at jx-j = 25
          enddo
      enddo

cjbjul31 2002 write out the initial xro vs psi

      if (kzs == 1) return
      do j=2,jx-1
          r1=ranf(b1)
          r2=ranf(b1)
          r3=ranf(b1)
          r4=ranf(b1)
          do i=2,ix-1
              if (psih > psi(j)) then
                  rbf=1./(r(i,j)*b(i,j))
                  kzsp=0
              else
                  r1=ranf(b1)
                  r2=ranf(b1)
                  r3=ranf(b1)
                  r4=ranf(b1)
                  kzsp=1
              endif
              k1=i-1+(j-2)*(ix-2)
              k2=j-1+(jx-2)*(i-2)
              k=.5*(1+isw)*k1+.5*(1-isw)*k2
              cosx=cos(.5*Pi*(z(i)*kzsp)/zedge)
              xroo(k)=ex0*cosx*rbf*(r1+r2-1.)+ex1*cos(.5*Pi*(z(i))/zedge)
              xioo(k)=ex0*cosx*rbf*(r3+r4-1.)+ex1*cos(.5*Pi*z(i)/zedge)
c     xioo(k)=cos(theta0)*xroo(k)
          enddo
      enddo

      return
      end # initial
