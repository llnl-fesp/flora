# $Id: math.m,v 1.6 2003/08/11 17:54:48 bulmer Exp $
 
      subroutine mymove(a,b,len)
      implicit none

      real a(*),b(*)
      integer len,i

      do i=1,len
          a(i)=b(i)
      enddo

      return
      end # mymove

      subroutine myzero(len,con,a,inc)
      implicit none

      real a(*),con
      integer len,i,inc

      do i=1,len,inc
        a(i)=con
      enddo

      return
      end # myzero

      function quad1(f2,f3,x2,x3)
      implicit none

      real quad1
      real f2,f3,x2,x3,d,a,c

c.....This function calculates the quadrature from 0 to x2,  
c..... assuminng a quadratic integrand of the form a*x**2+c 
cjb   called by  equiltm 
c.....flrd11, improved the accurracy of calculating line integrals
c..... flute1,flute2,flute3,rhoave,xxxave,yyyave,droave by assuming
c..... a quadratic integrand from z=0 to z(1).
 
      d=x2**2-x3**2
      a=(f2-f3)/(3.*d)
      c=(x2**2*f3-x3**2*f2)/d
      quad1=a*x2**3+c*x2

      return
      end # quad1

      subroutine simps(fin,fout,knn,df)
      implicit none

c.....simpsons rule quadratures. knn must be odd   <<<doesnt say over 0toknn

      integer knn,l
      real fin(*),fout,df,se,so
cjb freis      real ssum
cjb      nse=(knn-1)/2
cjb freis     se=ssum(nse,fin(2),2)
cjb        se=0.
cjberr        do l = 2,2+nso-1,2    clearly wrong, nso not even defined yet
cjberr with this nso changed to nse, results for amaxabar
cjberr              from this  f77 doloop version agrees  w f90 SUM version
cjb        do l = 2,2+nse-1,2
cjb          se = se +fin(l)
cjb        enddo
cjb      nso=nse-1
cjb freis     so=ssum(nso,fin(3),2)
cjb        so=0.
cjb        do l = 3,3+nso-1,2
cjb           so = so +fin(l)
cjb         enddo
cjb all above looks wrong  bc of my misinterp of old ssum routine
cjb  what i did in the do loops clearly dont go over entire extent
cjb  really dont need nse,nso defined at all 
cjb  just use   evenpts  2,knn-1,2 oddpts 3,knn-2,2

      se=0.
      do l=2,knn-1,2
         se=se+fin(l)
      enddo
      so=0.
      do l=3,knn-2,2
           so=so+fin(l)
      enddo

      fout=df/(3.*(knn-1))*(fin(1)+fin(knn)+4.*se+2.*so)

cjb internal pts wtd  w   4/3*se +2/3 so   ~(knn-1)/2 [4/3 +2/3]/(knn-1) = 1
cjb  so overall normalization looks ~ok
cjb   if do full counting including endpts and nso=nse-1    assuming allwts=1
cjb      numerator = 1 +1 + [(knn-1)/2]*4  +[(knn-1)/2 -1]*2
cjb                = 2    + (knn-1)*2      +(knn-1)*1 -2
cjb                =        (knn-1)*3     exactly balanced out by denom

cjb  still, clearly the evenpts wtd ~twice as heavy as odd pts
cjb   suspect this is inappropriate even if formal accuracy ok
cjb   ie trapezoidal rule probably better
cjb  in any case seems pretty clear this is a quadratures over entire {1,knn}
cjb   ie not over half the range

      return
      end # simps

      subroutine amax2(a,im,jm,amax,imx,jmx) # Max. of 2D array
      implicit none

      # Find the maximum value of 2D array: amax=max(a(1:im,1:jm)) and the
      # corresponding indices (imx,jmx).

      integer im,jm,imx,jmx,i,j
      real a(im,jm),amax

      amax=a(1,1)
      imx=1
      jmx=2
      do j=1,jm
          do i=1,im
              if (a(i,j) > amax) then
                  amax=a(i,j)
                  imx=i
                  jmx=j
              endif
          enddo
      enddo

      return
      end # amax2

      subroutine spline(x,y,n,bc1,bcn,yp,ypp,yppp,yint) # Cubic splines (from Corsica)
      implicit none
      integer n
      real x(n), y(n), yp(n), ypp(n), yppp(n), yint(n), bc1, bcn
      integer nbeg, nn, i
      real sig, sag, sog, yppi, xn1, rat, bc11, one, two, p, qn, un, bcnn

      nbeg=1
      nn=n
      if (bc1 == 777) then # not-a-knot
          nbeg=2
          sig=(x(2) - x(1))/(x(3) - x(2))
          sag=sig
          ypp(2)=(sig - 1)/(2 + sig)
          yppp(2)=6*((y(3) - y(2))/(x(3) - x(2)) - (y(2) - y(1)) 
     *            /(x(2) - x(1)))/((x(3) - x(2))*(1 + sig)*(2 + sig))
      elseif (bc1 == 555) then # natural, 2nd derivative equals zero
          ypp(1)=0
          yppp(1)=0
      elseif (bc1 == -777) then # third derivative equals zero
          ypp(1)=1
          yppp(1)=0
      elseif (bc1 == -666) then # periodic
          nn=n-1
          xn1=x(n-1) - 2*Pi
          sag=(x(1)-xn1)/(x(2)-xn1)
          sog=(x(n-1) - x(n-2))/(x(n) - x(n-2))
          ypp(1)=(sag - 1)/(sog + 1)
          yppp(1)=6*((y(2) - y(1))/(x(2) - x(1)) - (y(1) - y(n-1)) \
                 /(x(1)-xn1))/(x(2) - xn1)/(sog + 1)
          yint(1)=1.0/(sog + 1)
      else # bc1 is the value of the first derivative 
          if (bc1 == -111) then # extrapolate for 1st derivative
              rat=(x(3) - x(1))/(x(2) - x(1))
              bc11=((y(2) - y(1))*rat - (y(3) - y(1))/rat)/(x(3) - x(2))
          else
              bc11=bc1
          endif
          ypp(1)=-0.5
          yppp(1)=(3.0/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - bc11)
      endif
      one=0
      two=2
      do i=nbeg+1,n-1
          if (i == n-1 .and. bc1 == -666) then
              one=1
              two=2 - sag
          endif
          sig=(x(i) - x(i-1))/(x(i+1) - x(i-1))
          p=sig*ypp(i-1) + two
          ypp(i)=(sig - 1)/p
          yppp(i)=(6*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1)) 
     *            /(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*yppp(i-1))/p
          if (bc1 == -666) then
              yint(i)=(one - sig*yint(i-1))/p
          endif
      enddo

      if (bcn == 777) then # not-a-knot
          sig=1.0/sig - 1
          qn=(ypp(n-2) - 1)*sig - 1
          un=-yppp(n-2)*sig
      elseif(bcn==555) then # natural, 2nd derivative equals zero
          un=0
          qn=0
      elseif (bcn == -777) then # third derivative equals zero
          un=0
          qn=-1.
      elseif (bcn == -666) then # periodic
          ypp(n-1)=yppp(n-1)
          yp(n-1)=yint(n-1)
      else # bc1 is the value of the first derivative
          if (bcn == -111) then # extrapolate for 1st derivative
              rat=(x(n-2) - x(n))/(x(n-1) - x(n))
              bcnn=((y(n-1) - y(n))*rat - (y(n-2) - y(n))/rat)/(x(n-2) - x(n-1))
          else
              bcnn=bcn
          endif
          qn=0.5
          un=3.0/(x(n) - x(n-1))*(bcnn - (y(n) - y(n-1))/(x(n) - x(n-1)))
      endif
      if (bcn .ne. -666) then
          ypp(n)=(un - qn*yppp(n-1))/(qn*ypp(n-1) + 1)
      endif
      do i=nn-1,nbeg,-1
          if (bcn == -666) then
              yp(i)=ypp(i)*yp(i+1) + yint(i)
          endif
          ypp(i)=ypp(i)*ypp(i+1) + yppp(i)
      enddo
      if (nbeg == 2) then
          ypp(1)=ypp(2) + sag*(ypp(2) - ypp(3))
      endif
      if (bcn == -666) then
          yppi=((1 - sog)*ypp(1) + sag*ypp(n-1))/(1 + (1 - sog)*yp(1) \
              + sag*yp(n-1))
          do i=1,n-1
              ypp(i)=ypp(i) - yppi*yp(i)
          enddo
          ypp(n)=ypp(1)
      endif
      do i=1,n-1
          yp(i)=(y(i+1) - y(i))/(x(i+1) - x(i)) - (2*ypp(i) + ypp(i+1)) \
               *(x(i+1) - x(i))/6
          yppp(i)=(ypp(i+1) - ypp(i))/(x(i+1) - x(i))
      enddo
      if (bcn == -666) then
          yppp(n)=yppp(1)
      else
          yppp(n)=yppp(n-1)
      endif
      yp(n)=yp(n-1) + 0.5*(ypp(n) + ypp(n-1))*(x(n) - x(n-1))

      if (yint(1) .ne. -999) then # evaluate integral
          call spline_integral(x,y,yp,ypp,yppp,yint,n)
      endif

      return
      end # spline

      subroutine spline_integral(x,y,yp,ypp,yppp,yint,n) # Spline integral (Corsica)
      implicit none
      integer i, n
      real x(n), y(n), yp(n), ypp(n), yppp(n), yint(n)

      do i=1,n-1
          yint(i+1)=(x(i+1) - x(i))*(y(i) + (x(i+1) - x(i))/2*(yp(i) + 
     *               (x(i+1) - x(i))/3*(ypp(i) + (x(i+1) - x(i))/4*yppp(i))))
      enddo
      yint(1)=0
      do i=2,n
          yint(i)=yint(i-1) + yint(i)
      enddo

      return
      end # spline_integral

      subroutine interp(x,y,n,bc1,bcn,x2,y2,y2p,m) # Spline interpolation (Corsica)
      implicit none
      Use(Splines)
      integer in, n, m, lin, i, l
      real x(*), y(*), x2(*), y2(*), y2p(*), bc1, bcn

      in=n
      if (n < 0) in=-n 
      if (n .ne. 0 .and. (in .ne. mspin .or. m .ne. mspout)) then # Allocate storage
          mspin=in
          mspout=m
          call gchange("Splines",0)
      endif
      
      if (n > 0) then # Calculate increments
          if (x(2) > x(1)) then
              lin=1
          elseif (x(2) < x(1)) then
              lin=-1
          endif
          do i=1,mspout
              do l=1,mspin-1
                  if (lin > 0) then
                      if (x2(i) - x(l) <= 0 .or. 
     *                    (x2(i) - x(l))*(x2(i) - x(l+1)) < 0) then
                          break
                      endif
                  elseif (lin == -1) then
                      if (x2(i) - x(l) >= 0 .or. 
     *                    (x2(i) - x(l))*(x2(i) - x(l+1)) < 0) then
                          break
                      endif
                  endif
              enddo
              lsp(i)=l
          enddo
      endif

      do i=1,mspout
          l=lsp(i)
          dsp(i)=x2(i) - x(l)
      enddo

      if (n .ne. 0) then
          call spline(x,y,in,bc1,bcn,csp1,csp2,csp3,csp4)
      endif

      do i=1,mspout
          l=lsp(i)
          y2(i)=y(l) + dsp(i)*(csp1(l) + dsp(i)/2*(csp2(l) + \
                dsp(i)*csp3(l)/3))
          y2p(i)=csp1(l) + dsp(i)*(csp2(l) + dsp(i)*csp3(l)/2)
      enddo

      return
      end # interp

      subroutine dydx5(x,y,n,yp1,ypn) # Evaluate end-point derivatives
      implicit none
      integer n
      real x(n), y(n), yp1, ypn

      integer i
      real h, c(5)
      data c /-25.0, 48.0, -36.0, 16.0, -3.0/

      h=x(2) - x(1)
      yp1=0
      ypn=0
      do i=1,5
          yp1=yp1 + c(i)*y(i)
          ypn=ypn - c(i)*y(n-i+1)
      enddo
      yp1=yp1/12.0/h
      ypn=ypn/12.0/h

      return
      end # dydx5

