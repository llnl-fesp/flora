# $Id: rightvec.m,v 1.5 2003/06/04 17:00:30 bulmer Exp $

      subroutine rightvec
      implicit none

c..... calculates right hand side vector for both equations,
c..... rhs1(k)=2*a2*xr(n)-a2*xr(n-1)+b1*(xi(l)-xi(n-1)) , and
c..... rhs2(k)=2*a3*xi(n)-a2*xi(n-1)+b1*(xr(l)-xr(n-1)) .

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
      Use(Fstor)
      Use(Matrix)

      integer m,mdel,koff1,koff2,koff,k,mbar,m1

      call myzero(kxx,0.,rhs1,1)
      call myzero(kxx,0.,rhs2,1)
cjb  what is  sscal?  some system call,  not in vdf  not in undefined symbols
cjb  just some vector 0 init?
cjb   sscal is in code.map  part of libblas.a

cjb  where kxx set?:  declrd  not intlzd in vdf

      do m=1,9
          mdel=(m-1)/3
          m1=m-1
          koff1=-mdel*5+m+(mdel-1)*ix
          koff2=(m-((mdel+1)*3-1))*jx+1-2*m1+7*mdel
          koff=.5*(1+isw)*koff1+.5*(1-isw)*koff2
          do k=1,kxx
              rhs1(k)=rhs1(k)+2.*a3(k,m)*xro(k+koff)-a2(k,m)*xroo(k+koff)
              rhs2(k)=rhs2(k)+2.*a3(k,m)*xio(k+koff)-a2(k,m)*xioo(k+koff)
          enddo
          if (m == 2 .or. m == 5 .or. m == 8) then
              do k=1,kxx
                  mbar=m-1-(m/4)*2
                  rhs1(k)=rhs1(k)+fac1*b1(k,mbar)*(xiol(k+koff)-xioo(k+koff))
                  rhs2(k)=rhs2(k)+fac2*b1(k,mbar)*(xrol(k+koff)-xroo(k+koff))
              enddo
          endif
      enddo

      return
      end # rightvec
