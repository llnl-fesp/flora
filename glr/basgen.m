# $Id: basgen.m,v 1.13 2004/06/30 19:08:45 byers Exp $

      integer function glrgen()
      implicit none

      # Called via ctl command 'generate'

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_3)

      Use(Fstor)

      integer optlist(32),nopt
      real bmx1_,bmx2_,bmx3_,d1trap_,ppas2_,ppas3_

      # Save initial values of input quantities that get modified
      bmx1_=bmx1
      bmx2_=bmx2
      bmx3_=bmx3
      d1trap_=d1trap
      ppas2_=ppas2
      ppas3_=ppas3

      totalcpu=0
      iscreen=max0(1,(nmax/10))
      t=0
      isw=1
      if (ix > jx) isw=-1

      call storage
      call baswline(STDOUT,'finished storage')

      call inputtm
      call baswline(STDOUT,'finished inputtm')
      call pob2poz
      call baswline(STDOUT,'finished pob2poz')
      call constant
      call baswline(STDOUT,'finished constant')
      call tmcon2
      call baswline(STDOUT,'finished tmcon2')
      call equiltm
      call baswline(STDOUT,'finished equiltm')
      call f1to11 
      call baswline(STDOUT,'finished f1to11')
      call amat
      call baswline(STDOUT,'finished amat')
      call comat
      call baswline(STDOUT,'finished comat')


      call initial
      call baswline(STDOUT,'finished initial')
      call mymove(xrol,xro,kxx)
      call baswline(STDOUT,'finished mymovexro')
      call mymove(xiol,xio,kxx)
      call baswline(STDOUT,'finished mymovexio')
      if (kzs == 1) then
          call mymove(xroo,xro,kxx)
          call baswline(STDOUT,'finished mymovexroo')
          call mymove(xioo,xio,kxx)
          call baswline(STDOUT,'finished mymovexioo')
      endif
      call wrabar
      call energy
      call baswline(STDOUT,'finished energy, last one in glrgen')

      # Restore initial values of input quantities for reentry
      bmx1=bmx1_
      bmx2=bmx2_
      bmx3=bmx3_
      d1trap=d1trap_
      ppas2=ppas2_
      ppas3=ppas3_

      glrgen = OK
      return
      end # glrgen

      integer function glrexe(optlist,nopt)
      implicit none

      # Single iteration of main (time) loop.

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_5)
      Use(Fstor)
      Use(Glrexec)
      Use(Matrix)
      Use(TMInput)  

      integer optlist(32),nopt
      real stepcpu,dumbs,second
      external second
      integer l
      character*80 msg

      startime=second(dumbs)
      n=n+1
      if (mod(n,iscreen) == 0) then
          write(msg,'(a,I10,a,e14.6)') 'timestep ',n, ' ; t = ',t
          call baswline(STDOUT,msg)
          write(msg,'(a,e14.6)') '(xro(1)**2+xio(1)**2)**.5       ', \
                    (xro(1)**2+xio(1)**2)**.5
          call baswline(STDOUT,msg)
          write(msg,'(a,e14.6,a,e14.6,a)') ' cumulative loop-CPU = ', \
                    totalcpu,' ; ',stepcpu,' prev. step'
          call baswline(STDOUT,msg)
      endif
      t=t+dt
      time(n)=t
      fac1=-1.0/dt
      fac2=1.0/dt
      do l=0,lmax
          call rightvec  
          call mymove(rhs1b(1:kxx,1),rhs1,kxx)
          call dgbtrs('N',kxx,kbw,kbw,1,abar,lda,ipvt,rhs1b,kxx,info)
          call mymove(rhs1,rhs1b(1:kxx,1),kxx)
          call mymove(xrol,rhs1,kxx)
          call mymove(rhs2b(1:kxx,1),rhs2,kxx)
          call dgbtrs('N',kxx,kbw,kbw,1,abar,lda,ipvt,rhs2b,kxx,info)
          call mymove(rhs2,rhs2b(1:kxx,1),kxx)
          call mymove(xiol,rhs2,kxx)
          fac1=-0.5/dt
          fac2=0.5/dt
      enddo

      if (mod(n,iscreen) == 0) then
          call baswline(STDOUT,'finished loop 90 w dgbtrs, glrexe')
      endif

      call mymove(xioo,xio,kxx)
      call mymove(xio,xiol,kxx)
      call mymove(xroo,xro,kxx)
      call mymove(xro,xrol,kxx)

      xrtime(n)=xro(kplot)
      if (mod(n,ndiag) == 0) call diagno
      if (((mod(n,nen)) == 0) .or. (ltt .ne. 0)) call energy

      endtime=second(dumbs)
      stepcpu=endtime-startime
      totalcpu=totalcpu+stepcpu
      if (mod(n,iscreen) == 0) then
          call baswline(STDOUT,'end of glrexe')
      endif

      if (n < nmax)  then
        glrexe = OK
      else
        glrexe = DONE
      endif
      return
      end # glrexe

      integer function glrfin(optlist,nopt)
      implicit none

      # Finish everything up; called by ctl command 'finish'.

      Use(ArraySizes)
      Use(Const_1)

      integer optlist(32),nopt
      real dumbs,second
      character*80 msg
      external second

      totalcpu=second(dumbs)

      write(msg,'(a,e14.6,a,i10)') 'total CPU for time-loop = ', \
                totalcpu,'    nmax = ',nmax
      call baswline(STDOUT,msg)

      glrfin = OK
      return
      end # glrfin

      subroutine wrabar # Write selected abar elements
      implicit none

      Use(ArraySizes)
      Use(Matrix)

      external dasum
      real dasum
      real amaxabar,genasum
      integer info,ibrx,jbrx
      character*80 msg

      call amax2(abar,lda,kxx,amaxabar,ibrx,jbrx)

      write(msg,'(a,e14.6)') 'before call dgbtrf, amaxabar= ', amaxabar
      call baswline(STDOUT,msg)
      write(msg,'(a,2i10)') 'ibrx,jbrx=', ibrx,jbrx
      call baswline(STDOUT,msg)
      call baswline(STDOUT,'abar(ibrx-2:ibrx+1,jbrx)')
      write(msg,'(4e14.6)') abar(ibrx-2,jbrx),abar(ibrx-1,jbrx), \
                abar(ibrx,jbrx),abar(ibrx+1,jbrx)
      call baswline(STDOUT,msg)
      call baswline(STDOUT,'abar(ibrx,jbrx-2:jbrx+1)')
      write(msg,'(4e14.6)') abar(ibrx,jbrx-2),abar(ibrx,jbrx-1), \
                abar(ibrx,jbrx),abar(ibrx,jbrx+1)
      call baswline(STDOUT,msg)
      genasum=dasum(50,abar(99,1),1) # Why abar(99,1) ?
      write(msg,'(a,e14.6)') 'genasum= dasum(50,abar(99,1),1)= ',genasum
      call baswline(STDOUT,msg)
      call dgbtrf(kxx,kxx,kbw,kbw,abar,lda,ipvt,info)
      call baswline(STDOUT,'finished dgbtrf(kxx,kxx,..)')

      call amax2(abar,lda,kxx,amaxabar,ibrx,jbrx)
      
      write(msg,'(a,e14.6)') 'after call dgbtrf, amaxabar= ', amaxabar
      call baswline(STDOUT,msg)
      write(msg,'(a,2i10)') 'ibrx,jbrx=', ibrx,jbrx
      call baswline(STDOUT,msg)
      call baswline(STDOUT,'abar(ibrx-2:ibrx+1,jbrx)')
      write(msg,'(4e14.6)') abar(ibrx-2,jbrx),abar(ibrx-1,jbrx), \
                abar(ibrx,jbrx),abar(ibrx+1,jbrx)
      call baswline(STDOUT,msg)
      call baswline(STDOUT,'abar(ibrx,jbrx-2:jbrx+1)')
 
      write(msg,'(4e14.6)') abar(ibrx,jbrx-2),abar(ibrx,jbrx-1), \
                abar(ibrx,jbrx),abar(ibrx,jbrx+1)
      call baswline(STDOUT,msg)

      return
      end # wrabar

      subroutine baswline(unit,msg)
      integer unit
      character(len=*) msg
      print *,msg
      return
      end

      subroutine fatal(msg)
      character(len=*) msg

      print *,"FATAL: ",msg
      return
      end
     
