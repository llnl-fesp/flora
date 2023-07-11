# $i: equil.m,v 1.8 2003/06/04 16:55:38 bulmer Exp $

      subroutine tmcon2
      implicit none

c.....calculates equilibruim quantities for tandem mirror
 
      Use(ArraySizes)
      Use(Coils)
      Use(Const_2)
      Use(Const_3)
      Use(Const_4)
      Use(Const_7)
      Use(Const_8)
      Use(GFiducials)
      Use(Fstor)
      Use(TMInput)

      integer i,j,isubpx
      real rp3,alsit,con1,con2,con3,con4,con5,con6,don1,don2
      real don3,p30,wp2,wp3,psi0r,g2p,g2f,rtemp,rpxi2,ter1,ter2
      real bstr2,pfloat,pcen1,altemp,betemp,denfac,ep1,fp1,gp1
      real ep2,fp2,gp2,ep3,fp3,gp3,factr,factr1,et1,ft1,gt1
      real denfac0,et0,ft0,gt0,eon1,eon2,eon3,at00,bt00,ct00
      real bstr,faccb,denfac1

cjb put alsiks and bvoutks into vdf;  check how alsi used elsewhere
      real alsitks,con1ks,con2ks,con3ks,don1ks,don2ks,don3ks
      real rtempks,rpxi2ks,ter1ks,ter2ks
      real bstr2ks,pfloatks
      real eon1ks,eon2ks,eon3ks,p2filt
      character*80 msg

c.....calculate heaveside step functions h12, hzt0, h34
      call myzero(ix,0.,h12,1)
      call myzero(ix,0.,hzt0,1)
      call myzero(ix,0.,h34,1)
      call myzero(ix,0.,htrans,1)
      call myzero(ix,0.,hzp0,1)
      call myzero(ix,0.,hzp1,1)
      call myzero(ix,0.,hzp2,1)
      call myzero(ix,0.,hzp3,1)
      call myzero(ix,0.,hzt1,1)
      call myzero(ix,0.,hzt2,1)
      call myzero(ix,0.,hzt3,1)

      do i=1,ix
          if (z(i) <= z1c) hzt0(i)=1.
          if (z(i) >= ztran) htrans(i)=1.
          if ((z(i) >= z1c) .and. (z(i) <= z3) .and. (bvac(i) <= bv3)) 
     *        h34(i)=1.
      enddo
 
      psim=psi0*(1.-p2wide)
      psime=psi0e*(1.-p2ewide)
c.....calculate heaveside step functions hp0(j), hp3(j)
      call myzero(jx,0.,hp3,1)
      call myzero(jx,0.,hp0,1)
      call myzero(jx,0.,hpm,1)
      call myzero(jx,0.,hpme,1)

 
      do j=1,jx
          if (psi(j) <= psi0) hp0(j)=1.
          if (psi(j) <= psi3) hp3(j)=1.
          if (psi(j) <= psim) hpm(j)=1.
          if (psi(j) <= psime) hpme(j)=1.
      enddo
 
      rp3=bv3/(+bmn1)
      rpx=bvx3/(bmn2)
      bm1=rp1*( bmn2)
      alsit=((1.-1./rp1**2)/(1.-1./rpx**2)*(rp1/rpx))**2
      alsi=alsit+.5*(1.-alsit)*pfudge

      write(msg,'(a,3e14.6)') 'rpx,bm1,rp1,alsit,alsi',rpx,bm1,rp1
      call baswline(STDOUT,msg)
      write(msg,'(22x,2e14.6)') alsit,alsi
      call baswline(STDOUT,msg)

cjb mimic sloshing rpx,bm1 for ks
cjb instead of using bvx3 for the wider comp of slsks   read in rpxks
      if (rpxks == 0) rpxks=bvx3/(bvac(ix))
      bm1ks=rp1ks*( bvac(ix))
      alsitks=((1.-1./rp1ks**2)/(1.-1./rpxks**2)*(rp1ks/rpxks))**2
      alsiks=alsitks+.5*(1.-alsitks)*pfudgeks
      write(msg,'(a,3e14.6)') 'rpxks,bm1ks,rp1ks,alsitks,alsiks', \
                rpxks,bm1ks,rp1ks
      call baswline(STDOUT,msg)
      write(msg,'(32x,2e14.6)') alsitks,alsiks
      call baswline(STDOUT,msg)
 
c......calculate heaveside functions for passing components and
c...... trapped components in region 1 and 2, and region 3 and 4
      do i=1,ix
          if ((z(i) < z1c)) hzp0(i)=1
          if ((z(i) >= z1c) .and. (z(i) < z1min)) hzp1(i)=1
          if ((z(i) >= z1min) .and. (z(i) < z2c)) hzp2(i)=1
          if ((z(i) >= z2c) .and. (z(i) < z2min)) hzp3(i)=1
          if ((z(i) > z1c) .and. (z(i) <= z2c) .and. (bvac(i) 
     *           <= bvx2)) hzt1(i)=1
          if ((z(i) > z2ct) .and. (z(i) <= z3c) .and. (bvac(i)
     *           <= bvx3)) hzt2(i)=1
          if ((z(i) > z2ct) .and. (z(i) <= z3c) .and. (bvac(i) 
     *           <= bm1)) hzt3(i)=1

cjb       mimic last 2 for ks region
          bvoutks=rpxks*bvac(ix)

cjb       if((z(i).gt.z3c).and.(z(i).le.zmax).and.(bvac(i).le.bvx3))
          if ((z(i) > z3c) .and. (z(i) <= zmax) .and. (bvac(i) 
     *           <= bvoutks)) hzt2ks(i)=1
          if ((z(i) > z3c) .and. (z(i) <= zmax) .and. (bvac(i) 
     *           <= bm1ks)) hzt3ks(i)=1
      enddo
c..... calculate pressure coeff. including electron ring (pr)
cjb con1,2,3 for sloshing in plug
      con1=1./((1.-1./rpx**2)**2*bvx3**4)
      con2=con1*(-2*bvx3**2)
      con3=-.5*con2*bvx3**2
cjb sametype con1,2,3 for ks region

cjb      con1ks= 1./((1.-1./rpxks**2)**2*bvx3**4)
cjb      con2ks=con1ks*(-2*bvx3**2)
cjb      con3ks=-.5*con2ks*bvx3**2
      con1ks=1./((1.-1./rpxks**2)**2*bvoutks**4)
      con2ks=con1ks*(-2*bvoutks**2)
      con3ks=-.5*con2ks*bvoutks**2



cjb con456 for pring
      con4=pring/((1.-1./rp3**2)**2*bv3**4)
      con5=con4*(-2*bv3**2)
      con6=-.5*con5*bv3**2

cjb don1,2,3 for sloshing in plug
      don1=1./((1.-1./rp1**2)**2*bm1**4)
      don2=don1*(-2*bm1**2)
      don3=-.5*don2*bm1**2
cjb sametype don1,2,3 for ks region
      don1ks=1./((1.-1./rp1ks**2)**2*bm1ks**4)
      don2ks=don1ks*(-2*bm1ks**2)
      don3ks=-.5*don2ks*bm1ks**2
cjb       write(*,*)'con1,con2,con3',con1,con2,con3
cjb       write(*,*)'don1,don2,don3',don1,don2,don3

cjb       write(*,*)'con1ks,con2ks,con3ks',con1ks,con2ks,con3ks
cjb       write(*,*)'don1ks,don2ks,don3ks',don1ks,don2ks,don3ks

c..... calculate constants for p3 pressure component
cjb  p30 not in vdf 
      data p30/.01/

      wp2=p2wide*.5*psi0
      wp2e=p2ewide*.5*psi0e
      wp3=wp2/psim
      psi0r=psi0/psim
      p3b=p30
      p3c=(.5*(1.-(tanh(2.))**2)/wp3-p3b*p2wide*(2.-p2wide))*.5/p2wide
      p3d=-(p3b+2.*p3c*psi0r)/(3.*psi0r**2)
      p3a=-(p3b*psi0r+p3c*psi0r**2+p3d*psi0r**3)
 
c.....calculate constants for pe1
      g2p=-.5*(1.-(tanh(-2.))**2)/wp2e
      g2f=.5*(1.-tanh(-2.))
      ae1=pe10
      ce1=psime*g2p-g2f+ae1
      be1=2.*(g2f-ae1)-psime*g2p
      be1=psislp
      de1=psime*g2p+be1-2.*(g2f-ae1)
      ce1=-g2p*psime-2.*be1+3.*(g2f-ae1)

c......readjust peak plug pressure for sloshing
      rtemp=(bm1/bvx3)**2
      rpxi2=1./rpx**2
      ter1=(rtemp-rpxi2)**2
      ter2=(1.-rpxi2)**2
      bstr2=bvx3**2*(-ter1+alsi*rtemp*ter2)/(-ter1+alsi*ter2)
      pfloat=((con1-don1*alsi)*bstr2**2+(con2-don2*alsi)*bstr2 \
            +con3-don3*alsi)*(1.+p3a*dip)

      write(msg,'(a,3e14.6)') 'rtemp,rpxi2,ter1,ter2,bstr2,pfloat', \
                rtemp,rpxi2,ter1
      call baswline(STDOUT,msg)
      write(msg,'(34x,3e14.6)') ter2,bstr2,pfloat
      call baswline(STDOUT,msg)
      if (psloshin == 0) then
          pfloat=(1.+p3a*dip)
          write(msg,'(a,e14.6)') 'after psloshin=0test, pfloat=',pfloat
          call baswline(STDOUT,msg)
      endif

cjb      rtempks=(bm1ks/bvx3)**2
      rtempks=(bm1ks/bvoutks)**2
      rpxi2ks=1./rpxks**2
      ter1ks=(rtempks-rpxi2ks)**2
      ter2ks=(1.-rpxi2ks)**2
cjb      bstr2ks=bvx3**2*(-ter1ks+alsiks*rtempks*ter2ks)/
cjb     .        (-ter1ks+alsiks*ter2ks)
      bstr2ks=bvoutks**2*(-ter1ks+alsiks*rtempks*ter2ks)/ \
              (-ter1ks+alsiks*ter2ks)

      pfloatks=((con1ks-don1ks*alsiks)*bstr2ks**2 \
              +(con2ks-don2ks*alsiks)*bstr2ks+con3ks-don3ks*alsiks)

      write(msg,'(a,2e14.6)') 'rtempks,rpxi2ks,ter1ks,ter2ks,bstr2ks,pfloatks', \
                    rtempks,rpxi2ks
      call baswline(STDOUT,msg)
      write(msg,'(46x,2e14.6)') ter1ks,ter2ks
      call baswline(STDOUT,msg)
      write(msg,'(46x,2e14.6)') bstr2ks,pfloatks
      call baswline(STDOUT,msg)
      if (poshinks == 0)then
          pfloatks=1.
          write(msg,'(a,e14.6)') 'after poshinks=0test, pfloatks=',pfloatks
          call baswline(STDOUT,msg)
      endif

cjb  now all pslosh,psloshe,nslosh adj with /pfloat
      pslosh=psloshin/pfloat
      psloshe=psloshen/pfloat
      nslosh=nsloshin/pfloat

      psloshks=poshinks/pfloatks
      nsloshks=noshinks/pfloatks


cjb also adjusts ppas123,dpas1   with /pfloat
cjb also adjusts d1trap with /pfloat  but not p1trap
cjb perhaps idea is:  if accept need for adj on sloshing species
cjb   then also need it for any species in contact with the sloshing
cjb   so passing into ppas3 in direct contact w sloshing
cjb   and then via ppas2 passing contact w choke d1trap?  but p1trap not adj
cjb   note pcent not adusted

cjb  if run w ncoil=2 i think d1trap,p1trap not used  and no passing in plug
cjb   the code will still adj ppas1 with pfloat
cjb  if run w ncoil=2 so no d1trap, and w no passing only pcenter
cjb   then these additional adjustments not made
      ppas1=ppas1/pfloat
      ppas2=ppas2/pfloat
      ppas3=ppas3/pfloat
      d1trap=d1trap/pfloat
      dpas1=dpas1/pfloat

      psls1=pslosh+psloshe
      pcen1=pcenter+pcentee
 
c....calculate passing pperp constants
      cp0=ppas1
      bp1=(ppas1-ppas2)*bmax**2*2./(bmax**2-bmn1**2)**2
      ap1=-bp1/(2.*bmax**2)
      cp1=ppas1-bp1*bmax**2-ap1*bmax**4
      altemp=bp1*(1.-(bmn1/bmax)**2)/(2.*(bmn1**2))
      bp2=(ppas3-ppas2-altemp*(bvx2**4-bmn1**4))/(bvx2**2-bmn1**2 \
         -.5*(bvx2**4-bmn1**4)/bmn1**2)
      ap2=altemp-bp2*.5/bmn1**2
      cp2=ppas2-ap2*bmn1**4-bp2*bmn1**2
      betemp=4.*ap2*bvx2**2+2.*bp2
      ap3=(-ppas3+.5*betemp*(bvx2**2-bmn2**2))/(bvx2**2-bmn2**2)**2
      bp3=.5*(betemp-4.*bvx2**2*ap3)
      cp3=-ap3*(bmn2**4)-.5*bmn2**2*(betemp-4.*bvx2**2*ap3)
c..... calculate passing density constants in region 1 2 3
      denfac=0.
      if (ppas1 > 0) denfac=dpas1/ppas1
      ep1=ap1*denfac
      fp1=bp1*denfac
      gp1=cp1*denfac
      ep2=ap2*denfac
      fp2=bp2*denfac
      gp2=cp2*denfac
      ep3=ap3*denfac
      fp3=bp3*denfac
      gp3=cp3*denfac
c.....calculate ptrapped constants in region 1 and 2 and 0
      factr=p1trap/(1.-(bmn1/bvx2)**2)**2
      at1=factr/bvx2**4
      bt1=-2.*factr/bvx2**2
      ct1=factr
      factr1=pcen1/(1.-(bv0/bmax)**2)**2
      at0=factr1/bmax**4
      bt0=-2.*factr1/bmax**2
      ct0=factr1
      at2=con1*psls1
      bt2=con2*psls1

cjb at2,bt2,at3,bt3 not used
cjb  ct2 used in dt2; ct3 used in dt3
cjb  add  ct2ks, ct3ks  mimic ct2,ct3
      ct2=con3*psls1
      ct2ks=con3ks*psloshks

      at3=don1*alsi*psls1
      bt3=don2*alsi*psls1

      ct3=don3*alsi*psls1
      ct3ks=don3ks*alsiks*psloshks



c.....calculate trapped density constants in region 1 2 and 0
      denfac1=0.
      if (p1trap > 0) denfac1=d1trap/p1trap
      et1=at1*denfac1
      ft1=bt1*denfac1
      gt1=ct1*denfac1
      denfac0=0.
      if (pcen1 > 0) denfac0=ncenter/pcen1
      et0=at0*denfac0
      ft0=bt0*denfac0
      gt0=ct0*denfac0
      do i=1,ix
          eon1=con1-don1*alsi*hzt3(i)
          eon2=con2-don2*alsi*hzt3(i)
          eon3=con3-don3*alsi*hzt3(i)      
cjb      if(hzt2(i)=1.)write(*,*)'i,eon1,eon2,eon3',i,eon1,eon2,eon3

          eon1ks=con1ks-don1ks*alsiks*hzt3ks(i)
          eon2ks=con2ks-don2ks*alsiks*hzt3ks(i)
          eon3ks=con3ks-don3ks*alsiks*hzt3ks(i)
cjb       if(hzt2ks(i)=1.)
cjb      . write(*,*)'i,eon1ks,eon2ks,eon3ks',i,eon1ks,eon2ks,eon3ks

          abf(i)=con4*h34(i)
          bbf(i)=con5*h34(i)
          cbf(i)=con6*h34(i)
 
          abp(i)=eon1*psls1*hzt2(i)+at0*hzt0(i)+ap1*hzp1(i)+ap2*hzp2(i) 
     *           +ap3*hzp3(i)+at1*hzt1(i)+eon1ks*psloshks*hzt2ks(i)
          bbp(i)=eon2*psls1*hzt2(i)+bt0*hzt0(i)+bp1*hzp1(i)+bp2*hzp2(i) 
     *           +bp3*hzp3(i)+bt1*hzt1(i)+eon2ks*psloshks*hzt2ks(i)
          cbp(i)=eon3*psls1*hzt2(i)+ct0*hzt0(i)+cp1*hzp1(i)+cp2*hzp2(i) 
     *           +cp3*hzp3(i)+cp0*hzp0(i)+ct1*hzt1(i)+eon3ks*psloshks*hzt2ks(i)

cjb      if(dosubp. and. i.ge.50) then
          if (dosubp .and. i > ibvmax(2)) then
cjb dont do p2filt in plug          if (dosubp .and. i > ibvmax(1)) then
cjb this ifloop is only for the p2filt,  really dont need this for plug 
              p2filt=.5*(1.-tanh(alfsubp*(i-isubp)))
              abp(i)=p2filt*abp(i)
              bbp(i)=p2filt*bbp(i)
              cbp(i)=p2filt*cbp(i)
              jdfz(i)=p2filt*jdfz(i)
              jd2fz(i)=p2filt*jd2fz(i)
              jperpz(i)=p2filt*jperpz(i)
              jparz(i)=p2filt*jparz(i)

              isubpx=(ix+isubp)/2+1
              if (i >= isubpx) then
                  abp(i)=0.
                  bbp(i)=0.
                  cbp(i)=0.
                  jdfz(i)=0.
                  jd2fz(i)=0.
                  jperpz(i)=0.
                  jparz(i)=0.
              endif
          endif


          ebp(i)=eon1*nslosh*hzt2(i)+et0*hzt0(i)+ep1*hzp1(i)+ep2*hzp2(i) 
     *           +ep3*hzp3(i)+et1*hzt1(i)+eon1ks*nsloshks*hzt2ks(i)
          fbp(i)=eon2*nslosh*hzt2(i)+ft0*hzt0(i)+fp1*hzp1(i)+fp2*hzp2(i) 
     *           +fp3*hzp3(i)+ft1*hzt1(i)+eon2ks*nsloshks*hzt2ks(i)
          gbp(i)=eon3*nslosh*hzt2(i)+gt0*hzt0(i)+gp1*hzp1(i)+gp2*hzp2(i) 
     *           +gp3*hzp3(i)+gt1*hzt1(i)+dpas1*hzt0(i)+eon3ks*nsloshks*hzt2ks(i)

          at00=at0
          bt00=bt0
          ct00=ct0
          if (pcen1 .ne. 0) then
              at00=at0*pcenter/pcen1
              bt00=bt0*pcenter/pcen1
              ct00=ct0*pcenter/pcen1
          endif
          abq(i)=eon1*pslosh*hzt2(i)+at00*hzt0(i)+(ap1*hzp1(i)+ap2*hzp2(i) 
     *           +ap3*hzp3(i)+at1*hzt1(i))*.5+eon1ks*psloshks*hzt2ks(i)  
          bbq(i)=eon2*pslosh*hzt2(i)+bt00*hzt0(i)+(bp1*hzp1(i)+bp2*hzp2(i) 
     *           +bp3*hzp3(i)+bt1*hzt1(i))*.5+eon2ks*psloshks*hzt2ks(i)  
          cbq(i)=eon3*pslosh*hzt2(i)+ct00*hzt0(i)+(cp1*hzp1(i)+cp2*hzp2(i) 
     *           +cp3*hzp3(i)+cp0*hzp0(i)+ct1*hzt1(i))*.5 
     *           +eon3ks*psloshks*hzt2ks(i)

          if (dosubp .and. i >= 50) then
              p2filt=.5*(1.-tanh(alfsubp*(i-isubp)))
              abq(i)=p2filt*abq(i)
              bbq(i)=p2filt*bbq(i)
              cbq(i)=p2filt*cbq(i)
              if (i >= isubpx) then
                  abq(i)=0.
                  bbq(i)=0.
                  cbq(i)=0.
              endif
          endif

cjb apr24 2002   this .5 factor on ks sloshing  is inconsistent 
cjb   the abq.. coeffs are for ions-only flrterms
cjb   eon1*pslosh*hzt2  cpt eon1*psls1*hzt2 in abp -- ionsonly
cjb   via at00,  at0*hzt0 in abp reduced by pcenter/pcen1--ionsonly
cjb    (ap1,2,3  +at1)  .5 factor   for pas1,2,3 and trap   ionsonly te=ti
cjb  to be fully parallel to plugsloshing should treat exactly same
cjb    take psloshks out of the .5* term
cjb     --- get ks equiv of  psls1=pslosh+psloshe
cjb     perhaps same if use psloshe=0
cjb  this discrepancy no effect if sf68=0
      enddo 

c..... constants for ering pperp(b) long-thin
      if (long .ne. 0) then
          bstr=bvx2*fring
          faccb=1./(bvx2**2-bstr**2)**2
          do i=1,ix
              if (bvac(i) >= bstr) then
                  abf(i)=-pring*faccb*h34(i)
                  bbf(i)=2.*pring*faccb*h34(i)*bstr**2
                  cbf(i)=pring*bvx2**2*(bvx2**2-2.*bstr**2)*faccb*h34(i)
              else
                  cbf(i)=pring*h34(i)
                  bbf(i)=0.
                  abf(i)=0.
              endif
          enddo
      endif

      return
      end # tmcon2

      subroutine equiltm
cjb feb 6 2004  start newplug model.
cjb  prepare code for using mma orbit code output 6barrays
cjb   using same concept as presently doing for ks:
cjb   mma orbit code picks desirable f(th) imagined to be 'injected'
cjb   at the plugmid;  orbit code then generates pperp[b/bmidpl]
cjb    and mma also generates 5 other f[b/bmidpl]
cjb  That operation is independent of betavalue chosen in flora.
cjb  The 6 f[b] arrays are read into flora at the beginning.
cjb  The beta value comes in when solving for btot in flora
cjb   in an iterative process, 
cjb   which also produces jperpz, normalized to max=unity,
cjb   and pperp[z]then comes from poshinks*jperpz


cjb  start preparation for newplug code by stripping out older cjb cmts

cjb  still a mess in code for errprl
cjb  still needs delb<epsb done better or rejected
cjb  still a problem i think on rterm,droternm, droave... diagnostics
cjb  still needs freis bsq1 v bsq2  in the r,rzz 2nd half of equiltm
cjb        modified the way i did for first half




      implicit none

c.....equilibrium for tandem mirror electron ring plug

      Use(ArraySizes)
      Use(GConst)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_2)
      Use(Const_3)
      Use(Const_4)
      Use(Const_6)
      Use(Const_7)
      Use(Const_8)
      Use(GFiducials)
      Use(Fstor)

      integer j,i,k,is,is1,ia,isubpx
      real rzzr,xxxabs
      real epsb,wp2,psibar,psir,psire,pe1,pe2,u1,u2,u3,at,bt
      real ct,bv,ordera2,bsq3,t1,radic2,round,radic,bsq1,bsq2,bsq
      real arg1,arg2,dp1,dp2,dp3,dt0,dt1,dt2,dt3,ppt,dter1,dter2
      real b2,b4,dter,ppart,ppert,pperte,fac3,dp2db,dp2dbe
      real dp3db,pharg,dphidpsi,dpdpsi,omegci,unit,omegstr,omeggb
      real ppertr,delb,abar,bbar,dter4,dbar,factor,sum1,sum2,sum3
      real sum4,sum5,dps,psib,psik,psikr,psikre,p1ek,p2ek,p3k
      real bk,capf,dcapi1,dcapi2,dcapi3,dcapi4,dcapi5,capi1
      real capi2,capi3,capi4,capi5,vv1,vv2,eterm,rterm,ringj,ans1
      real ans2,ans3,ans4,ans5,omegwkb,trad
      real pans1r,pans2r,pansjb
      real capfnew,dcap5new,cap5new,sum5new
      real         dcap2new,cap2new,sum2new
      real         dcap3new,cap3new,sum3new
      real         dcap1new,cap1new,sum1new
      real bsqjb,bkjb,roldpg,rnewks,rzzoldpg,rzznewks
      real ppertold,ppartold,dp3dbold
      integer iia
      real sumpostint
      real dt2ks,dt3ks
cjb add c2tks,c23ks to vdf  to mimic ct2,ct3
      real quad1
cjb these data vars not in vdf  so must be local 
      real ss

cjb add local vars for pressure filter
      real p2filt
      integer l
      character*100 msg

cjb jan29 2004 
cjb changed *80 to *100 to prevent some runtime 'toolongrecord'errors
cjb actually equil.mwr50 got around that error by adding addl writes
cjb  so this change is for future protection

      character*80 msg4(4)

      ss= (4.*Pi)**.5

      data epsb/1.e-4/
c..... loop 10, calculate p2 p1 and b and dp2/dpsi and dp1/dpsi
      wp2=p2wide*.5*psi0
      psibar=-be1*.5/ce1
      psistr=-ce1/(3.*de1)
      psibar=psistr+sqrt(psistr**2-be1/(3.*de1))
      p1max=ae1+be1*psibar+ce1*psibar**2+de1*psibar**3
 
      do 10 j=1,jx
          psir=psi(j)/psim
          psire=psi(j)/psime
          pe1=(ae1+be1*psire+ce1*psire**2+de1*psire**3)*hpme(j)
          pe2=.5*(1.-tanh((psi(j)-psi0e)/wp2e))*(1.-hpme(j))
          p1(j)=(pe1+pe2)/p1max
          dp1dpsi(j)=(de1*3.*psire**2 
     *               +be1+2.*ce1*psire)/(psime*p1max)*hpme(j) 
     *               -.5*(1.-(tanh((psi(j)-psi0e)/wp2e))**2) 
     *               /(wp2e*p1max)*(1.-hpme(j))

          p2t(j)=(1.-tanh((psi(j)-psi0)/wp2))*.5
      
          dp2dpst(j)=-.5*(1.-(tanh((psi(j)-psi0)/wp2))**2)/wp2
cjb replace tanh psi model  with linear decrease, const dpdpsi

cjb  p2t(j)=(1.-tanh((psi(j)-psi0)/wp2))*.5
cjb  psi=0  .5(1-th(-psi0/(.5psi0*p2wide)))
cjb         .5(1+th(2/p2wide))   for p2wide<<1  th->1
cjb  psi=0 ->.5(1+1) = 1

cjb  psi=psi0  .5(1-th(0)) = .5

cjb  psi=2*psio  .5(1-th(psi0/(.5*psi0*p2wide)))
cjb              .5(1-th(2/p2wide))  ->.5(1-1)=0  for p2wide<<1

cjb replace this w linear  psi=0 p2t=1, psi=psi0 p2t=.5  psi=2psi0 p2t=0
cjb  i think this rqrs psi0rel=.5
cjb      if(doflat) then
cjb        p2t(j)= 1.-.5*psi(j)/psi0
cjb        dp2dpst(j)=-.5/psi0
cjb      endif

cjb to get const dpdpsi/p  rigidrotor omstr   assuming rho proport to p
cjb  for alfrigid=1, this will be down only to e-2 at psi=2psi0, maybe ok
          if (doflat) then
              p2t(j)= exp(-alfrigid*psi(j)/psi0)  
              dp2dpst(j)=-exp(-alfrigid*psi(j)/psi0)*alfrigid/psi0
          endif

          p3(j)=(p3a+p3b*psir+p3c*psir**2+p3d*psir**3)*hp0(j)
          dp3dpsi(j)=(p3b+2.*p3c*psir+3.*p3d*psir**2)*hp0(j)/psim
          p2(j)=p2t(j)+p3(j)*dip
          dp2dpsi(j)=dp2dpst(j)+dp3dpsi(j)*dip
          hflr(j)=1.
 
c..... set p2 and p3 constant at large psi
          if (p2(j) <= p2flag) then
              p2(j)=p2floor
              p1(j)=p1floor
              dp2dpsi(j)=0.
              dp1dpsi(j)=0.
              hflr(j)=0.
          endif
          do 11 i=1,ix
c     if(abp(i).eq.0)then
c     b(i,j)=bvac(i)
c     go to 11
c     end if
              u1=p2(j)*abp(i)+p1(j)*abf(i)
              u2=p2(j)*bbp(i)+.5+p1(j)*bbf(i)
              u3=p2(j)*cbp(i)-bvac(i)**2*.5+p1(j)*cbf(i)


              if (i==1 .and. dowriteb) then
                  write(msg,'(a,2i5,3e14.6)') 'i,j,u1,u2,u3       ', \
                            i,j,u1,u2,u3
                  call baswline(STDOUT,msg)
              endif
    
c..... at,bt,ct coefficients for low density expansion
              at=u1
              bt=u2-.5
              ct=u3+bvac(i)**2*.5
              if (i == 1 .and. dowriteb) then
                  write(msg,'(a,2i5,3e14.6)') 'i,j,at,bt,ct       ', \
                            i,j,at,bt,ct
                  call baswline(STDOUT,msg)
              endif


              bv=bvac(i)
              ordera2=4.*ct*bt+4.*(bt*bv)**2+8.*at*ct*bv**2+12.*at*bt*bv**4 \
                     +8*at**2*bv**6
              bsq3=bv**2-2.*(ct+bt*bv**2+at*bv**4)+ordera2
              if (i >=44 .and. i<=65 .and. j==1 .and. dowriteb) then
                  write(msg,'(a,2i5,3e14.6)') 'i,j,bv,ordera2,bsq3', \
                            i,j,bv,ordera2,bsq3
                  call baswline(STDOUT,msg)
                  write(msg,'(a,2i5,2e14.6)') 'i,j,abs(ordera2/bsq3),epsp', \
                            i,j,abs(ordera2/bsq3),epsp
                  call baswline(STDOUT,msg)

              endif

              if (abs(ordera2/bsq3) <= epsp) then
                  if (bsq3 < 0) then
                      write(msg,'(a,4i10)')'neg bsq3: i,j,ordera2,bsq3', \
                                i,j,ordera2,bsq3
                      call baswline(STDOUT,msg)
                      call baswline(STDOUT,'at,bt,ct,bv**2,-2.*(ct+bt*bv**2+at*bv**4)')
                      write(msg,'(41x,2e14.6)') at,bt
                      call baswline(STDOUT,msg)
                      write(msg,'(41x,2e14.6)') ct,bv**2
                      call baswline(STDOUT,msg)
                      write(msg,'(41x,e14.6)') -2.*(ct+bt*bv**2+at*bv**4)
                      call baswline(STDOUT,msg)
                  endif
                  b(i,j)=sqrt(bsq3)
                  if (j == 1 .and. dowriteb) then
                      write(msg,'(a,i10,2e14.6)') 'bsq3: i,bvac(i)*ss,b(i,1)*ss', 
     *                           i,bvac(i)*ss,b(i,1)*ss
                      call baswline(STDOUT,msg)
                  endif
                  if (j == 2) bj2old(i)=sqrt(bsq3)
              else
                  
                  t1=-u2*.5/u1
                  radic2=(t1**2-u3/u1)
                  if (j==1 .and. i >= 50  .and. i<=65 .and. dowriteb) then
                      write(msg,'(a,2i10,2e14.6)') 'i,j,t1,radic2', \
                                i,j,t1,radic2
                      call baswline(STDOUT,msg)
                  endif


                  round=abs(radic2)/(t1**2+abs(u3/u1))
                  if (i == 1 .and. dowriteb) then
                      write(msg,'(a,2i10,e14.6)') 'i,j,round',i,j,round
                      call baswline(STDOUT,msg)
                  endif


                  if (round <= 1.e-6) radic2=0.
                  if (radic2 < 0) call errorend(i,j)
                  radic=sqrt(radic2)
                  bsq1=t1+radic
                  bsq2=t1-radic
                  if (j==1 .and. i >= 44 .and. i<=65 .and. dowriteb) then
                      write(msg,'(a,2i10,2e14.6)') 'i,j,t1,radic', 
     *                           i,j,t1,radic
                      call baswline(STDOUT,msg)
                      if(bsq1>0 .and. bsq2>0) then
                         write(msg,'(a,2i10,2e14.6)') 'i,j,bsq1**.5*ss,bsq2**.5*ss', 
     *                           i,j,bsq1**.5*ss,bsq2**.5*ss
                         call baswline(STDOUT,msg)
                      endif
                  endif
                  bsq=bsq1
cjb                  if (t1 > radic) bsq=bsq2
cjb replace latter with added decision re  bsq1 >< bvacsq
cjb this seems to prevent large oddball b,rzz,pperp   
cjb   for combo plug sloshing plus largebeta
cjb Still need to add same thing in the r,rzzloops later in this subr

                  if (t1 > radic .and. bsq1>bvac(i)**2) bsq=bsq2

                  if (j==1 .and. i >= 44 .and. i<=65  .and. dowriteb) then
                      write(msg,'(a,2i10,3e14.6)')'i,j,bsq1,bsq2,bsq', 
     *                           i,j,bsq1,bsq2,bsq
                      call baswline(STDOUT,msg)
                  endif

cjb save old b(i,2)
                  if (j == 2) bj2old(i)= sqrt(bsq)


                  if (bsq < 0) call errorend(i,j)
                  b(i,j)=sqrt(bsq)
                  if (j == 1 .and. dowriteb) then
                      write(msg,'(a,i10,2e14.6)') 'bsq12:  i,bvac(i)*ss,b(i,1)*ss', 
     *                           i,bvac(i)*ss,b(i,1)*ss
                      call baswline(STDOUT,msg)
                  endif
              endif

cjb need new ks p(b) model
              p2filt=1.
              isubpx=(ix+isubp)/2+1
              if (dosubp .and. i >= ibvmax(2)) then
                  p2filt=.5*(1.-tanh(alfsubp*(i-isubp)))      
                  if (i >= isubpx) p2filt=0.
              endif 
cjbrepl              if (i > ibvmax(2)) then
              if (i > ibvmax(1)) then

cjbrepl                  bsq= bvac(i)**2-2*poshinks*jperpz(i)*p2(j)
                  bsq= bvac(i)**2-2*posh(i)*jperpz(i)*p2(j)

                  b(i,j)=sqrt(bsq)
              endif


 11       continue
 10   continue
cjb endof 10,11 j,i loop   bsq->b(i,j)   both freis and jperpz coding 
cjb  jperpz code overwrites freis for i>ibvmax(1)   prev i>ibvmax(2)



 
c......calculate phi, the electric potential
      do i=1,ix
          arg1=((z(i)-z1)/(z1-z0))**2*(-xpot)*(1.-hzt0(i))
          arg2=wpot*((z(i)-z0)/(z0-z2))**2
          phi1(i)=phice*exp(arg1)+phipl/cosh(arg2)
      enddo
      do j=2,jx
          phi2(j)=(1.-(psi(j)*hp3(j)/psi3)**ypot)*hp3(j)
      enddo

c.....calculate pperp,  ppar and db/dpsi
c.....first calculate d coef. for ppar
      dp1=(ap1*bmax**2*4./3.+2.*bp1)*bmax
      dp2=(ap2-ap1)*bmn1**3/3.+(bp2-bp1)*bmn1+dp1+(cp1-cp2)/bmn1
      dp3=(ap3-ap2)*bvx2**3/3.+(bp3-bp2)*bvx2+dp2+(cp2-cp3)/bvx2

      write(msg,'(a,3e14.6)') 'dp1,dp2,dp3',dp1,dp2,dp3
      call baswline(STDOUT,msg)
      dt0=-8.*ct0/(bmax*3.)
      dt1=-8.*ct1/(bvx2*3.)
      dt2=-8.*ct2/(3.*bvx3)
      dt3=-8.*ct3/(3.*bm1)

cjb add dt2ks, dt3ks  mimic dt2,dt3
      dt2ks=-8.*ct2ks/(3.*bvoutks)
      dt3ks=-8.*ct3ks/(3.*bm1ks)

      write(msg,'(a,4e14.6)') 'dt0,dt1,dt2,dt3',dt0,dt1,dt2,dt3
      call baswline(STDOUT,msg)

      ppt=1./(1.-1./rp1**2)**2
      dter1=(-8./3-4.*(bvx3/bmn2)**3/3.+4.*(bvx3/bmn2))/bmn2
      dter2=-alsi*ppt*dter1*psls1
cjb dter1 not used;  dter2 not used except in nonsense zerosetting
cjb  ppt used only in dter2

cjb      write(msg,'(a,3e14.6)') 'ppt,dter1,dter2',ppt,dter1,dter2
cjb      call baswline(STDOUT,msg)




      do 20 i=1,ix
          do 21 j=1,jx
cjb i,j loop for pperp,ppar and related (i,j) grid arrays
cjb  both freis coding and jperpz coding, jperpz overwrites freis i>ibvmax(1)

              if (dter2 == 0) dter2=0.

              b2=b(i,j)**2
              b4=b2**2

              dter=dp1*hzp1(i)+dp2*hzp2(i)+dp3*hzp3(i)+dt1*hzt1(i) 
     *             +dt2*hzt2(i)-dt3*hzt3(i)+dt0*hzt0(i) 
     *             +dt2ks*hzt2ks(i)-dt3ks*hzt3ks(i)

cjb need to get rid of hardwired >= 50
              if (dosubp .and. i >= 50) then
                  p2filt=.5*(1.-tanh(alfsubp*(i-isubp)))
                  dter=dter*p2filt
                  if (i >= isubpx) then
                      dter=0.
                  endif
              endif

cjb add diagnostic array for dter
              dterjb(i)=dter

              ppart=(-abp(i)/3.*b4-bbp(i)*b2+cbp(i)+dter*b(i,j))
              ppert=(abp(i)*b4+bbp(i)*b2+cbp(i))
              pperte=(abq(i)*b4+bbq(i)*b2+cbq(i))
cjb retain these oldforms outside of newoldif   so can get pperpold,pparold
              if (j == 2) then
                  ppertold=ppert
                  ppartold=ppart
                  pperpold(i)=p2(2)*ppert
                  pparold(i)=p2(2)*ppart
              endif
              if (j == 25) then
                 pperpold25(i)=p2(25)*ppert
              endif

cjbrepl              if (i > ibvmax(2)) then
              if (i > ibvmax(1)) then

cjbrepl                  ppart=poshinks*jparz(i)
cjbrepl                  ppert=poshinks*jperpz(i)
                  ppart=posh(i)*jparz(i)
                  ppert=posh(i)*jperpz(i)

                  pperte=ppert
cjb      keeping same symbols should transmit better to code below
              endif

              pperp(i,j)=p2(j)*ppert
              ppar(i,j)=p2(j)*ppart



              if (pperp(i,j) < 0) then
                  write(msg,'(a,2i10,2e14.6)') 'i,j,pperp(i,j),p2(j),ppert', 
     *                       i,j,pperp(i,j),p2(j)
                  call baswline(STDOUT,msg)
                  write(msg,'(25x,e14.6)') ppert
                  write(msg,'(a,3e14.6)') 'abp(i),bbp(i),cbp(i),b4,b2', 
     *                       abp(i),bbp(i),cbp(i)
                  call baswline(STDOUT,msg)
                  write(msg,'(25x,2e14.6)') b4,b2
                  call baswline(STDOUT,msg)
                  call kaboom(0)
              endif
              qub(i,j)=  (b2-ppar(i,j)+pperp(i,j))/b(i,j)
              
              qubold(i)=(bj2old(i)**2-pparold(i)+pperpold(i))/bj2old(i)

cjb will have to read in separate rho model or for now just use rho ~ pperp

cjbrepl              if (i > ibvmax(2)) then
              if (i > ibvmax(1)) then

cjbrepl                  rho(i,j)=amass*(p2(j)**.5*noshinks*jperpz(i)+ 
                  rho(i,j)=amass*(p2(j)**.5*nosh(i)*jperpz(i)+ 
     *                      ncenter*cold*exp(-alfcold*(jx-j)))
              else
                  rho(i,j)=amass*(p2(j)**.5*(ebp(i)*b4+fbp(i)*b2+gbp(i))+ 
     *                      ncenter*cold*exp(-alfcold*(jx-j)))
              endif


cjb why using p2(j)**.5????  just making the tanh shape steeper yet?

cjb for doflat=true  make rho proportional to p2   ie take outthe sqrt
cjb for somereason that 'doflat' test applied only to freisplug
cjb   ie  not  jperpz ks
cjb    We havent used freis ks for a long time.

cjb  the following  elseif looking goofy:  if doflat false (default) then
cjb   dont set the rho at all!!?? 
cjb      well then the immediately preceding settng of rho takes over.
cjb need to clean this up an make the coding transparent
cjb  as it stands  following code always resets the rho for ksjperpz
cjb    but only if doflat=true does the following reset the rho for freisplug

cjb  might well want to just settle on one p2(j) model for rho
cjb   we could then take doflat out of this rho setting
cjb   need to check whereelse  doflat is used

cjbrepl              if (i > ibvmax(2)) then
              if (i > ibvmax(1)) then

cjbrepl                  rho(i,j)=amass*(p2(j)*noshinks*jperpz(i)+ 
                  rho(i,j)=amass*(p2(j)*nosh(i)*jperpz(i)+ 
     *                     ncenter*cold*exp(-alfcold*(jx-j)))
              elseif (doflat) then
                  rho(i,j)=amass*(p2(j)*(ebp(i)*b4+fbp(i)*b2+gbp(i))+ 
     *                     ncenter*cold*exp(-alfcold*(jx-j)))
              endif


cjb apr24 2002  put in new cold model variable in psi, maybe in z
cjb  trying to control the spurious mode seen localizing to high psi
cjb  E^(-0.3*(50 - j)) ~ .1 at j42
cjb  E^(-0.1*(50 - j)) ~ .1 at j25
cjb exp(-alfcold*(50-1)) =    9.52181D-01 at j=1   for alfcold=1e-3




              fac1=dp1dpsi(j)*(abf(i)*b4+bbf(i)*b2+cbf(i))
              fac2=dp2dpsi(j)*(abp(i)*b4+bbp(i)*b2+cbp(i))
              fac3=p2(j)*(2.*abp(i)*b2+bbp(i))+p1(j)*(2.*abf(i)*b2+bbf(i))
              dbdpsi(i,j)=-(fac1+fac2)/(b(i,j)*(1.+fac3*2.))
              if (j == 2) dbdpsiold(i)=dbdpsi(i,j)
              
cjbrepl              if (i > ibvmax(2)) then
              if (i > ibvmax(1)) then

                  fac1=0.
cjbrepl                  fac2=dp2dpsi(j)*poshinks*jperpz(i)
                  fac2=dp2dpsi(j)*posh(i)*jperpz(i)

cjb  prob want fac2 w p2filt also,  doit via p2filt on jperpz


cjbrepl                  fac3=p2(j)*poshinks*jdfz(i)/bvac(isubp)**2
cjbreplagain                  fac3=p2(j)*posh(i)*jdfz(i)/bvac(isubp)**2
                  fac3=p2(j)*posh(i)*jdfz(i)/bvacnorm(i)**2

                  dbdpsi(i,j)=-(fac1+fac2)/(b(i,j)*(1.+fac3*2.))
              endif
 
              dp2db=b(i,j)*(abp(i)*b2*2.+bbp(i))*2.
cjb              4*abp*b^3 +2*bbp*b
cjb             this is d/db of pperp, wo p2 factor
cjb             will have to transmit this from mma  ddbprpz(i)
              dp2dbe=b(i,j)*(abq(i)*b2*2.+bbq(i))*2.
              dp3db=dp2db-4.*abp(i)*b(i,j)**3/3.-2.*bbp(i)*b(i,j)+dter
cjb     looks like d/db  (pperp+ppar)    so need d/db of jparz also

              if (j == 2) then
                  dp3dbold=dp3db
              endif
              
cjbrepl              if (i > ibvmax(2)) then
              if (i > ibvmax(1)) then

cjbrepl                  dp2db=poshinks*ddbprpz(i)/bvac(isubp)
cjbreplagain             dp2db=posh(i)*ddbprpz(i)/bvac(isubp)
                  dp2db=posh(i)*ddbprpz(i)/bvacnorm(i)

                  dp2dbe=dp2db
cjbrepl                  dp3db=dp2db+poshinks*ddbparz(i)/bvac(isubp)
cjbreplagain             dp3db=dp2db+posh(i)*ddbparz(i)/bvac(isubp)
                  dp3db=dp2db+posh(i)*ddbparz(i)/bvacnorm(i)

              endif
cjb if i>ibvmax(1) this last overwrites the old version


              qubv(i,j)=(-dp2dpsi(j)*(ppert+ppart)-p2(j)*dp3db*dbdpsi(i,j))

              if (j == 2) then
                  qubvold(i)=(-dp2dpsi(2)*(ppertold+ppartold)- 
     *                        p2(2)*dp3dbold*dbdpsiold(i))
              endif

              pharg=psi(j)*hp3(j)
              dphidpsi=phi1(i)*ypot*(pharg/psi3)**(ypot-1.)/(-psi3)*hp3(j)
              epsi(i,j)=-dphidpsi

              dpdpsi=pperte*dp2dpsi(j)+p2(j)*dp2dbe*dbdpsi(i,j)

              omegci=ECHARG*b(i,j)/(amass*cee)
              unit=1./(sqrt(4.*Pi))
              omegstr=-unit*b(i,j)*dpdpsi/(omegci*rho(i,j))
              omeggb=unit*pperte*p2(j)*dbdpsi(i,j)/(omegci*rho(i,j))
              omegexb=+cee*dphidpsi

cjb doflat nearly const-in-z bvac   linear-in-psi p2 decr  tie exb to omstr
cjb aug5 2002  this should be default for omegexb
cjb  using only for doflat=true is really dangerous
cjb   we could easily try our tanh radial model (from doflat=false)
cjb   and end up with nonsense for omegexb

cjb      if(doflat) omegexb=exbratio*omegstr

              omegexb=exbratio*omegstr

cjb              if (i == 1) then
cjb                  write(msg,'(a,i10,3e14.6)') 'j,omegexb,omegstr,omeggb', 
cjb     *                       j,omegexb,omegstr,omeggb
cjb                  call baswline(STDOUT,msg)
cjb              endif
              if (i == 1) then
                  omstri1(j)=omegstr
                  omgbi1(j)=omeggb
              endif
              if (i == 34) omstri34(j)=omegstr
              if (i == 66) omstri66(j)=omegstr
              if (i == 148) omstr148(j)=omegstr
              
              if (j == 1)then
                  omebj1(i)=omegexb
                  omstrj1(i)=omegstr
              endif

              xxx(i,j)=hflr(j)*rho(i,j)*(omeggb-omegstr+2.*omegexb)*sf6
              yyy(i,j)=-hflr(j)*rho(i,j)*(omegexb+omeggb)*(omegexb-omegstr)*sf8
c....special yyy for cold halo case (flrm6)
c     yyy(i,j)=p2(j)*sf8

 
              ppertr=(abf(i)*b4+bbf(i)*b2+cbf(i))
              pperps(i,j)=pperp(i,j)+p1(j)*ppertr
              pperpe(i,j)=p1(j)*ppertr
cjb when p1=0  and abf,,=0 ppertr=pperpe=0   and pperps = pperp
cjb  pperps used only in errprp
cjb  i think can let these stand for new ks p(b)




 21       continue
 20   continue
cjb endof 20,21 ij loop   pperp,ppar,rho and related zgrid arrays



c......check for negative pressure, and terminate if necessary
      do j=1,jx
          do i=1,ix
              if (pperp(i,j) < 0 .or. pperpe(i,j) < 0) then
                  write(msg,101) i,j,pperp(i,j)
 101              format('terminated bc neg press,i,j,pperp',2i5,e12.4)
                  call fatal(msg)
              endif
          enddo
      enddo
 
c....calculate diagnostic on perpendicular pressure balance
cjb this errprp should be ok for either freis or ksjperpz
      do i=1,ix
          do j=1,jx
             errprp(i,j)=(b(i,j)**2-bvac(i)**2+2.*pperps(i,j))/bvac(i)**2
          enddo
      enddo


cjb following parallel pressure balance diagnostic seems only for freismodel
cjb   need equivalent for ks jperpz model also
cjb there seem to be other errors/inconsistencies wrt 'factor'

c....calculate diagnostic on parallel pressure balance
      do j=2,jx
          do i=2,ix-1
              delb=b(i+1,j)-b(i-1,j)
              if (abs(delb) < epsb) then

cjb epsb local, not in common or glr.v  unlike epsp
cjb epsb data set .01 in this subr:  appropriate value needs checking

                  abar=-abp(i)*5./3.
cjb abar is localscalar var;  danger to be misinterpreted as matrix abar(i,j)
                  bbar=-3.*bbp(i)
                  dter4=dp1*hzp1(i)+dp2*hzp2(i)+dp3*hzp3(i)+dt1*hzt1(i) 
     *                  +dt2*hzt2(i)-dt3*hzt3(i)+dt2ks*hzt2ks(i) 
     *                  -dt3ks*hzt3ks(i)
cjb  dter4 should be same as dter used prev for ppart
cjb  ppart=(-abp(i)/3.*b4-bbp(i)*b2+cbp(i)+dter*b(i,j)); ppar(i,j)=p2(j)*ppart
cjb 2nd line of dter4  seems to be missing +dt0*hzt0(i) term

                  dbar=dter4*2.
                  factor=-(abar*b(i,j)**4+bbar*b(i,j)**2+cbp(i)+dbar*b(i,j))

cjb  outside '-' sign  on factor looks wrong, cpt 'else' branch below
cjb  also this analytical form of factor= ddb(ppar*b)   missing p2(j)
cjb    --this says should see gross error anyplace this smalldelb branch used?
                  write(msg,'(a,2i10,e14.6)') 'small delb, i,j,delb',i,j,delb
                  call baswline(STDOUT,msg)

              else
                  factor=(ppar(i+1,j)*b(i+1,j)-ppar(i-1,j)*b(i-1,j))/delb
              endif

cjb  add firehose check d/db (Q/b) > 0    -> d/db qub >0    CFN(9)
cjb  dont use outside multiplier of b^2 as in CFN

              dqubdb(i,j)=(qub(i+1,j)-qub(i-1,j))/delb
cjb  might well need small delb 'analytical' branch  as for errprl

              errprl(i-1,j-1)=(pperp(i,j)-2.*ppar(i,j)+factor)/bvac(i)**2*2
          enddo
      enddo

cjb  this errprl  w the -2*ppar coupled w the  factor= d/db (ppar*b)  
cjb    looks different than but is equivalent to
cjb    the ppar pressure balance eqn in  CFN(7)  using d/db (ppar/b)



cjb why this i-1,j-1 ?    in glr.v errprl(IZX2,JRX-1)  -->  148,49 
cjb  this is consistent w doloop indices  i=2,ix-1    j=2,jx
cjb  the i-finite-difference rqrs this, but could still have used full dims;
cjb  means errprl(i,j) does not refer to same loc as b(i,j) or as errprp(i,j)
cjb  very awkward, very prone to misinterpretation

cjb it seems likely this implies same trtment for  xro(148*48)
cjb   ie  xro(1:48) for psiscan  implies need  rho,r,b( ,2:49)
cjb       xro(1:148 effectively) for  z scan  implies need  rho,r,b(2:149, )










c.....cal. r and rzz*r
      rzzrmax=0.
      do 40 i=1,ix
          sum1=0.
          sum2=0
          sum3=0.
          sum4=0.
          sum5=0.
          sum1new=0.
          sum2new=0.
          sum3new=0.
          sum5new=0.

          do 41 j=2,jx
              dps=dpsi(j)
              psib=psi(j-1)
              if (j == 2) then
                  dps=psi(2)
                  psib=0.
              endif
cjb dpsi(2) = 0, from way it is calcd
cjb any reason not to just correct dpsi(2) ??


              if (i <= 2) then
                  call myzero(ksimp,0.,hpk0,1)
                  call myzero(ksimp,0.,hpkm,1)
                  call myzero(ksimp,0.,hpkme,1)
                  do k=1,ksimp
                      psik=psib+dps*(k-1)/(ksimp-1)
                      if (psik < psi0) hpk0(k)=1.
                      if (psik < psim) hpkm(k)=1.
                      if (psik < psime) hpkme(k)=1.
                  enddo
c..... cal. p1 (p1k) at intermediate points in dpsi interval
                  do  k=1,ksimp
                      psik=psib+dps*(k-1)/(ksimp-1)
                      psikr=psik/psim
                      psikre=psik/psime
                      p1ek=(ae1+be1*psikre+ce1*psikre**2+de1*psikre**3)*hpkme(k)
                      p2ek=.5*(1.-tanh((psik-psi0e)/wp2e))*(1.-hpkme(k))
                      p1k(k,j)=(p1ek+p2ek)/p1max
                      p3k=(p3a+p3b*psikr+p3c*psikr**2+p3d*psikr**3)*hpk0(k)
                      p2k(k,j)=(1.-tanh((psik-psi0)/wp2))*.5+p3k*dip

cjb apr602 correct p2k to use the exp form if doflat true
                      if (doflat) p2k(k,j)= exp(-alfrigid*psik/psi0)
                  enddo
              endif

              do 32 k=1,ksimp
                  if (p2k(k,j) <= p2flag) then

                      p2k(k,j)=p2floor  
                      p1k(k,j)=p1floor
                  endif
                  u1=p2k(k,j)*abp(i)+p1k(k,j)*abf(i)
                  u2=p2k(k,j)*bbp(i)+p1k(k,j)*bbf(i)+.5
                  u3=p2k(k,j)*cbp(i)+p1k(k,j)*cbf(i)-bvac(i)**2*.5
c..... at,bt,ct coefficients for low density expansion
                  at=u1
                  bt=u2-.5
                  ct=u3+bvac(i)**2*.5
                  bv=bvac(i)
                  ordera2=4.*ct*bt+4.*(bt*bv)**2+8.*at*ct*bv**2 
     *                   +12.*at*bt*bv**4+8*at**2*bv**6
                  bsq3=bv**2-2.*(ct+bt*bv**2+at*bv**4)+ordera2
cjb       if((ordera2/bsq3).le.epsp) then
                  if (abs(ordera2/bsq3) <= epsp) then
cjb      this branch if  either bsq3<0 or ordera2<0   one only
                  bk=sqrt(bsq3)
              else
                  t1=-u2*.5/u1
                  radic2=(t1**2-u3/u1)
                  round=abs(radic2)/(t1**2+abs(u3/u1))
                  if (round <= 1.e-6) radic2=0.
                  radic=sqrt(radic2)
                  bsq1=t1+radic
                  bsq2=t1-radic
                  bsq=bsq1
                  if (t1 > radic .and. bsq1>bvac(i)**2) bsq=bsq2
                  bk=sqrt(bsq)
              endif


              p2filt=1.
              if (dosubp .and. i >= 50) then
                  p2filt=.5*(1.-tanh(alfsubp*(i-isubp)))      
                  if (i >= isubpx) p2filt=0.
              endif 


cjbrepl              bsqjb=bvac(i)**2-2*poshinks*jperpz(i)*p2k(k,j)
              bsqjb=bvac(i)**2-2*posh(i)*jperpz(i)*p2k(k,j)


cjb p2filt in jperpz 

              bkjb=sqrt(bsqjb)

              capf=1.+2.*p1k(k,j)*(2.*abf(i)*bk**2+bbf(i))+2.*p2k(k,j)* 
     *              (2.*abp(i)*bk**2+bbp(i))

cjb  capf  = 1. +2.*p2k*(2.*abp*bk**2+bbp)    ignoring p1k eringterm
cjb  capf  = 1. +2*dpperp/dbsq

cjbrepl              capfnew=1.+2.*p2k(k,j)*poshinks*jdfz(i)/bvac(isubp)**2
cjbreplagain              capfnew=1.+2.*p2k(k,j)*posh(i)*jdfz(i)/bvac(isubp)**2
cjb must need different bvac(isubp) for plug vs ks  as well as posh
              capfnew=1.+2.*p2k(k,j)*posh(i)*jdfz(i)/bvacnorm(i)**2


cjb do i need posh(i) set in ccell  for capfnew?  prev capfnew=0 if jdfz=0

              deli1(k)=1./bk
              deli2(k)=1./(bk**3*capf)
              deli3(k)=deli2(k)/(capf*bk**2)
              deli4(k)=p1k(k,j)/(capf*bk)**3
              deli5(k)=p2k(k,j)/(capf*bk)**3


              del1new(k)=1./bkjb
              del2new(k)=1./(capfnew*bkjb**3)
              del3new(k)=del2new(k)/(capfnew*bkjb**2)
              del5new(k)=p2k(k,j)/(capfnew*bkjb)**3  

 32           continue         


cjb    deli2 =     1/[bk^3*(1.+2.*p2k*(2.*abp*bk**2+bbp))] 
cjb                            dpperp/dbsq = p2k*(2*abp*bksq +bbp)
cjb                        (1.+2*...) = (1+2*dpperp/dbsq)

cjb    deli3 = deli2/[bk^2*(1.+2.*p2k*(2.*abp*bk**2+bbp))]
cjb          =     1/[bk^5*(1.+2.*p2k*(2.*abp*bk**2+bbp))^2]

cjb    deli5 = p2k  /[bk^3*(1.+2.*p2k*(2.*abp*bk**2+bbp))^3]
cjb    deli5  leading p2k is just the dimless p2(j)  no p(b) in there
cjb   later on  capi5 is multiplied by  abp(i)
cjb   capi1,capi2,capi3 _not_ multiplied later by any pressure coeff(i)

cjb                (1/2bk)*dpperp/dbk =p2k*(1/2bk)*[4*abp*bk^3 +2*bbp*bk]
cjb                                             = p2k*(2*abp*bk^2 +bbp)
cjb                                             = dpperp/dbsq

cjb                              d2pperp/dbsq^2 = p2k*2*abp

cjb     dp/db    = p2k[4abp*bk^3 +2bbp*bk]
cjb     d2p/db^2 = p2k(12abp*bk^2 +2*bbp)
cjb    {d2p/db^2 -(1/b)dp/dbk} = p2k*(12abp*bk^2 +2*bbp -[4*abp*bk^2 +2*bbp])
cjb                            = p2k*  8abp*bk^2

cjb     so d2p/dbsq^2 = p2k*2*abp = p2k*8abp*bk^2/(4bk^2)
cjb                               = { }/(4bk^2)

cjb     also, rzz3 factor  p2k*8*abp = 4*d2p/dbsq^2  <<< checks w TeX notes???


cjb  dp/dbsq = dp/(2b db) = (1/2b) dp/db
          call simps(deli1,dcapi1,ksimp,dps)
          call simps(deli2,dcapi2,ksimp,dps)
          call simps(deli3,dcapi3,ksimp,dps)
          call simps(deli4,dcapi4,ksimp,dps)
          call simps(deli5,dcapi5,ksimp,dps)

          call simps(del1new,dcap1new,ksimp,dps)
          call simps(del2new,dcap2new,ksimp,dps)
          call simps(del3new,dcap3new,ksimp,dps)
          call simps(del5new,dcap5new,ksimp,dps)

       
          capi1=dcapi1+sum1
          capi2=dcapi2+sum2
          capi3=dcapi3+sum3
          capi4=dcapi4+sum4
          capi5=dcapi5+sum5

          cap1new=dcap1new+sum1new
          cap2new=dcap2new+sum2new
          cap3new=dcap3new+sum3new
          cap5new=dcap5new+sum5new

cjb  to manage needing old form for ccplg and newform for ks
cjb  retain old coding here also;  assign oldvalue and newvalue to tempscalars
cjbold
          roldpg=sqrt(2.*capi1)
cjbksnew
          rnewks=sqrt(2.*cap1new)

cjbrepl          if (i > ibvmax(2)) then
          if (i > ibvmax(1)) then
              r(i,j)= rnewks
          else      
              r(i,j)=roldpg
          endif
cjb this might be faulty?  does rzz have same problem?
cjb  i think problem is the rnewks picks up only vac terms in ccpg
cjb    this means of course should not use rnewks when in ccpg
cjb     except if understand that rnewks in ccpg will use bvac
cjb  Here the above iftest uses  only rnewks in ks and only roldpg in pgcc

cjb  it seems cap1new must be same as capi1 in ccpg  NO
cjb   NO, cap1new uses bvacinpg ie noselfb in pg
cjb  this procedure must extend to all the capi vs capnew

cjb   but if they are really separate in i  then why isnt this ok as is?
cjb   ie if these r calcs are totally indep in each i no linkage to lower i
cjb  then we should be ok.


cjb  r= sqrt[2*simps(1/bk)]    why cant we calc 2nd deriv of this?  toonoisy?



          vv1=(bvac(i)*dbvdz(i))**2
          vv2=bvac(i)*d2bvdz2(i)+dbvdz(i)**2

cjbold
          rzzoldpg=-vv1*capi2**2/r(i,j)**3+3.*vv1*capi3/r(i,j) 
     *             -vv2*capi2/r(i,j)+8.*bvac(i)**2*dbvdz(i)**2*(abf(i)*capi4 
     *             +abp(i)*capi5)/r(i,j)

cjbnewks
          rzznewks=-vv1*cap2new**2/r(i,j)**3+3.*vv1*cap3new/r(i,j) 
     *             -vv2*cap2new/r(i,j)+ 
     *              4.*bvac(i)**2*dbvdz(i)**2* 
     *              (posh(i)*jd2fz(i)/bvacnorm(i)**4*cap5new)/r(i,j)

cjbreplagain                   (posh(i)*jd2fz(i)/bvac(isubp)**4*cap5new)/r(i,j)
cjbrepl                   (poshinks*jd2fz(i)/bvac(isubp)**4*cap5new)/r(i,j)

cjb  to manage needing old form for ccplg and newform for ks
cjb  retain old coding here also;  assign oldvalue and newvalue to tempscalars
cjbrepl          if (i > ibvmax(2)) then
          if (i > ibvmax(1)) then
              rzz(i,j)=rzznewks
          else      
              rzz(i,j)=rzzoldpg
          endif


cjb for comparison w old model in ks, at j=2 calc old rzz pieces
          if (j == 2) then
              rold(i)=sqrt(2*capi1)
              rzz1a(i)=-vv1*capi2**2/rold(i)**3
              rzz1b(i)=3*vv1*capi3/rold(i)
              rzz2(i)=-vv2*capi2/rold(i)
              rzz3(i)=8*bvac(i)**2*dbvdz(i)**2*(abp(i)*capi5)/rold(i)
              rzzold(i)=rzz1a(i)+rzz1b(i)+rzz2(i)+rzz3(i)

cjb     p2k in capi5            from deli5
cjb     1/[bk^3 (1+2dp/dbsq)^3] from deli5
cjb     8*p2k*abp = 4* d2p/dbsq^2

cjb   only rzz3 has abp(i) factor,
cjb   ie rzz1a,rzz1b,rzz2 have no xtra  pperp(i) coeff


cjb  nolonger use rnew except for rd2new
              rnew(i)=sqrt(2*cap1new)


              rzz1anew(i)=-vv1*cap2new**2/r(i,j)**3
              rzz1bnew(i)=3*vv1*cap3new/r(i,j)
              rzz2new(i)=-vv2*cap2new/r(i,j)
              rzz3new(i)=4*bvac(i)**2*dbvdz(i)**2* 
     *                     (posh(i)*jd2fz(i)/bvacnorm(i)**4*cap5new)/r(i,j)

cjbreplagain              (posh(i)*jd2fz(i)/bvac(isubp)**4*cap5new)/r(i,j)
cjbrepl                   (poshinks*jd2fz(i)/bvac(isubp)**4*cap5new)/r(i,j)


cjb     8*p2k*abp = 4* d2p/dbsq^2   p2k factor inside capi5    
cjb  now 4*poshinks*jd2fz(i) takes place of  8*abp(i)   p2k factor in cap5new

              rzznew(i)= rzz1anew(i)+rzz1bnew(i)+rzz2new(i)+rzz3new(i)
cjb now rzznew(i) should be identical to rzz(i,2)


cjb  sumof 1anew+1bnew+2new has vac rzz  in ccplug
cjb   this means for combined model  cant simply add oldcode +new code
cjb   we could add oldcode plus new code if use heaviside fns 
cjb   here instead trying if tests on r(i,j) and on rzz(i,j)
cjb        perhaps clumsy but should work.

              icapi2(i)=capi2
              icapi3(i)=capi3
              icapi5(i)=capi5
              icap2new(i)=cap2new
              icap3new(i)=cap3new
              icap5new(i)=cap5new
 

          endif


cjb  vv1 = bv^2*dbvdz^2           dims: b^4/L^2
cjb  vv2 = bv*d2bvdz2 + dbvdz^2   dims: b^2/L^2 + b^2/L^2

cjb  rzz = -vv1*capi2^2/r^3  +vv1*3.*capi3/r
cjb        -vv2*capi2  /r
cjb        +8bvac^2*abp*capi5                      
cjb                 ^^^ ~dpperp/db^4  p2(j)or p2k factor is in capi5
cjb                      implies restricted to b^4 pperp model?

cjb  rzz =       -vv1  /{r^3*b^6*(1+2*dp/dbsq)^2}   dims:b^4/L^2*1/(b^6*L^3)
cjb                                                      1/(b^2 *L^5) ??????
cjb              +vv1*3/{r  *b^5*(1+2*dp/dbsq)^2}   dims:b^4/L^2*1/(b^5*L)
cjb                                                      1/(b   *L^3) ????
cjb              -vv2  /{r  *b^3*(1+2*dp/dbsq)}     dims: b^2/L^2 *1/(b^3*L)
cjb                                                      1/(b   *L^3) 
cjb  +8bvac^2*dp/db^4 /{    b^3*(1+2*dp/dbsq)^3}/r   dims b^2 *p/b^4 /b^3/r
cjb                                                      1/(b^3   L)

cjb if the simps terms all need a  L^2*b  (from dps)
cjb   1st vv1 term ~ capi2^2   needs dps^2         1/L
cjb                                                1/L
cjb                                                1/L
cjb                                                L/b^2
cjb  only 1st,2nd,3rd terms consistent w dims of rzz ~ 1/L
cjb   adding dbvdz(i)**2 to 4th term  supplys b^2/L^2  --> net of  1/L  ok


          if (dorzzlb) rzz(i,j)=-r(i,j)/(lb*lb)

          sum1=+capi1
          sum2=capi2
          sum3=capi3
          sum4=capi4
          sum5=capi5
          sum1new=cap1new
          sum2new=cap2new
          sum3new=cap3new
          sum5new=cap5new

          rzzr=r(i,j)*rzz(i,j)

          rzzrmax=max(rzzrmax,rzzr)
 41       continue
 40   continue
cjb endof 40,41 i,j loop for  r,rzz





      do i=1,ix
          r(i,1)=r(i,2)
      enddo

cjb 2nd finite diff 
      do i=3,isubp
cjb rd2 should now have correct total rold in ccpg and rnew in ks
          rd2diff(i)=(r(i+1,2)-2*r(i,2)+r(i-1,2))/dz(i)**2
cjb rnewd2 now should be same in ks as rd2  but just vac in ccpg  OK
          rnewd2(i)=(rnew(i+1)-2*rnew(i)+rnew(i-1))/dz(i)**2
cjb  roldd2 should have old r thruout
          roldd2(i)=(rold(i+1)-2*rold(i)+rold(i-1))/dz(i)**2
      enddo








c..... calculate diagnositc quantities flute1, flute2 and flute3
c..... and xxxfreq
      do 60 j=2,jx-1
          is=2
          is1=0
          if (mod(ix-2,2) == 0) is1=1

cjb is1=1 for evenvals of ix-2
cjb iloop is  i=2,ix-2
cjb  ia = i -2+1 = i-1    or ia=1,ix-3

          do i=is,ix-1-is1
              ia=i-is+1
cjb eterm only used further on  *ringj   so can ignore for  new p(b) ?
              eterm=dp2dpsi(j)*(abp(i)*b(i,j)**4+bbp(i)*b(i,j)**2+cbp(i)) 
     *              +p2(j)*(4*abp(i)*b(i,j)**3+2*bbp(i)*b(i,j))*dbdpsi(i,j)
              if (doflat) then
cjb this branch for p2(j) = exp(-alfrigid*psi)  _and_  rho ~ p2
                  rterm=0.
                  if (dp2dpsi(j) .ne. 0) then
                      rterm=dp2dpsi(j)*(ebp(i)*b(i,j)**4+fbp(i)*b(i,j)**2 
     *                      +gbp(i))+p2(j)*(4*ebp(i)*b(i,j)**3 
     *                      +2*fbp(i)*b(i,j))*dbdpsi(i,j)
                  endif
              else
cjb this branch original; uses tanh model for p2(j)  _and_ rho ~ p2**.5
cjb   effect 1st 2 stmnts:    rterm=0 if dp2dpsi=0
cjb   not sure why nec, prevents using 2nd term(~dbdpsi) if dp2dpsi=0
cjb    dp2dpsi can =0 via p2floor   
                  rterm=0.
                  if (dp2dpsi(j) .ne. 0) then
                      rterm=dp2dpsi(j)*(ebp(i)*b(i,j)**4+fbp(i)*b(i,j)**2 
     *                      +gbp(i))/(2.*p2(j)**.5) 
     *                      +p2(j)**.5*(4*ebp(i)*b(i,j)**3 
     *                      +2*fbp(i)*b(i,j))*dbdpsi(i,j)
                  endif
              endif
cjb rterm used in droterm   which is used in simps integral for  ans5 -> droave
cjb   droave used in trad for gamwkbsq
cjb seems like we definitely need new p(b) code for rterm
              droterm(ia)=rterm*amass/b(i,j)**2/uuz(i)
              ringj=dp1dpsi(j)*(abf(i)*b(i,j)**4+bbf(i)*b(i,j)**2+cbf(i)) 
     *              +p1(j)*(4*abf(i)*b(i,j)**3+2*bbf(i)*b(i,j))*dbdpsi(i,j)
              dflute1(ia)=+rzz(i,j)*qubv(i,j)/(r(i,j)*b(i,j)**2)/uuz(i) 
     *                    +eterm*ringj/b(i,j)**3/uuz(i)
              dflute2(ia)=yyy(i,j)/((r(i,j)*b(i,j))**2*b(i,j))/uuz(i)
              dflute3(ia)=rho(i,j)/((r(i,j)*b(i,j))**2*b(i,j))/uuz(i)
              dflute4(ia)=xxx(i,j)/((r(i,j)*b(i,j))**2*b(i,j))/uuz(i)
cjb for equivalent of Post stability integral, want ~ dflute1(ia)
cjb  but replace qubv ~ dpdpsi  by  
cjb           dimless ptot= (pperp(i,2)+ppar(i,2))/(pperp(137,2)+ppar(137,2))
cjb is it clear we want ccell and plug also normalized to ptot(inj=137)?
cjb  if we really want equiv of Post 'unit pressure' instab integral
cjb  then we need dflute1a w qubv->1
          enddo

cjb at end of iloop 61  ia = ix-3
cjb the following simps calls are being done successively for each j
cjb  which means eg that ans1 from j=2 is overwritten by ans1 from j=3,...

cjb there does not appear to be any sum or integral over j  ie over psi

          call simps(dflute1,ans1,ia,du)
          call simps(dflute2,ans2,ia,du)
          call simps(dflute3,ans3,ia,du)
          call simps(dflute4,ans4,ia,du)
          call simps(droterm,ans5,ia,du)
          flute1(j)=ans1*(ia-1)+quad1(dflute1(1),dflute1(2),z(2),z(3))
cjb wo the boundary quad1 term    flute1(j) = ans1*(ia-1)  = ans1*(ix-4)

          flute2(j)=(ans2*(mm**2-1)+ans1)*(ia-1)
          flute3(j)=-ans1/ans3
cjb wo bdy quad1 term flute3 = -simps(dflute1)/simps(dflute3)
cjb  dflute1 ~ rzz dp/dpsi/(r*b^2)  ~ rzz  vtsq drho/dpsi/(r*b^2)
cjb  dflute3 ~ rho        /(r^2*b^3)
cjb  so ~ dflute1/dflute3 ~vtsq*rzz * r*b (1/rho)drho/dpsi
cjb  drho/dpsi = drho/dr *dr/dpsi     dpsi/dr = r*b   so drho/dr =r*b drho/dpsi
cjb     ~ dflute1/dflute3      ~vtsq*rzz/L_rho    1/L_rho = (1/rho)drho/dr
cjb  ok consistent w  gamsq;   the division by rho is necessary

cjb  what looks like if dont divide the flute1 by flute3 until after integrals
cjb   simps(dflute1) ~ rzz vtsq  rho/psip/(rb^2)
cjb   simps(dflute3) ~ rho/(r^2*b^3)

cjb using bcohen r^2b constancy  so replace 1/b^2 ->r^4  in simps(dflute1)
cjb  then net facor is  ~ r^3 rzz  like Post
cjb  so it looks like post scaling is something like  simps(dflute1)
cjb  ie before dividing by simps(dflute3)    using constancy of psip


          flute3(j)=-flute1(j)/(ans3*(ia-1)+quad1(dflute3(1),dflute3(2), 
     *                z(2),z(3)))
          rhoave(j)=ans3*(ia-1)+quad1(dflute3(1),dflute3(2),z(2),z(3))
          xxxave(j)=ans4*(ia-1)+quad1(dflute4(1),dflute4(2),z(2),z(3))
          yyyave(j)=ans2*(ia-1)+quad1(dflute2(1),dflute2(2),z(2),z(3))
          droave(j)=ans5*(ia-1)+quad1(droterm(1),droterm(2),z(2),z(3))
          xxxfreq(j)=dt*xxxave(j)/rhoave(j)
 60   continue
cjb endof 60jloop  and internal i->ia loop
cjb     doing more simpson stuff for diagnostics



cjb  looks like the instabintegral wo the division by the rhointegral is flute1
cjb  if I want the cumulative integral in z  need to work with dflute1
cjb  maybe can get the cumintegral from simps(dflute1,ans1,iaincr,du)
cjb   it looks like that simps call would have to be inside the iloop




cjb for equivalent of Post stability integral, want dpost2~ dflute1(ia)
cjb  but replace qubv ~ dpdpsi  by  
cjb           dimless ptot= (pperp(i,2)+ppar(i,2))/(pperp(137,2)+ppar(137,2))
cjb is it clear we want ccell and plug also normalized to ptot(inj=137)?

cjb  if we really want equiv of Post 'unit pressure' instab integral
cjb  then we need dpost1 ~ dflute1a w qubv->1

cjb  rather than imbedding more j arrays to get Post equivalents
cjb  instead put in  a separate ix loop   like loop 61 above

      is=2
      is1=0
      if (mod(ix-2,2) == 0) is1=1

cjb is1=1 for evenvals of ix-2
cjb iloop is  i=2,ix-2   why ix-2?
cjb  ia = i -2+1 = i-1    or ia=1,ix-3   odd num, rqrd for simps?

      do 961 i=is,ix-1-is1
          ia=i-is+1
          pnorm=1.

          if (pperp(isubp,2)+ppar(isubp,2) > 0) then
              pnorm=2*(pperp(isubp,2)+ppar(isubp,2))
          endif

          ptot(i)=(pperp(i,2)+ppar(i,2))/pnorm

cjb prefer pnorm = 2*ptot(isubp)   the tanh pfilt cuts to .5 at isubp
cjb  using isubp takes guesswork out of establishing equiv to 137


cjb  using  rzz*r^3 integrand     dpost1r  is unit press  dpost2r is ptot wtd
          dpost1r(ia)=        rzz(i,5)*r(i,5)**3/uuz(i)*1.e-4
          dpost2r(ia)=ptot(i)*rzz(i,5)*r(i,5)**3/uuz(i)*1.e-4
          dpostjb(ia)=(pperp(i,2)+ppar(i,2))* 
     *                 rzz(i,5)*r(i,5)**3/uuz(i)*1.e-4
          dpostjbo(ia)=(pperpold(i)+pparold(i))* 
     *                  2.64575*rzzold(i)*(2.64575*rold(i))**3/uuz(i)*1.e-4

cjb      dpostjbo(ia)=(pperpold(i)+pparold(i))*
cjb     .                    rzzold(i,5)*rold(i,5)**3/uuz(i)*1.e-4
cjb                    uhhcant dothis j=5 on old
cjb   maybe could just use 4.69*r(i,2)?  but what about rzz(i,5}?
cjb with betslsh=1.e-6 
cjb r(34,2:6)          shape: (1,5)
cjb  (2)           1.77553D+00   3.07530D+00   3.97020D+00   4.69760D+00
cjb  (6)           5.32658D+00    4.7+5.3 /2 = 5.0
cjb rold(34) =    1.77553D+00
cjb   so if  try this  kind ofthing factor needed on r(j=2) =4.6976/1.7755

cjb rzzold(34) =   -3.14011D-03
cjb rzz(34,2:6)        shape: (1,5)
cjb (2)          -3.14011D-03  -5.43884D-03  -7.02151D-03  -8.30796D-03
cjb (6)          -9.42034D-03
cjb 8.30796D-03/3.14011D-03
cjb 8.30796D-03/3.14011D-03 =    2.64575D+00
cjb 4.6976/1.7755
cjb 4.6976/1.7755 =    2.64579D+00    ok same factor on both

cjb used j=5 bc r(35,5) =    4.70277D+00 
cjb             r(35,6) =    5.33244D+00
cjb could do better i guess if use avg of j=5,j=6 

cjb Post uses r0=5cm = .05m  his entire calc is in mks
cjb   our  rzz*r^3 is a net  L^2   so should need a factor 10^-4
cjb  seems to approx work:

cjb   post1r = 7.52374D+03    7.52e3/10^6 = .0075   not far from Post .005
cjb   rationale is  rzz*r^3 *dz = r^4/z^2 *dz = r^4/z = L^3  or 10^6

cjb dpost1r(147) = 1.73051D+02   173/10^4 = .0173   not far from Post .013
cjb  rationale  is  rzz*r^3 = r^4/z^2 = L^2  or 10^4  

 961  continue
cjb end of 961 iloop using j=2 for poststab integrals




      call simps(dpost1r,pans1r,ia,du)
      call simps(dpost2r,pans2r,ia,du)
      call simps(dpostjb,pansjb,ia,du)

      post1r=pans1r*(ia-1)*1.e-2
      post2r=pans2r*(ia-1)*1.e-2
      postjb=pansjb*(ia-1)*1.e-2

cjb check that simpson integral w sum over z
      post1sum=0.
      post2sum=0.
      postjbsum=0.
      postjbsumo=0.

cjb      do iia=1,147
      do iia=1,ix-3
          post1sum=post1sum+dpost1r(iia)*du*1.e-2
          post2sum=post2sum+dpost2r(iia)*du*1.e-2
          postjbsum=postjbsum+dpostjb(iia)*du*1.e-2
          postjbsumo=postjbsumo+dpostjbo(iia)*du*1.e-2
      enddo

cjb this ignores flora quad1 treatment using z(2),z(3)  implement later
      write(msg,'(a,e14.6)') 'post1r unitpressure instab integral', post1r
      call baswline(STDOUT,msg)
      write(msg,'(a,2e14.6)') 'post2r ptot/ptotinj instab integral, pnorm ', 
     *           post2r,pnorm
      call baswline(STDOUT,msg)

cjb   the instab integrand?  is it simply dpost1,2(ia)?
cjb    ia=1,ix-3 corresponds to i,2,148   altho the final 2  not important?

cjb   if want an integrand over just the ks region:
cjb      either mess w the i -> ia indexing  
cjb             ia needed to start from 1 simps?  not clear,  iamax=ix-3 why?
cjb      or     modify the integrands  with a  stepfn filter 
cjb   it might be better to get rid of the simps thing for all z integrals
cjb   i dont think there is subdivsion of deltaz going on here like there
cjb   is on the simps treatment of psi integrals for r,rzz, etc

cjb  to do my own integrals over ia    need  du equiv for dz

      post1cc=0.
      post2cc=0.
      postjbcc=0.
      do iia=1,ibvmax(1)-1
          post1cc=post1cc+dpost1r(iia)*du*1.e-2
          post2cc=post2cc+dpost2r(iia)*du*1.e-2
          postjbcc=postjbcc+dpostjb(iia)*du*1.e-2
      enddo

      post1pg=0.
      post2pg=0.
      postjbpg=0.

cjb                   bvac max i51 nearly same i50
      do iia=ibvmax(1),ibvmax(2)-1    
          post1pg=post1pg+dpost1r(iia)*du*1.e-2
          post2pg=post2pg+dpost2r(iia)*du*1.e-2
          postjbpg=postjbpg+dpostjb(iia)*du*1.e-2
      enddo



cjb  for the ptot ints probably better to rely on full ia,1,147
cjb  for unitp ints probably best to use isubp  or ia=isubp-1
cjb     this is particularly so if we identify isubp with Postzmax

      p1kssubp=0.
      do iia=ibvmax(2),isubp-1
          p1kssubp=p1kssubp+dpost1r(iia)*du*1.e-2
      enddo


cjb now repeat but extend integral to maxz; use maxz used by simps
cjb  note p2ksmax still uses     2ptotisubp for renormalization
      p1ksmax=0.
      p2ksmax=0.
      pjbksmax=0.
      pjbksmaxo=0.

      do iia=ibvmax(2),ix-3
          p1ksmax=p1ksmax+dpost1r(iia)*du*1.e-2
          p2ksmax=p2ksmax+dpost2r(iia)*du*1.e-2
          pjbksmax=pjbksmax+dpostjb(iia)*du*1.e-2
          pjbksmaxo=pjbksmaxo+dpostjbo(iia)*du*1.e-2
      enddo


cjb using ptot normalized w ptot(137) means ks integral wont change w betks
cjb    instead the plug and cc pieces will reduce with increased betks.
cjb    That is nonintuitive and prone to error for interpretation.

cjb there is no absolute need to normalize ptot at all,
cjb    except for comparison to Post mma where he normalizes to 2 or 1 at inj.
cjb we could also calc postintegrals  without the normalization,
cjb    perhaps more useful in some sense, certainly more like flute3
cjb    To do this need an unnormalized ptot  for a new dpostjb(ia)


cjb  set up beta array both normalized, unnormalized
cjb   for normalized use pperpnorm=2x pperpisubp   2x bc pfilt  just like pnorm
cjb  be aware of pperpks model sensitivity to pfudgeks  
cjb     betaks plot much different eg pfudgeks = 1.9995  vs 1.99
cjb  even for Post p(b), pperp approaches 0+  at injection
      pperpinj=2*pperp(isubp,2)
      btksinj=2*pperpinj/bvac(isubp)**2
      pparinj=2*ppar(isubp,2)

      do i=2,ix
          beta(i)=2*pperp(i,2)/bvac(i)**2
          betanorm(i)=beta(i)/btksinj 
          betaold(i)=2*pperpold(i)/bvac(i)**2
      enddo
      btksinjo=4*pperpold(isubp)/bvac(isubp)**2


      write(msg,'(a,2e14.6)') 'btksinj=2*pperpinj/bvac(isubp)**2,btksinjo ',
     *                          btksinj,btksinjo
      call baswline(STDOUT,msg)

      write(msg,'(a,3e14.6)') 'postjbpg,pjbksmax,pjbksmaxo ',
     *                          postjbpg,pjbksmax,pjbksmaxo
      call baswline(STDOUT,msg)

      sumpostint=postjbpg+pjbksmax
      write(msg,'(a,3e14.6)') 'sumpostint,postjbsum,postjb ',
     *                          sumpostint,postjbsum,postjb 
      call baswline(STDOUT,msg)

       

      grow=0.
      growmax=0.
      xfreqmax=0.
      do j=2,jx-1
          if (flute3(j) > 0) grow=grow+flute3(j)
          growmax=max(growmax,flute3(j))

          xxxabs=abs(xxxfreq(j))
          xfreqmax=max(xfreqmax,xxxabs)
      enddo
      write(msg4,102) grow,growmax,xfreqmax,rzzrmax,p2(jx)
 102  format('grow=',e14.6,3x,'growmax=',e14.6,/'maximum xxxfreq*dt= ', 
     *        e14.6,4x,'(flr real frequency times dt)',/, 
     *        'maximum rzz*r= ',e14.6 ,/'relative wall pressure p2(jx)=',e14.6)
      do l=1,4
          call baswline(STDOUT,msg4(l))
      enddo
                 
c..... local growth rate for high mm, (wkb approx. )
      do j=2,jx-1
          omegwkb=mm*.5*xxxave(j)/rhoave(j)
          trad=.25*mm**2*(xxxave(j)**2+4.*rhoave(j)*yyyave(j)) 
     *         -rhoave(j)**2*flute3(j)-droave(j)*yyyave(j)
          if (trad <= 0) then
              gamwkb(j)=sqrt(-trad)/rhoave(j)
              omeg1wkb(j)=omegwkb
              omeg2wkb(j)=0.
          else
              omeg1wkb(j)=omegwkb+sqrt(trad)/rhoave(j)
              omeg2wkb(j)=omegwkb-sqrt(trad)/rhoave(j)
          endif
      enddo

      return
      end # equiltm

      subroutine errorend(i,j)
      implicit none

c**** This sub. is called if a negetive square root is generated in
c**** equiltm during the calculation of b.
c**** A diagnostic message is sent out and the problem is terminmated.

      Use(ArraySizes)
      Use(Coils)
      Use(Const_3)
      Use(Const_4)

      integer i,j,l
      real zbad,psibad
      character*15 cc1
      character*80 msg(11)

      zbad=z(i)
      psibad=psi(j)
      if (zbad < z1c) then
          cc1='center cell'
      elseif ((zbad < z2c) .and. ncoil == 3) then
          cc1='choke cell'
      elseif ((zbad < z2c) .and. ncoil == 2)then
          cc1='end plug cell'
      elseif (zbad < z3c) then
          cc1='end plug cell'
      endif
 
      write (msg,101) cc1,i,j,zbad,psibad
  101 format('Inconsistent equilibrium attempted.',/5x, 
     * 'Total perpendicular pressure is too large for the vacuum',/, 
     * 5x,'B field in',/, 
     * 15x,a,/5x,'at',/15x,'i= ',i3,/15x,'j= ',i3,/15x,'z= ',e14.6,/, 
     * 15x,'psi= ',e14.6,/5x,'Try again with appropriately reduced ', 
     * 'betas.',/5x,'Sorry for the inconvenience.')
      do l=1,11
          call baswline(STDOUT,msg(l))
      enddo

      call kaboom(0)

      end # errorend







