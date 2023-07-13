# $Id: setup.m,v 1.15 2004/09/21 19:23:51 byers Exp $


      subroutine inputtm
      implicit none

      # Calculates grid quantities needed by sub. grid.

      Use(ArraySizes)
      Use(Btable)
      Use(Coils)
      Use(Const_1)
      Use(Const_1a)
      Use(Const_2)
      Use(Const_3)
      Use(Const_4)
      Use(Const_6)
      Use(Const_7)
      Use(Const_8)
      Use(GFiducials)
      Use(Splines)
      Use(TMInput)
      Use(Fstor)

      integer i, ii, iks
      real dbudz, d2budz2, bum1, bum2
      real z1bum, z2bum, psimax, pmag1, pmag0, pmag2, pmagend
      real bscale, bc1, bcn, zks2(NINPUT), bks2(NINPUT)

c..... boundary conditions are set as follows:
c.....   at z=z0 (i=1), fi1=-1. implies x=0.
c.....                  fi1=1. implies slope=0.
c.....   at z=zmax (i=ix), fizx=-1. implies x=0
c.....                     fizx=1. implies slope=0.
c.....   at psi=psi0 (j=1), fji=-1. implies x=0.
c.....                      fj1=1. implies slope=0.
c.....   at psi=max (j=jx), fjrx=-1. implies x=0.
c.....                       fjrx=1. implies slope=0.
c.....set the proper b.c. at psi=0
      if (fj1 == 0) then
          if (mm == 1)then
              fj1=1.0
          else
              fj1=-1.0
          endif
      endif

c.....transform input z into physical units
      z0=z0rel*zmax
      z1=z1rel*zmax
      z2=z2rel*zmax
      z3=z3rel*zmax
      if (long .ne. 0) z3=z2c
      ztran=ztrans*zmax
      ltran=ltrans*zmax
      bm=bmg/sqrt(4.*Pi)
      bcen=bceng/sqrt(4.*Pi)
      bmx1=bmx1/sqrt(4.*Pi)
      bmx2=bmx2/sqrt(4.*Pi)
      bmx3=bmx3/sqrt(4.*Pi)
      bs(1)=bmx1
      bs(2)=bmx2
      bs(3)=bmx3
      zs(1)=z1c
      zs(2)=z2c
      zs(3)=z3c
c.....calculate bmax, psimax, du, dv, then call grid
      call bvcal(bs,zs,ass,als,z1c,1,ztran,bcen,ltran,bmax,dbudz,d2budz2, 
     *            bum1,bum2,z1bum,z2bum,ncoil)
      psimax=rw**2*bmax*.5
      if (rw1 == 0) rw1=rw
      psiw=rw1**2*bmax*.5
      psi0=psimax*psi0rel
      psih=psimax*psihrel
      dpsih=psimax*dpsihrel
      psi0e=psimax*psi0erel
      psi3=psimax*psi3rel
      v(jx)=psimax
      u(ix)=zmax
      du=u(ix)/(ix-1.5)
      dv=v(jx)/(jx-1.5)
      call grid
c..... calculate bvac, and first and second z derivatives
      call bvcal(bs,zs,ass,als,z,ix,ztran,bcen,ltran,bvac,dbvdz,d2bvdz2, 
     *            bmn1,bmn2,z1min,z2min,ncoil)
      call bvcal(bs,zs,ass,als,0.0,1,ztran,bcen,ltran,bv0,dbudz,d2budz2, 
     *            bum1,bum2,z1bum,z2bum,ncoil)
      call bvcal(bs,zs,ass,als,z3,1,ztran,bcen,ltran,bv3,dbudz,d2budz2, 
     *            bum1,bum2,z1bum,z2bum,ncoil)
      call bvcal(bs,zs,ass,als,z2c,1,ztran,bcen,ltran,bvx2,dbudz,d2budz2, 
     *            bum1,bum2,z1bum,z2bum,ncoil)
      call bvcal(bs,zs,ass,als,z3c,1,ztran,bcen,ltran,bvx3,dbudz,d2budz2, 
     *            bum1,bum2,z1bum,z2bum,ncoil)

c.... Get indices of fiducial points
      call fiducials

c.... Insert Bvac from table, if available
      if (nks > 0) then
          ! Translate and scale B(z) table to match bvac at last mirror peak
          iks = ibvmax(2)
          bscale = bvac(iks)/(bks(1)/sqrt(4*Pi))
          do i = 1,nks
              zks2(i)=zks(i) - zks(1) + z(iks)
              bks2(i)=bscale*bks(i)/sqrt(4*Pi)
          enddo
          ! Adjust nks for table entry consistent with FLORA zmax
          do i = 1,nks
              if (zks2(i) > zmax) break
          enddo
          nks = i - 1
          ! Evaluate end-point derivatives for BCs
          call dydx5(zks2,bks2,nks,bc1,bcn)
          ! Interpolate table data onto bvac(z)
          call interp(zks2,bks2,nks,bc1,bcn,z(iks:ix),bvac(iks:ix), 
     *                 dbvdz(iks:ix),ix-iks+1)
          ! Evaluate B' and B'' numerically everywhere
          call spline(z(2:ix),bvac(2:ix),ix-1,dbvdz(2),bcn,dbvdz(2:ix), 
     *                 d2bvdz2(2:ix),d3bvdz3(2:ix),bint(2:ix))
      endif

cjb subtract variable fraction of bvac(ix) from entire bvac array
cjb   intent is to increase the mirror ratio for the ks species
cjb this shouldnt affect 1st or 2nd derivatives 
cjb        unless they are somehow normalized
      if (dobvsub) then
          bsub=frbvac*bvac(ix)
          bmax=bmax-bsub
          bv0=bv0-bsub
          bv3=bv3-bsub
          bvx2=bvx2-bsub
          bvx3=bvx3-bsub
          do ii=1,ix
              bvac(ii)=bvac(ii)-bsub
          enddo
      endif

c.....transform input constants from relative to physical units
      pmag1=(bmn1)**2*.5
      pmag0=(bv0)**2*.5
      pmag2=bmn2**2*.5
      pmagend=bvac(ix)**2*.5

      z2ct=z2c
c.....special normalization for 2 region problem (ncoil=2)
      if (ncoil == 2) then
          pmag2=pmag1
          bvx3=bvx2
          bmn2=bmn1
          z3c=z2c
          z2ct=z1c
          betpas1=0.
      endif
      psloshin=(betslsh)*pmag2
cjb      poshinks=(betslsks)*pmag2
cjb for ks should use the low b at zmax
      poshinks=(betslsks)*pmagend

cjbfeb12 2004  Note pmagend uses bvac(ix)
cjb  could redefine pmagend using bvac(isubp);
cjb  but we have been running ks w  bvac(ix) in pmagend
cjb  for so long now we probably want to keep it.
cjb  And in fact there is a prescription re betaks using bvac(ix) 

cjb  If we wanted to change so that both plug,ks used pmag2
cjb  then the interpretations of betslsh,betslsks same
cjb  and that could be useful. 
cjb  As it is now w bvac(ix) betslsks not intutitive
 
      pcenter=(betcent)*pmag0
      psloshen=betslse*pmag2
      pcentee=betcene*pmag0
      pring=betring*pmag1
      ppas1=betpas1*pmag0
      ppas2=ppas1*ppas2
      ppas3=ppas1*ppas3
      p1trap=betrap*pmag1
      phice=phicen*psloshin/(nsloshin*ECHARG)
      phipl=phiplg*psloshin/(nsloshin*ECHARG)

      do ii=ibvmax(1),ibvmax(2)
        posh(ii)=psloshin
        nosh(ii)=nsloshin
        bvacnorm(ii)= bvac((ibvmax(1)+ibvmax(2))/2)*bvadj
      enddo
      do ii=ibvmax(2)+1,ix
        posh(ii)=poshinks
        nosh(ii)=noshinks
        bvacnorm(ii)= bvac(isubp)
      enddo

cjb need dummy nonzero value for bvacnorm in cc: rzznewks uses even in cc?
      do ii=1,ibvmax(1)-1
        bvacnorm(ii)=9999999999.
      enddo

cjb posh(ii) set to different constants in plug vs ks
cjb   using an array here only to simplify coding in equil.m
cjb   that now uses same jperpz coding for plug and ks  as repl for freis.
cjb    see subroutine pob2poz also.

      return
      end # inputtm

      subroutine constant
      implicit none       

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_1a)
 
      integer ip,jp,kp1,kp2

      gam1=0.25*(3*bias+1)
      gam2=0.25*(1-bias)
      ip=0.5*(ix-2)
      jp=0.5*(jx-2)
      kp1=ip-1+(jp-2)*(ix-2)
      kp2=jp-1+(jx-2)*(ip-2)
      kplot=0.5*(1+isw)*kp1+0.5*(1-isw)*kp2
      if (kplotm .ne. 0) kplot=kplotm

      return
      end # constant

      subroutine grid
      implicit none

c.....relates physical grid z,psi to computational grid u,v (equally
c.....spaced ). uses input fpsi, fv, fz, fu and azm, apsim .

      Use(ArraySizes)
      Use(Const_1)
      Use(Const_4)

      integer i,j
      real zzp,psip

      xv=log(fpsi)/log(fv)
      xu=log(fz)/log(fu)

      zzp=0.
      psip=0.
      do i=1,ix
          u(i)=u0+du*(i-1.5)
      enddo
      u(1)=-u(1)
      azm=u(ix)**(1.-xu)
      do i=1,ix
          uuz(i)=u(i)**(1.-xu)/(xu*azm)
          z(i)=azm*u(i)**xu
           uuzh(i)=(u(i)+.5*du)**(1.-xu)/(xu*azm)
          dz(i)=z(i)-zzp
          zzp=z(i)
      enddo

      zedge=azm*.5*(u(ix)**xu+u(ix-1)**xu)
      do j=1,jx
          v(j)=v0+(j-1.5)*dv
      enddo
      v(1)=-v(1)
      apsim=v(jx)**(1.-xv)
      do j=1,jx
          vpsi(j)=v(j)**(1.-xv)/(apsim*xv)
          vpsih(j)=(v(j)+.5*dv)**(1.-xv)/(apsim*xv)
          psi(j)=apsim*v(j)**xv
          dpsi(j)=psi(j)-psip
          psip=psi(j)
      enddo

c..... calculates v at plasma edge
      vw=(psiw/apsim)**(1./xv)
      dvin=(vw-v(jx-1))/dv
      dvout=(v(jx)-vw)/dv

      return
      end # grid

      subroutine fiducials # Evaluate bvac fiducial points
      implicit none

      Use(ArraySizes)
      Use(Coils)
      Use(Const_8)
      Use(GFiducials)

      integer i

      nregions=ncoil # Must be identical, for now

      # Region 1
      ibvmin(1)=1
      do i=3,ix
          if (bvac(i) < bvac(i-1)) then
              ibvmax(1)=i-1
              break
          endif
      enddo

      # Region 2
      do i=ibvmax(1)+1,ix
          if (bvac(i) > bvac(i-1)) then
              ibvmin(2)=i-1
              break
          endif
      enddo
      do i=ibvmin(2)+1,ix
          if (bvac(i) < bvac(i-1)) then
              ibvmax(2)=i-1
              break
          endif
      enddo

      if (nregions == 2) return

      # Region 3
      do i=ibvmax(2)+1,ix
          if (bvac(i) > bvac(i-1)) then
              ibvmin(3)=i-1
              break
          endif
      enddo
      do i=ibvmin(3)+1,ix
          if (bvac(i) < bvac(i-1)) then
              ibvmax(3)=i-1
              break
          endif
      enddo

      return
      end # fiducials

      subroutine pob2poz # Map p(B) to p(z)
      implicit none

      Use(ArraySizes)
      Use(Const_2)
      Use(Const_8)
      Use(Fstor)
      Use(GFiducials)

      integer jbt,iter,ib,iplugmid,jbt0
      real dbt,bbsq,bt,btnorm,dbt0,bt0,jperpzbt,jperpzbt0,bvmaxadj

cjb 1st do original ksloop, should be same as previously
cjb  keeping this loop identical to old should at least prevent errors
cjb  to ks region
cjb once plugloop below tested we might combine into one loop

cjb we need separate jperpb etc files for the plug


cjb      do ib=ibvmax(2)+1,ix
cjb modify loop to not use jperpz i> isubp+
cjb  might just want to apply the p2filt in here instead of in equil.m

      do ib=ibvmax(2)+1,isubp+delsubp

        bt=bvac(ib)/bvac(isubp)
        bt0=bt
        if (bt < bvmaxp) then

          do iter=1,100
              jbt0=npob*bt0/bvmaxp
              dbt0=npob*bt0/bvmaxp-jbt0

              jbt=npob*bt/bvmaxp
              dbt=npob*bt/bvmaxp-jbt

cjb  if bt/bvmaxp < 1/npob  then jbt=0 and will pick up junk
cjb  if bt/bvmaxp = npob-del  then jbt= npob-1   so never hit jbt=npob

cjb   jbt= 1 +npob*bt/bvmaxp
cjb           bt ->0  then jbt =1
cjb           bt <bvmaxp/npob then jbt=1
cjb           bt/bvmaxp= 1-eps =  1+ npob -eps*npob  can never ==1+npob

cjbsep1804  that 1+npob seems to have subtle but serious bug
cjb         so revert to old for now, comment out the 1+npob stuff
cjb              jbt0=1+npob*bt0/bvmaxp
cjb              dbt0=1+npob*bt0/bvmaxp-jbt0

cjb              jbt=1+npob*bt/bvmaxp
cjb              dbt=1+npob*bt/bvmaxp-jbt


cjb              jperpz(ib)=jperpb(jbt)+dbt*(jperpb(jbt+1)-jperpb(jbt))

           jperpzbt0=jperpb(jbt0)+dbt0*(jperpb(jbt0+1)-jperpb(jbt0))
           jperpzbt=jperpb(jbt)+dbt*(jperpb(jbt+1)-jperpb(jbt))

              jperpz(ib)=alfks*jperpzbt +(1.-alfks)*jperpzbt0

              bbsq=bvac(ib)**2-2*poshinks*jperpz(ib)
              bt0=bt
              bt=sqrt(bbsq)/bvac(isubp)

              if(ib==isubp) then
                 btisubp(iter)=bt
                 jpzsubp(iter)=jperpz(isubp)
              endif
              if(ib==isubp-1) then
                 btisubpm1(iter)=bt
                 jpzsubpm1(iter)=jperpz(isubp-1)
              endif

              if(ib==isubp-2) then
                 btisubpm2(iter)=bt
                 jpzsubpm2(iter)=jperpz(isubp-2)
              endif

              if(ib==isubp-3) then
                 btisubpm3(iter)=bt
                 jpzsubpm3(iter)=jperpz(isubp-3)
              endif
            
          enddo
          jparz(ib)=jparb(jbt)+dbt*(jparb(jbt+1)-jparb(jbt))
          jdfz(ib)=jdfb(jbt)+dbt*(jdfb(jbt+1)-jdfb(jbt))
          jd2fz(ib)=jd2fb(jbt)+dbt*(jd2fb(jbt+1)-jd2fb(jbt))
          ddbprpz(ib)=ddbprpb(jbt)+dbt*(ddbprpb(jbt+1)-ddbprpb(jbt))
          ddbparz(ib)=ddbparb(jbt)+dbt*(ddbparb(jbt+1)-ddbparb(jbt))
        endif
      enddo


cjb 2nd,repeat for plug, needs bvadj,bvmaxpg  in glr.v  and in input
cjb uses psloshin instead of poshinks
cjb  note we are not using the posh(i) array here ir in ks loop above,
cjb       which we might do if combine into one loop
cjb test plug first at lower betslsh ~.1
cjb   then do betslsh=.5  with added underrelaxtion for iter

cjb not sure if this simple linear interpolation will work in 
cjb    nonmonotonic bplug(z)

cjb need a separate set of 'b' arrays for plug  jperpbpg, etc


cjb feb14 2004  impose underrelaxation for iteration

cjb  model after mma d2fplug.mnew :
cjb alf=.5
cjb Do[bt79=Table[(bvii[[i]]^2 -UnitStep[bvmax-bvii[[i]]]*2*poshinks/bvac140^2*
cjb   (alf*jperp[btot[[i]]] +(1-alf)*jperp[btoto[[i]]]))^.5,{i,43,120}];
cjb    btoto=btot;
cjb    btot = Join[bv78,bt79];  (* has bvii vals in 1-78  ie not zeros *)
cjb    btoti65[[iter]]=btot[[65]];
cjb    btoti82[[iter]]=btot[[82]]
cjb   (* test early iter ListPlot[bt79,PlotLabel->"bt79 43,120 in iterloop"]*)
cjb         ,{iter,1,200}]
cjbplbtoti65=ListPlot[btoti65, PlotLabel->"btot[65] vs iter",PlotJoined->True]

      if(doplugpob) then

        iplugmid=(ibvmax(1)+1+ibvmax(2))/2
        btnorm=bvac(iplugmid)*bvadj

        bvmaxadj=bvmaxpg/bvadj

        do ib=ibvmax(1)+1,ibvmax(2)
          bt=bvac(ib)/btnorm
          bt0=bt

cjb            if (bt < bvmaxpg) then
            if (bt < bvmaxadj) then


cjb alfpg initilzd in glr.v  /0.5/
              do iter=1,100
       
cjb              jbt0=npob*bt0/bvmaxpg
cjb              dbt0=npob*bt0/bvmaxpg-jbt0
cjb              jbt=npob*bt/bvmaxpg
cjb              dbt=npob*bt/bvmaxpg-jbt


              jbt0=npob*bt0/bvmaxadj
              dbt0=npob*bt0/bvmaxadj-jbt0
              jbt=npob*bt/bvmaxadj
              dbt=npob*bt/bvmaxadj-jbt


cjb           jperpz(ib)=jperpbpg(jbt)+dbt*(jperpbpg(jbt+1)-jperpbpg(jbt))

           jperpzbt0=jperpbpg(jbt0)+dbt0*(jperpbpg(jbt0+1)-jperpbpg(jbt0))
           jperpzbt=jperpbpg(jbt)+dbt*(jperpbpg(jbt+1)-jperpbpg(jbt))

              jperpz(ib)=alfpg*jperpzbt +(1.-alfpg)*jperpzbt0

cjb update bt0 only just before bt is to be updated
              bt0=bt
              bbsq=bvac(ib)**2-2*psloshin*jperpz(ib)
              bt=sqrt(bbsq)/btnorm

              if(ib==iplugmid) bti82(iter)=bt
            enddo
        
            jparz(ib)=jparbpg(jbt)+dbt*(jparbpg(jbt+1)-jparbpg(jbt))
            jdfz(ib)=jdfbpg(jbt)+dbt*(jdfbpg(jbt+1)-jdfbpg(jbt))
            jd2fz(ib)=jd2fbpg(jbt)+dbt*(jd2fbpg(jbt+1)-jd2fbpg(jbt))
            ddbprpz(ib)=ddbprpbpg(jbt)+dbt*(ddbprpbpg(jbt+1)-ddbprpbpg(jbt))
            ddbparz(ib)=ddbparbpg(jbt)+dbt*(ddbparbpg(jbt+1)-ddbparbpg(jbt))
          endif
        enddo
      endif
c234567890123456789012345678901234567890123456789012345678901234567890

cjb maybe we could test this plug iteration loop in mma?
cjb  ie at least to see that the method doesnt glitch near the minimum
cjb  the d2f.m codes still only have bvac for z300
cjb   but that shouldnt prevent testing the idea
cjb  I have already tested the mma plugiteration, and that works;
cjb    what is in question is this cruder linear interpolation
cjb    and especially its possible dependence on monotonic bvac(i) as in ks


      return
      end # pob2poz

      subroutine storage # Allocate and clear storage
      implicit none

      Use(ArraySizes)
      Use(Matrix)

      # Verify array sizes
      ksimp=2*(ksimp/2)+1 # need an odd number for Simpson

      # Determine secondary array sizes
      kxx=(ix-2)*(jx-2)
      kbw=min(ix-1,jx-1)
      lda=3*kbw+1

      # Allocate global dynamic arrays
      call gchange("Const_4",0)
      call gchange("Const_5",0)
      call gchange("Const_7",0)
      call gchange("Const_8",0)
      call gchange("GEnergy",0)
      call gchange("Fstor",0)
      call gchange("Glrexec",0)
      call gchange("Matrix",0)

      # Initialize arrays
      call myzero(kxx*9,0.0,a1(1,1),1)
      call myzero(lda*kxx,0.0,abar(1,1),1)

      return
      end # storage
