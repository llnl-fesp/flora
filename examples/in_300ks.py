# 
# Namelist2bas version @[#]$RCSfile: main.c,glr.v $ $Revision: 1.20 $
# Report bugs to namelist2bas@alvin.llnl.gov
# 
#
# Beginning of namelist: $noplot
#

from flora import *
#
# End of namelist: $noplot
#
#
# Beginning of namelist: $now1
#
glr.fi1 = 1
glr.fizx = 1
glr.fjrx = -1
glr.fj1 = 0
glr.bias = 0.25
glr.mm = 1
glr.ndiag = 60
glr.fu = 0.5
glr.fv = 0.5
glr.fpsi = 0.5
glr.fz = 0.5
glr.nen = 5
###aname = [ "passing ", "only    ", "flr     ", "stm50b" ] 
###   aname = [ "pfudg1.9", "only    ", "flr     ", "colde-6" ] 
glr.sf6 = 1
glr.sf8 = 1
glr.ex0 = 1
glr.ex1 = 0
glr.nmax = 800
glr.lmax = 4
glr.kzs = 1
#
# End of namelist: $now1
#
#
# Beginning of namelist: $curve
#
glr.betrap = 0
glr.ncenter = 4e+11
glr.nsloshin = 4e11
glr.noshinks=4e11
glr.betpas1 = 0.0000
glr.ppas2 = 0.98
glr.ppas3 = 0.96
glr.dpas1 = 0
glr.d1trap = 100000
glr.betslsh = .25
glr.betslsks= .25
glr.betslse = 0
glr.betcent = .10
glr.betcene = 0
glr.cold = 0.01
glr.dip = 0
glr.z3rel = 0.94
glr.psihrel = 1.6
glr.dpsihrel = 0.03
glr.psi0rel = 0.5
glr.rw1 = 11.99
glr.rw = 12
glr.rp1 = 2.
glr.rp1ks=100.
##glr.p2floor = 0.04
glr.p2floor=0.
glr.betring = 0
glr.psi3rel = 1.05
glr.p2wide = 0.2
glr.long = 1
glr.fring = 0.5
##glr.p2flag = 0.04
glr.p2flag = 0.00
glr.pe10 = 0
glr.p2ewide = 0.1
glr.psi0erel = 0.65
glr.ypot = 1
#
# End of namelist: $curve
#
#
# Beginning of namelist: $field
#
glr.phicen = 0
glr.phiplg = 0

glr.ltrans = 0.1
glr.ztrans = 1.2
### glr.bceng = 1000
glr.bceng=0.
glr.pfudge = 0
##glr.bmx1 = 3200
glr.bmx1 = 4000
glr.bmx2 = 4000
###glr.bmx3 = 6000
glr.bmx3=4000
glr.z1c = 34.2
glr.z2c = 99.5
glr.z3c = 160

glr.ncoil = 2
glr.ass[0]=25.71
glr.ass[1]=25.71
glr.als[0]=12.
glr.als[1]=10.

glr.dt = 3e-9
glr.betslsh=.3;glr.betslsks=3.e-4 
glr.betcent=1.e-7;glr.betcene=1.e-7
glr.rp1=2.
glr.rp1ks=100.
glr.pfudge=0.;glr.pfudgeks=1.9
glr.noshinks=4e5;glr.betslsks=0.;glr.cold=1.e-8;glr.betslsh=.01;
glr.doflat=True;glr.alfrigid=4.;glr.psi0rel=1.;glr.doxrlin=False;glr.doxrtanh=True
glr.doflat=False;glr.psi0rel=0.25;glr.p2wide=.9
glr.doxrlin=False;glr.doxrtanh=True

glr.exbratio=0.
glr.dt=1.e-7;glr.lmax=1;glr.nmax=3000
##check prev run glr.dt=5.e-7 600step with glr.dt=1e-7 3000step
#glr.betcent=1.e-4;glr.betcene=1.e-4
glr.betcent=1.e-128;glr.betcene=1.e-12

##glr.noshinks=4.e11; glr.betslsks=1.e-5   ##old dr8ks
##glr.noshinks=4.e9; glr.betslsks=.7e-3    ##old dr8ks
##glr.noshinks=4.e9; glr.betslsks=18.7668 << equiv old .7e-3 for new normaliz glr.betslsks
##[glr.bmn2/glr.bvac[150]]**2 =    2.68097D+04  use 2.6810e4 as conversion factor

glr.noshinks=4.e10; glr.betslsks=2.6810e0  ## equiv old glr.betslsks=1.e-4
glr.noshinks=4.e10;glr.betslsks=3.e0
glr.alfrigid=999.;glr.sf6=0;glr.sf8=0;glr.lmax=0
glr.dorzzlb=False
glr.cold=1e-6;glr.alfcold=.001;glr.dt=3.e-7
glr.betslsh=1.e-6
glr.pfudgeks=1.99

glr.noshinks=4e7;glr.betslsks=0.
glr.kplotm=2
glr.rp1ks=20.;glr.rpxks=30.
##glr.rp1ks=100.;glr.rpxks=250.
###glr.zmax=300.;ibvmax1=19;ibvmax2=51

###glr.zmax=500.;ibvmax1=12;ibvmax2=31
glr.zmax=300.
##glr.betslsks=3.e-4
glr.betslsks=30.; glr.betslsh=0.1
###glr.betslsks=glr.betslsks*1.12719
###ibksinit=68  no longer used
#
# End of namelist: $field
#

