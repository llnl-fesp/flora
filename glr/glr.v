glr

#$$$$$$$$$$$$$$$$$$$$$$$$$$  C O N V E N T I O N S  $$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Attention...
#       ?             : Item needs attention by developers.
#
# Attributes...
#       input         : Identifies FLORA input quantities.
#
# Unit specifications (in addition to physical units)...
#       dimensionless : Dimensionless quantity
#       flag          : Integral quantity or logical variable
#       none          : Default (provided by Basis if not in .v file)
#
# Comments...
#       Comment strings should not exceed 70 characters, otherwise they will be
#       truncated by Basis. Use multiple lines for long comments.
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#$$$$$$$$$$$$$$$$$$$$$$$$$$  P A R A M E T E R S  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
{

  # STATIC STORAGE LIMITS
    NCOIL=3		# Maximum number of coils (see NREGIONS).
    NREGIONS=NCOIL	# For now, these must be the same
    NINPUT=1000		# Maximum length of input arrays
}
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  R O U T I N E S  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
***** GCONST:
# CONSTANTS
ECHARG real /4.8e-10/
Pi     real /3.1415926535/
    # Ion charge (statcoulomb).
***** Flora_routines:

glrgen				subroutine
amat								subroutine
	# ?
amax2(a:real,im:integer,jm:integer,amax:real,
	imx:integer,jmx:integer)				subroutine
	# Given 2D array a(im,jm), amax2 finds its maximum (amax) and
	# corresponding indices (imx,jmx).
bcccal(bs:real,zs:real,ass:real,als:real,z:real,n:integer,
	bv:real,bvp:real,bvpp:real,b4:real,b5:real,z4:real,
	z5:real,nc:integer)					subroutine
	# ?
bvcal(bs:real,zs:real,ass:real,als:real,z:real,n:integer,
	zcc:real,bcc:real,znorm:real,bv:real,bvp:real,
	bvpp:real,b4:real,b5:real,z4:real,z5:real,
	nc:integer)						subroutine
	# ?
comat								subroutine
	# ?
constant							subroutine
	# ?
diagno								subroutine
	# ?
ddpsi(f1:real,f2:real,i:integer)				subroutine
	# ?
ddz(f1:real,f2:real)						subroutine
	# ?
dydx5(x:real,y:real,n:integer,yp1:real,ypn:real)		subroutine
	# Evaluate end-point derivatives
energy								subroutine
	# ?
equiltm								subroutine
	# ?
errorend(i:integer,j:integer)					subroutine
	# ?
f1to11								subroutine
	# ?
grid								subroutine
	# ?
initial								subroutine
	# ?
inputtm								subroutine
	# ?

interp(x:real,y:real,n:integer,bc1:real,bcn:real,x2:real,y2:real,y2p:real,m:integer) subroutine
	# Interpolates from y(x(1:n)) onto y2(x2(1:m)) using cubic splines.
	# The boundary conditions (see subroutine spline) are specified with
	# bc1 and bcn. The 1st derivative of y2, y2p(1:m), is also evaluated.

mymove(a:real,b:real,len:integer)				subroutine
	# ?
rightvec							subroutine
	# ?
simps(fin:real,fout:real,knn:integer,df:real)			subroutine
	# ?

spline(x:real,y:real,n:integer,bc1:real,bcn:real,yp:real,ypp:real,yppp:real,yint:real) subroutine
	# Evaluate derivatives and integrals of y(x), using cubic splines.
	# The input is contained in x(1:n) and y(1:n) and the boundary
	# conditions are specified with bc1 (1st point) and bcn (last point):
	#   bc1, bcn == -666 periodic;
	#   bc1, bcn ==  777 not-a-knot;
	#   bc1, bcn == -111 extrapolates for 1st deriv;
	#   bc1, bcn ==  555 natural 2nd deriv=0;
	#   bc1, bcn == -777 3rd deriv=0;
	# otherwise bc1 and bcn contain value of 1st derivative.
	# The output consists of the three derivatives: yp, ypp, yppp and
	# the integral yint (all of length 1:n).
	# The running integral is not computed if yint(1)=-999 on input.
	# Except for periodic splines, scalar yi=-999 will work.

tmcon2								subroutine
	# ?
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
***** ArraySizes:
ix		integer	/150/				+input
	# Number of axial grid points.
jx		integer /50/				+input
	# Number of radial (psi) grid points.
npob		integer /200/				+input
	# Number of p(B) points.
ksimp		integer /23/				+input
	# Number of points for Simpson's Rule integration.
nmax		integer	/5/		[dimensionless]	+input	$ WHY SO LOW?
	# Upper limit of time steps for problem.
nps		integer	/100/				+input
	# ?
neng		integer	/500/				+input
	# ?
kxx		integer /0/
	# Secondary array size, kxx=(ix-2)*(jx-2).
kbw		integer /0/
	# Secondary array size, kbw=min(ix-1,jx-1).
lda		integer /0/
	# Secondary array size, lda=3*kbw+1.


***** Fstor:

f1(ix,jx)	_real
f2(ix,jx)	_real
f3(ix,jx)	_real
f4(ix,jx)	_real
f5(ix,jx)	_real
f7(ix,jx)	_real
g1(ix,jx)	_real
g2(ix,jx)	_real
g3(ix,jx)	_real
g4(ix,jx)	_real
swg1		real	/1.0/		[dimensionless]	+input
	# Arbitrary scaling factor on curvature-drive term.
swg2		real	/1.0/		[dimensionless]	+input
	# Arbitrary scaling factor on some of the line-bending terms.
swg3		real	/1.0/		[dimensionless]	+input
	# Arbitrary scaling factor on some of the line-bending terms.
swg4		real	/1.0/		[dimensionless]	+input
	# Arbitrary scaling factor on some of the line-bending terms.

b(ix,jx)	_real
rho(ix,jx)	_real
qub(ix,jx)	_real
r(ix,jx)	_real
yyy(ix,jx)	_real
xxx(ix,jx)	_real

xioo(kxx)	_real
xio(kxx)	_real
xiol(kxx)	_real
xroo(kxx)	_real
xro(kxx)	_real
xrol(kxx)	_real

omstri1(jx)	_real
omgbi1(jx)	_real
omstri34(jx)	_real $ Specific to ix=150 ?
omstri66(jx)	_real $ Specific to ix=150 ?
omstr148(jx)	_real $ Specific to ix=150 ?
omebj1(ix)	_real
omstrj1(ix)	_real

exbratio	real	/-1.0/		[?]		+input
	# omegexb = exbratio*omegstr.
alfrigid	real	/1.0/		[?]		+input
	# larger for faster decay rigidrotor dpdpsi/p
###ibksinit	integer	/79/		[none]		+input
	# Index in z marking start of KS pressure
isubp		integer	/140/		[none]		+input
	# Index in z marking end of KS pressure
alfsubp		real	/0.5/
dosubp		logical	/.true./
doflat		logical	/.false./	[flag]		+input
	# ?
doxrlin		logical	/.false./	[flag]		+input
	# ?
doxrtanh	logical	/.true./	[flag]		+input
dowriteb	logical	/.false./
dorzzlb		logical	/.false./	[flag]		+input
	# ?
doplugpob       logical	/.true./  	[flag]		+input
        # if true iterate in plug for bt,jperpz  in pob2poz    

***** Matrix:
ipvt(kxx)	_integer
a1(kxx,9)	_real
a2(kxx,9)	_real
a3(kxx,9)	_real
b1(kxx,3)	_real
bc(jx)		_real
rhs1(kxx)	_real
rhs2(kxx)	_real
abar(lda,kxx)	_real


***** Const_1a:
iscreen		integer	/1/		[dimensionless]	+input
	# Tells how often to write timestep+ to screen.
mm		integer	/4/		[dimensionless] +input
	# Azimuthal mode number.
lmax		integer	/2/		[dimensionless]	+input
	# Iteration limit.
isw		integer
ndiag		integer	/100/		[dimensionless]	+input	$ USED?
	# Number of time steps between diagnostic plots.
kplot		integer
npm		integer
kplotm		integer	/0/		[dimensionless]	+input	$ USED?
	# Index of spatial location of time history plots, if ==0, center of
	# region is used.
kzs		integer	/1/		[flag]		+input
	# Flute mode initialization, ==1 & ex0==1 & ex1==1 sets initial
	# condition r B xro=0 and and r B xio=0.
n		integer
	# Time loop counter.
ltt		integer
nen		integer	/1/		[dimensionless]	+input
	# Number of time steps between energy checks.
lee		integer


***** Const_1:

startime	real [sec. CPU]  # Initial value from  'call second' (Cray).
endtime		real [sec. CPU]  # Final value from 'call second' (Cray).
totalcpu	real [sec. CPU]  # Total timeloop CPU from 'call second' (Cray/Sun)
t		real

gam1		real
gam2		real
fac1		real
fac2		real
bias		real	/0.5/		[dimensionless]	+input
	# Time centering parameter, ==0 fully centered, ==1 fully forward bias.

du		real	/1.0/
dv		real	/1.0/
dt		real	/1.0e-5/	[s]		+input
	# Time step.
ex0		real	/0.0/		[dimensionless]	+input
	# Initial perturbation coefficient, ==1 random perturbation.

ex1		real	/0.0/		[dimensionless]	+input
	# Initial perturbation coefficient, ==1 for cosine cosine distribution.
fi1		real	/1.0/		[flag]		+input $ WHY REAL?
	# Minimum z boundary condition, ==1 for 0 slope, ==-1 for 0 value.
fizx		real	/1.0/		[flag]		+input $ WHY REAL?
	# Maximum z boundary condition, ==1 for 0 slope, ==-1 for 0 value.
fj1		real	/1.0/		[flag]		+input $ WHY REAL?
	# Minimum psi boundary condition, ==1 for 0 slope, ==-1 for 0 value;
	# For fj1==0, mm==1 results in 0 slope, mm>=2 results in 0 value.
fjrx		real	/1.0/		[flag]		+input $ WHY REAL?
	# Maximum psi boundary condition, ==1 for 0 slope, ==-1 for 0 value.
fpsi		real	/0.5/		[dimensionless]	+input
	# Grid stretching parameter [0,1].
fu		real	/0.5/		[dimensionless]	+input
	# Grid stretching parameter [0,1].
fv		real	/0.5/		[dimensionless]	+input
	# Grid stretching parameter [0,1].
fz		real	/0.5/		[dimensionless]	+input
	# Grid stretching parameter [0,1].g
azm		real	/1.0/
apsim		real	/1.0/

u0		real
v0		real
amass		real	/3.34e-24/
omegexb		real
sf6		real	/1.0/		[dimensionless]	+input
	# Arbitrary scaling factor on gyroscopic FLR term xxx.
sf8		real	/1.0/		[dimensionless]	+input
	# Arbitrary scaling factor on quasi-elastic FLR term yyy.
zedge		real

xu		real
xv		real
vw		real
psiw		real
dvin		real
dvout		real



***** Const_2:
psi0rel		real	/0.50/		[dimensionless]	+input
	# Psi value relative to psimax at which ion radial pressure is half the
	# maximum.
z1rel		real	/0.5/
z2rel		real	/1.0/
z3rel		real	/0.625/		[?]		+input
	# Outer axial position where hot electrons go to zero (pring=0).
z0rel		real	/0.75/
nslosh		real
	# Peak sloshing ion number density.
nsloshks	real

bmg		real	/1.0e4/
ncenter		real	/1.0e13/	[1/cm^3]	+input
	# Center cell density.
pslosh		real
psloshks	real

pcenter		real
ztrans		real	/0.4/		[dimensionless]	+input
	# Z position normalized to zmax to middle of center cell field
	# transition region.
ltrans		real	/0.05/		[dimensionless]	+input
	# Scale length, normalized to zmax, for center cell vacuum field
	# transition region.g
bm		real
ltran		real
ztran		real
rp1		real	/2.0/		[dimensionless]	+input
	# Mirror ratio of the inner component of the sloshing profile
	# (UCID-20400, Sect. 6). 
rp1ks		real	/100.0/		[dimensionless]	+input
	# Mirror ratio of inner comp of KS sloshing profile.

bcen		real
pring		real
epsp		real	/1.0e-6/	[dimensionless]	+input
	# Approximately the minimum pressure (normalized to 1), below which B
	# is calculated by an expansion (UCID-20400, Sect. 3).
phicen		real	/0.1/		[?]		+input
	# Maximum center cell potential relative to plug ion tperp.
phiplg		real	/0.1/		[?]		+input
	# Maximum plug potential relative to plug ion tperp.
xpot		real	/1.0/		[dimensionless]	+input
	# Exponent coefficient in center cell potential.
ypot		real	/2.0/		[dimensionless]	+input
	# Power of polynomial in potential psi dependence.
wpot		real	/2.0/		[dimensionless]	+input
	# Exponent coefficient in end plug cell potential.
pfudge		real	/0.0/		[dimensionless]	+input
	# Parameter on sloshing shape.
pfudgeks	real	/0.0/		[?]		+input
	# ?
rpx		real
rpxks		real	/0.0/		[?]		+input
	# ?

phice		real
phipl		real
betslsh		real	/0.25/		[dimensionless]		+input
	# Peak sloshing ion perpendicular beta, wrt bmn2.
betslsks	real	/0.0/		[dimensionless]	+input
	# Peak KS sloshing ion perpendicular beta, wrt bmn2.

betcent		real	/0.10/		[dimensionless]		+input
	# Peak center cell ion perpendicular beta, ert bv0.
z0		real
z1		real
z2		real
z3		real
psi0		real
betcene		real	/0.1/		[dimensionless]		+input
	# Peak center cell electron perpendicular beta, wrt bv0.
betslse		real	/0.1/		[dimensionless]		+input
	# Peak warn sloshing electron perpendicular beta, wrt bmn2.
psloshe		real
pcentee		real
bmax		real
alsi		real
alsiks		real
bm1		real
bm1ks		real
bvoutks		real

psls1		real
cold		real	/1.0e-5/	[dimensionless]	+input
	# Global density minimum as fraction of ncenter.
alfcold		real	/1.0/		[?]		+input
	# exp(-alfcold(JRX-j)) psivariable model  alf 1.->.1 
p2wide		real	/0.1/		[dimensionless]	+input
	# Parameter inversely proportional to ramp width of p2t in the ion
	# radial pressure profile (UCID-20400, Sect. 8).
                       #  ion radial pressure profile
psi3rel		real	/1.05/		[dimensionless]	+input
	# Psi value relative to psimax where electric field=0.
psi3		real
p1max		real
bv0		real
bv3		real
bceng		real	/1.0e3/		[G]		+input
	# Center cell magnetic field.
psloshin	real
poshinks	real

psloshen	real
nsloshin	real	/2.0e13/	[1/cm^3]	+input
	# Peak plug cell density.
noshinks	real	/0.0/		[?]		+input
	# ?

p3a		real
p3b		real
p3c		real
p3d		real
psim		real
pe10		real	/0.0/		[?]		+input
	# Coefficient of hot electron radial pressure profile (UCID-20400,
	# Sect. 8)
ae1		real
be1		real
ce1		real
de1		real
psi0erel	real	/0.5/		[?]		+input
	# Psi value relative to psimax at which hot electron radial pressure is
	# half the maximum (eringpperp(pe2)=halfmax).
psi0e		real
psime		real
p2ewide		real	/0.2/		[?]		+input
	# Parameter inversely proportional to ramp width of p2e in the hot
	# electron radial pressure profile (UCID-20400, Sect. 8).
wp2e		real
p1floor		real	/0.0/		[?]		+input
	# Value to which p1 is set if p2 <= p2flag.
p2flag		real	/1.0e-3/	[?]		+input
	# Value of p2 at which p1 and p2 are set constant.
p2floor		real	/0.0/		[?]		+input
	# Value to which p2 is set if p2 <= p2flag.
fring		real	/0.9/		[?]		+input
	# Used for elongated hot electrons (UCID-20400, Sect. A.3).
dip		real	/1.0/		[?]		+input
	# Parameter in psi pressure profile (UCID-20400, Sect. 8).
psistr		real
psislp		real	/0.01/		[?]		+input
	# Coefficient of pe1 (UCID-20400, Sect. 8).
psihrel		real	/0.0/		[?]		+input
	# Psi value relative to psimax of halo.
psih		real
dpsihrel	real	/0.0/		[?]		+input
	# Psi width relative to psimax of transition to halo region.
dpsih		real


***** Coils: # Solenoid coil parameters
ncoil		integer	/2/		[none]		+input
	# Number of solenoids (also regions).
ass(NCOIL)	real	/NCOIL*0.0/	[cm]		+input
	limited(ncoil)
	# Solenoid radius.
als(NCOIL)	real	/NCOIL*0.0/	[cm]		+input
	limited(ncoil)
	# Solenoid coil length.
zs(NCOIL)	real
	limited(ncoil)
bs(NCOIL)	real
	limited(ncoil)
z1c		real	/0.0/		[cm]		+input
	# Axial position of choke coil.
z2c		real	/0.0/		[cm]		+input
	# Axial position of end plug inboard mirror coil.
z3c		real	/0.0/		[cm]		+input
	# Axial position of end plug outboard mirror coil.


***** GFiducials:
nregions	integer	/0/		[none]
	# Number of regions
ibvmin(NREGIONS)	integer		[none]
	limited(nregions)
	# Indices of B minima in z
ibvmax(NREGIONS)	integer		[none]
	limited(nregions)
	# Indices of B maxima in z


***** Const_3:
at0		real
bt0		real
ct0		real
cp0		real
ap1		real
ap2		real
ap3		real
at1		real
at2		real
at3		real
bp1		real
bp2		real
bp3		real
bt1		real
bt2		real
bt3		real
cp1		real
cp2		real
cp3		real
ct1		real
ct2		real
ct3		real
ct2ks		real
	# Mimic ct2 for outer kinetic stabilizer sloshing .
ct3ks		real
	# Mimic ct3  for outer kinetic stabilizer sloshing.
bmx1		real	/0.0/		[G]		+input
	# Magnetic field at choke coil solenoid (UCID-20400, Sect. 3).
bmx2		real	/0.0/		[G]		+input
	# Magnetic field at inboard end plug solenoid (UCID-20400, Sect. 3).
bmx3		real	/0.0/		[G]		+input
	# Magnetic field at outer end plug solenoid (UCID-20400, Sect. 3).
bmn1		real
bmn2		real
ppas1		real
ppas2		real	/0.0/		[?]		+input
	# Minimum passing pressure in choke cell, fraction of ppas1
	# (UCID-20400, Sect. 7).
ppas3		real	/0.0/		[?]		+input
	# Maximum passing pressure at inboard mirror of end plug cell.
p1trap		real
z1min		real
z2min		real
dpas1		real	/0.0/		[1/cm^3]	+input
	# Center cell passing density.
d1trap		real	/0.0/		[1/cm^3]	+input
	# Peak density in choke cell.
betrap		real	/0.0/		[none]		+input
	# Peak choke cell perpendicular beta.
betpas1		real	/0.0/		[none]		+input
	# Center cell passing beta.
bvx2		real
bvx3		real
z2ct		real

***** Const_4:
psi(jx)		_real
z(ix)		_real
u(ix)		_real
v(jx)		_real
dpsi(jx)	_real
dz(ix)		_real
vpsi(jx)	_real
uuz(ix)		_real
vpsih(jx)	_real
uuzh(ix)	_real


***** Const_5:
xrtime(nmax)	_real
xrspz(ix,nps)	_real
xrsppsi(jx,2*nps)	_real
time(nmax)	_real
xflute(ix,nps)	_real
enpot(neng)	_real
tenergy(neng)	_real
enkin(neng)	_real
timengy(neng)	_real
tenrel(neng)	_real
enbend(neng)	_real
encurve(neng)	_real
enflr(neng)	_real


***** Const_6:
lb		real	/1.0e9/
	# With doflat==true for fixed curv rzz =-r/lb^2
rw		real	/0.0/		[cm]		+input
	# Wall radius.
cee		real	/3.0e10/


***** Const_7:
h12(ix)		_real
hzt0(ix)	_real
h34(ix)		_real
abp(ix)		_real
bbp(ix)		_real
cbp(ix)		_real
abf(ix)		_real
bbf(ix)		_real
cbf(ix)		_real
hp3(jx)		_real
htrans(ix)	_real
abq(ix)		_real
bbq(ix)		_real
cbq(ix)		_real
ebp(ix)		_real
fbp(ix)		_real
gbp(ix)		_real
betring		real	/0.0/		[none]		+input
	# Peak hot electron perpendicular beta, wrt bmn1.
hp0(jx)		_real
hpm(jx)		_real
hpme(jx)	_real
hflr(jx)	_real
hzp0(ix)	_real
hzp1(ix)	_real
hzp2(ix)	_real
hzp3(ix)	_real
hzt1(ix)	_real
hzt2(ix)	_real
hzt3(ix)	_real
hzt2ks(ix)	_real
hzt3ks(ix)	_real
dterjb(ix)	_real


***** Const_8:
bvac(ix)	_real
dbvdz(ix)	_real
d2bvdz2(ix)	_real
d3bvdz3(ix)	_real
bint(ix)	_real
dp1dpsi(jx)	_real
p1(jx)		_real
p1k(ksimp,jx)	_real
hpk0(ksimp)	_real
deli1(ksimp)	_real
deli2(ksimp)	_real
deli3(ksimp)	_real
deli4(ksimp)	_real
rzz(ix,jx)	_real
rzz1a(ix)	_real
rzz1b(ix)	_real
rzz2(ix)	_real
rzz3(ix)	_real
rzz1anew(ix)	_real
rzz1bnew(ix)	_real
rzz2new(ix)	_real
rzz3new(ix)	_real
rzznew(ix)	_real
rzzold(ix)	_real

rd2diff(ix)	_real
rnew(ix)	_real
rnewd2(ix)	_real
rold(ix)	_real
roldd2(ix)	_real

icapi2(ix)	_real
icapi3(ix)	_real
icapi5(ix)	_real
icap2new(ix)	_real
icap3new(ix)	_real
icap5new(ix)	_real

dbdpsi(ix,jx)	_real
phi1(ix)	_real
phi2(ix)	_real
pperp(ix,jx)	_real
ppar(ix,jx)	_real
jperpz(ix)	_real
jperpb(NINPUT)	real
	limited(npob)
jparb(NINPUT)	real
	limited(npob)
jdfb(NINPUT)	real
	limited(npob)
jd2fb(NINPUT)	real
	limited(npob)
ddbparb(NINPUT)	real
	limited(npob)
ddbprpb(NINPUT)	real
	limited(npob)

jperpbpg(NINPUT)	real
	limited(npob)
jparbpg(NINPUT)	real
	limited(npob)
jdfbpg(NINPUT)	real
	limited(npob)
jd2fbpg(NINPUT)	real
	limited(npob)
ddbparbpg(NINPUT)	real
	limited(npob)
ddbprpbpg(NINPUT)	real
	limited(npob)


jparz(ix)	_real
jd2fz(ix)	_real
jdfz(ix)	_real
jbtotz(ix)	_real $ ? This variable is not used.
jmbselfz(ix)	_real $ ? This variable is not used.
ddbprpz(ix)	_real
ddbparz(ix)	_real
bvmaxp           real	/20.0/		[dimensionless]	+input
	# maxb for jperpz (previously hardwired to 20)
bvmaxpg          real /2.10/ $ maxb for jperpz for plug
bvadj            real /1.0/ $ < 1, selfb adjustment for btnorm in plug

posh(ix)        _real $poshinks i>ibvmax(2); psloshin in plug; notusedincc?
nosh(ix)        _real $noshinks i>ibvmax(2); nsloshin in plug; notusedincc?
bvacnorm(ix)    _real $bvac(isubp) i>ibvmax(2); bvac(iplugmid) in plug; dummyincc?
bti82(100)      _real $diagnostic for iteration in plug
alfpg            real /0.5/ $underrelax for iteration in plug
btisubp(100)    _real $diagnostic for iteration in ks
btisubpm1(100)  _real $diagnostic for iteration in ks
btisubpm2(100)  _real $diagnostic for iteration in ks
btisubpm3(100)  _real $diagnostic for iteration in ks
jpzsubp(100)    _real $diagnostic for iteration in ks
jpzsubpm1(100)  _real $diagnostic for iteration in ks
jpzsubpm2(100)  _real $diagnostic for iteration in ks
jpzsubpm3(100)  _real $diagnostic for iteration in ks
alfks            real /0.5/ $underrelax for iteration in ks
delsubp          integer /0/ $ks pob2poz loopindex = isubp+delsubp

bj2old(ix)	_real
pperpold(ix)	_real
pperpold25(ix)  _real
pparold(ix)	_real
betaold(ix)	_real
qubold(ix)	_real
qubvold(ix)	_real
dbdpsiold(ix)	_real

dflute3(ix)	_real
qubv(ix,jx)	_real
p2(jx)	_real
dp2dpsi(jx)	_real
dflute1(ix)	_real
dflute2(ix)	_real
flute1(jx)	_real
flute2(jx)	_real
flute3(jx)	_real
dpost1(ix)	_real $ ? This variable is not used.
dpost2(ix)	_real $ ? This variable is not used.
post1		real $ ? This variable is not used.
post2		real $ ? This variable is not used.
dpost1r(ix)	_real
dpost2r(ix)	_real
dpostjb(ix)	_real
dpostjbo(ix)	_real
post1r		real
post2r		real
postjb		real
post1sum	real
post2sum	real
postjbsum	real
postjbsumo	real

post1cc		real
post2cc		real
postjbcc	real
post1pg		real
post2pg		real
postjbpg	real
post1ks		real $ ? This variable is not used.
post2ks		real $ ? This variable is not used.
postjbks	real $ ? This variable is not used.
p1kssubp	real
p1ksmax		real
p2ksmax		real
pjbksmax	real
pjbksmaxo	real

pnorm		real
ptot(ix)	_real
beta(ix)	_real
pperpinj	real
btksinj		real
pparinj		real
pperpinjo	real
btksinjo        real

betanorm(ix)	_real
p2k(ksimp,jx)	_real
deli5(ksimp)	_real
del1new(ksimp)	_real
del2new(ksimp)	_real
del3new(ksimp)	_real
del5new(ksimp)	_real

pperps(ix,jx)	_real
errprp(ix,jx)	_real
errprl(ix-2,jx-1)	_real
dqubdb(ix,jx)	_real
xxxfreq(jx)	_real
pperpe(ix,jx)	_real
epsi(ix,jx)	_real
omeg1wkb(jx)	_real
omeg2wkb(jx)	_real
gamwkb(jx)	_real
dflute4(ix)	_real
rhoave(jx)	_real
xxxave(jx)	_real
yyyave(jx)	_real
p2t(jx)		_real
dp2dpst(jx)	_real
p3(jx)		_real
dp3dpsi(jx)	_real
hpkm(ksimp)	_real
hpkme(ksimp)	_real
droave(jx)	_real
droterm(ix)	_real
ering(ix,jx)	_real
grow		real
growmax		real
xfreqmax	real
rzzrmax		real


***** TMInput:
long		integer	/1/		[dimensionless]	+input
	# Switch which sets hot electron z-length as elongated (long==1) or
	# regular (long==0).
rw1		real	/0.0/		[cm]		+input
	# Equal, or slightly less (within one grid cell), to rw, if zero, gets
	# set to rw.
zmax		real	/0.0/		[cm]		+input
	# Maximum z of the domain.
info		integer
	# Returned by dgbtrf and dgbtrs.
bsub		real
	# Value subtracted from bvac.
frbvac		real	/0.0/
	# Allowed range is 0. to 1.-delta.
dobvsub		logical	/.false./
	# If true, subtract frbvac*bvac(ix).


***** GEnergy: # Really want local, but then can't allocate dynamically
drxi(ix,jx)     _real
drxio(ix,jx)    _real
drxr(ix,jx)     _real
drxro(ix,jx)    _real
dxiz2(ix)       _real
dxrz2(ix)       _real
dum1(ix)        _real
dum2(ix)        _real
psuma(jx,12)    _real
psumb(jx,4)     _real
quad(12)        real
quad1(4)        real
rxi(jx)         _real
rxr(jx)         _real
xig(ix)         _real
xigo(ix)        _real
xrg(ix)         _real
xrgo(ix)        _real



***** Glrexec: # Really want local, but then can't allocate dynamically
rhs1b(kxx,1)	_real
rhs2b(kxx,1)	_real



***** Splines: # Storage for subroutine interp

mspin		integer
mspout		integer
lsp(mspout)	_integer
dsp(mspout)	_real
csp1(mspin)	_real
csp2(mspin)	_real
csp3(mspin)	_real
csp4(mspin)	_real

***** Btable: # Storage for imported B in KS region
nks		integer		/0/			+input
	# Number of points in imported B in KS region
zks(NINPUT)	real			[cm]		+input
	# Z positions for imported B(z) in KS region
	limited(nks)
bks(NINPUT)	real			[dimensionless]	+input
	# B at zks for imported B(z) in KS region
	limited(nks)
zks2(NINPUT)	real			[cm]
bks2(NINPUT)	real			[dimensionless]
