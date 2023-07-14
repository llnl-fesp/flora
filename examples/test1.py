# flora -p case1 read case1.bas

# Exercise case1 input
from flora import *
from flora.flora_graphics import plots1

utl.probid="Kinetic Stabilizer, Case 1"
import in_iter
import jpostvbb
#read basitervbb.glr.b [NO LONGER USED]
glr.dowriteb=True
glr.glrgen()

#uncomment the next line for plots
#plots1()
