try:
   import IPython
   from IPython.terminal.prompts import Prompts,Token
   from IPython.terminal.embed import InteractiveShellEmbed
except:
   pass

try:
   from traitlets.config.loader import Config
except:
   pass

import sys,os,__main__

import numpy as np
from numpy import array,tanh,exp,arange
ArrayType = np.ndarray

if sys.hexversion >= 0x03000000:
    # --- With Python3, the so files of each Fortran package are imported
    # --- separately. The dlopen flag needs to be set so that cross references
    # --- among the packages can be satisfied.
    sys.setdlopenflags(os.RTLD_LAZY | os.RTLD_GLOBAL)

from . import floraC
from .floraC import *

#from Forthon import *

if sys.hexversion >= 0x03000000:

    from .glrpy import glr
    from .utlpy import utl
else:
    from glrpy import glr
    from utlpy import utl

import time
import os.path
import __main__
# import all of the neccesary packages



def gettypecode(x):
    return x.dtype.char

def oldnonzero(a):
    return a.nonzero()[0]

# Import the floraC shared object which contains all of UEDGE

# --- Set default runid to first filename in the command line, stripping off
# --- the .py suffix.
if sys.argv[0]:
    if sys.argv[0][-3:] == '.py':
        h, t = os.path.split(sys.argv[0][:-3])
        runid = t
        del h, t
    else:
        h, t = os.path.split(sys.argv[0])
        runid = t
        del h, t

try:

   class MyPrompt(Prompts):
     def in_prompt_tokens(self, cli=None):
         return [(Token.Prompt, 'Flora>>> ')]
     def out_prompt_tokens(self, cli=None):
         return [(Token.Prompt, 'Flora>>> ')]

   get_ipython
except:
   sys.ps1='Flora>>> '
else:
   ip = get_ipython()
   ip.prompts = MyPrompt(ip)


##############################################################################
######  Don't put anything below this line!!! ################################
##############################################################################
