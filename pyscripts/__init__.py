from .flora import *
from os import path
from pathlib import Path
try:
    from flora.__version__ import __version__
    from flora.__src__ import __src__
except:
    try:
        from __version__ import __version__
        from __src__ import __src__
    except:
        __version__ = 'unknown'
        __src__ = 'unknown'

#
# Load the startup file .florarc.py from cwd or home.
#
_homepath = path.expanduser('~')
_homefile = Path('{}/.florarc.py'.format(_homepath))
_localpath = path.expanduser('.')
_localfile = Path('{}/.florarc.py'.format(_localpath))

if path.exists(_localfile):
   with open(_localfile) as f:
      exec(open(_localfile).read())
elif path.exists(_homefile):
   with open(_homefile) as f:
      exec(open(_homefile).read())
      
