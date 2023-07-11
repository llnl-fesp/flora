# Flora Equilibrium and Stabiliry Code

This requires gfortran. It's been tested with Gfortran 11.2. It does use a compiler directive that needs a gfortran compiler >4.8.5. The floor has not been determined.

To build:

pip install forthon<br>
python setup.py build install

To run the code:

import flora<br>
flora.glr.glrgen()

To setup case set variables with:<br>
from flora import glr<br>
glr.varname1 = ....<br>
glr.varname2 = ....<br>

For a list of variables use <br>
glr.varlist()

Because of the way the external objects are created a dir(glr) will not reveal the variable names.

## Release 

Flora is released under an LGPL license.  For more details see the
NOTICE and LICENSE files.

``LLNL-CODE-431811``
------
--------
