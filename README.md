# Flora Equilibrium and Stability Code


To build:

This requires gfortran. It requires some compiler switches not available in gforran < 10. It's been tested with 10.0.1, 11.2, 11.3, and 12.2.

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

Running:

There are a couple of cases in the examples directory. Run with:

python test1.py <br>
python test2.py <br>

Compare with test1_ref.log and test2_ref.log. These logs are very old, compare for general agreement only. 

# Release 

Flora is released under an LGPL license.  For more details see the
NOTICE and LICENSE files.

``LLNL-CODE-431811``
------
--------
