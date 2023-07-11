# flora

This requires gfortran. It's been tested with Gfortran 11.2. It does use a compiler directive that needs a gfortran compiler >4.8.5. The floor has not been determined.

To build:

pip install forthon
python setup.py build install

To run the code:

import flora
flora.glr.glrgen()

To setup case set variables with:
from flora import glr
glr.<varname1> = ....
glr.<varname2> = ....

For a list of variables use 
glr.varlist()

Because of the way the external objects are created a dir(glr) will not reveal the variable names.
