#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys
import os
import os.path
import string
import site
from Forthon.compilers import FCompiler
import getopt
import logging

version='1.0.0'

try:
    os.environ['PATH'] += os.pathsep + site.USER_BASE + '/bin'
    import setuptools
    import distutils
    from distutils.core import setup
    from distutils.core import Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
    from distutils.command.install import install
    from subprocess import call
    import numpy
except:
    raise SystemExit("Distutils problem")

optlist, args = getopt.getopt(sys.argv[1:], 'g', ['debug'])
machine = sys.platform
debug = 0
fcomp = None
parallel = 0
petsc = 0

for o in optlist:
    if o[0] == '-g':
        debug = 1
    if o[0] == '--debugg':
        debug = 1
        

fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)


class floraInstall(build):
    def run(self):
       install.run(self)
       logging.basicConfig(stream=sys.stderr,level=logging.INFO)
       log = logging.getLogger()
       log.info("test")
class floraBuild(build):
    def run(self):
        # with python2 everything is put into a single floraC.so file
        if sys.hexversion < 0x03000000:
            raise SystemExit("Python versions < 3 not supported")
        else:
            call(['make','-i', '-f','Makefile.Forthon'])
            build.run(self)


class floraClean(build):
    def run(self):
        if sys.hexversion < 0x03000000:
            raise SystemExit("Python versions < 3 not supported")
        else:
            call(['make', '-f', 'Makefile.Forthon', 'clean'])

florapkgs = ['glr', 'utl']


def makeobjects(pkg):
    return [pkg+'_p.o', pkg+'pymodule.o']


floraobjects = []

# add here any extra dot o files other than pkg.o, pkg_p.o


if sys.hexversion < 0x03000000:
    raise SystemExit("Python versions < 3 not supported")
else:
    dummydist = Distribution()
    dummydist.parse_command_line()
    dummybuild = dummydist.get_command_obj('build')
    dummybuild.finalize_options()
    builddir = dummybuild.build_temp

floraobjects = map(lambda p: os.path.join(builddir, p), floraobjects)

library_dirs = fcompiler.libdirs
libraries = fcompiler.libs

with open('pyscripts/__version__.py','w') as ff:
    ff.write("__version__ = '%s'\n"%version)
with open('pyscripts/__src__.py','w') as ff:
    ff.write("__src__ = '%s'\n"%os.getcwd())

define_macros=[("WITH_NUMERIC", "0"),
               ("FORTHON_PKGNAME", '\"floraC\"'),
               ("FORTHON","1")]

# check for readline
rlncom = "echo \"int main(){}\" | gcc -x c -lreadline - "
rln = os.system(rlncom)
if rln == 0: 
   define_macros = define_macros + [("HAS_READLINE","1")]
   os.environ["READLINE"] = "-l readline"
   libraries = ['readline'] + libraries


setup(name="flora",
      version=version,
      author='Bruce Cohen',
      author_email="cohen1@llnl.gov",
      maintainer='Bill Meyer',
      maintainer_email='meyer8@llnl.gov',
      packages=['flora'],
      package_dir={'flora': 'pyscripts'},
      # include_package_data=True,
      scripts=['pyscripts/flora_bas2py'],
      ext_modules=[Extension('flora.floraC',
                             ['floraC_Forthon.c',
                              os.path.join(builddir, 'Forthon.c')],
                             include_dirs=[builddir, numpy.get_include()],
                             library_dirs=library_dirs,
                             libraries=libraries,
                             define_macros=define_macros,
                             extra_objects=floraobjects,
                             extra_link_args=['-g','-DFORTHON'] +
                             fcompiler.extra_link_args,
                             extra_compile_args=fcompiler.extra_compile_args
                             )],

      cmdclass={'build': floraBuild, 'clean': floraClean},
      test_suite="pytests",
      install_requires=['forthon'],
      # note that include_dirs may have to be expanded in the line above
      classifiers=['Programming Language :: Python',
                   'Programming Language :: Python :: 3']

      )
