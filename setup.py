
from setuptools import setup, Extension, find_packages

#############################################################################
### check for Python.h
#############################################################################
from distutils.command.config import config as _config

# a faux config command that checks for the Python.h header file
class config(_config):
  def run(self):
    self.check_python_dev()
    _config.run(self)

  def check_python_dev(self):
    from distutils import sysconfig
    ok = self.check_header('Python.h',include_dirs=[sysconfig.get_python_inc()])
    if not ok:
      from distutils.errors import DistutilsPlatformError
      errmsg = ("The compiler cannot find the 'Python.h' header file.\n"
                "Please check your configuration to see if you have python-dev installed.")
      raise DistutilsPlatformError(errmsg)

# a faux build_ext command that calls the faux config command
from distutils.command.build_ext import build_ext as _build_ext
class build_ext(_build_ext):
  def run(self):
    self.run_command('config')
    _build_ext.run(self)

#############################################################################
### define the C extension for the original mfinder code
#############################################################################
mfinder = Extension('_mfinder',
                    sources=['/'.join(['pymfinder','mfinder',f]) \
                    for f in ['clustering.c',
                              'globals.c',
                              'grassberger.c',
                              'hash.c',
                              'list.c',
                              'mat.c',
                              'metropolis.c',
                              'motif_ids.c',
                              'output.c',
                              'permutation.c',
                              'prob.c',
                              'random.c',
                              'results.c',
                              'role.c',
                              'stubs.c',
                              'switches.c',
                              'wrapper.c',
                              #'mfinder.i',     # not required unless the user modifies the swig interface
                              'mfinder_wrap.c', # preferred since the user does not need to use swig directly
                              ]
                            ],
                    define_macros=[('UNIX', None),],
                    extra_compile_args = ["-O3",],
                    )

#############################################################################
#### the pymfinder setup
#############################################################################
setup(
    name = "pymfinder",
    version = "1.0",
    description = "Python wrapper for mfinder 1.2",
    author = "Daniel B. Stouffer",
    author_email = "daniel.stouffer@canterbury.ac.nz",
    url = 'http://github.com/stoufferlab/pymfinder',
    packages = find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    ext_package = 'pymfinder.mfinder',
    ext_modules = [mfinder,],
    cmdclass={'build_ext': build_ext,
              'config': config,
              },
    test_suite = 'tests.test_pymfinder',
)
