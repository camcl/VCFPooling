from distutils.core import setup
from Cython.Build import cythonize

setup(
     ext_modules=cythonize("cydecoder.pyx", language_level="3"),
)