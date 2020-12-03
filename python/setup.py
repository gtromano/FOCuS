from distutils.core import setup, Extension
from Cython.Build import cythonize

# you specify the c source file in the sources list
ext = Extension('FOCuS', sources = ['python_bindings.pyx', 'python_bindings.cpp'], language="c++")
setup(name="FOCuS", ext_modules = cythonize([ext], language_level = "3"))