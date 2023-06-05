from setuptools import setup
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import numpy

macros = [('_USE_MATH_DEFINES', None)]
compile_args = ['-Wno-cpp', '-std=c++17', '-fopenmp', '-Wno-format', '-DMS_WIN64']
link_args = ['-fopenmp']
link_args_static = [
    '-static-libgcc',
    '-static-libstdc++',
    '-Wl,-Bstatic,--whole-archive',
    '-lwinpthread',
    '-lgomp',
    '-Wl,--no-whole-archive',
]
# Disable NumPy warning
macros.append(('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'))
compiler_directives = {'binding': True, 'cdivision': True}

class Build(build_ext):

    def build_extensions(self):
        compiler = self.compiler.compiler_type
        for e in self.extensions:
            e.define_macros = macros
            e.extra_compile_args = compile_args
            if compiler == 'mingw32':
                e.extra_link_args = link_args_static 
            else:
                e.extra_link_args = link_args
        super(Build, self).build_extensions()

setup(
    ext_modules = cythonize("volume_cython_ex1.pyx"),
    cmdclass={'build_ext': Build},
    include_dirs=[numpy.get_include()]
)



