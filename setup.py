from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import sys
import setuptools

__version__ = '0.1a6'

base_path = os.path.dirname(__file__)


class GetIncludes(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    # If you need to link extra libraries or specify extra include directories
    # see https://docs.python.org/3/extending/building.html#building-c-and-c-extensions-with-distutils
    Extension(
        'pysolnp',
        ['python_solnp/pysolver.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            GetIncludes(),
            GetIncludes(user=True)
        ],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        compiler_type = self.compiler.compiler_type
        compiler_options = self.c_opts.get(compiler_type, [])
        link_options = self.l_opts.get(compiler_type, [])
        if compiler_type == 'unix':
            compiler_options.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            compiler_options.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                compiler_options.append('-fvisibility=hidden')
        elif compiler_type == 'msvc':
            compiler_options.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = compiler_options
            ext.extra_link_args = link_options
        build_ext.build_extensions(self)


with open("python_solnp/README.md", "r") as file:
    long_description = file.read()


setup(name='pysolnp',
      version=__version__,
      author='Krister Sune Jakobsson',
      author_email='krister.s.jakobsson@gmail.com',
      description='This provides the SOLNP optimizaiton Algorithm.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/KristerSJakobsson/cpp-solnp',
      license='Boost Software License',
      ext_modules=ext_modules,
      install_requires=['pybind11>=2.4'],
      setup_requires=['pybind11>=2.4'],
      cmdclass={'build_ext': BuildExt},
      zip_safe=False,
      include_package_data=True,
      classifiers=[
          "Programming Language :: Python",
          "Programming Language :: C++",
          "License :: OSI Approved :: Boost Software License 1.0 (BSL-1.0)",
          "Operating System :: MacOS :: MacOS X",
          "Operating System :: Microsoft :: Windows",
          "Operating System :: POSIX :: Linux"
      ]
      )
