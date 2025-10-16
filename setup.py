import os
import re
import io
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

__version__ = "2025.10.17"


# Get requirements from requirements-dev.txt file
here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, 'requirements-dev.txt'), "r", encoding="utf-8") as f:
    requirements_dev = f.read().replace('==', '>=').splitlines()

# Load readme file for Python package
with io.open(os.path.join(here, "python_solnp", 'README.md'), "r", encoding="utf-8") as file:
    long_description = file.read()


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1)
            )
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(
            self.get_ext_fullpath(ext.name)))
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DBUILD_PYSOLNP=TRUE'
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                    cfg.upper(), extdir
                ),
                '-DPYBIND3=TRUE'  # Python 3.14t only builds with latest pybind on windows, but changing it breaks Macosx builds, so this is the workaround
            ]
            if platform.machine().lower() == "arm64":
                # For ARM64 builds, set the architecture accordingly
                cmake_args += ['-A', 'ARM64']
            elif sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            else:
                cmake_args += ['-A', 'Win32']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        if platform.system() == "Darwin":
            arch = os.environ.get("CIBW_ARCHS_MACOS")
            if arch:
                cmake_args.append(f"-DCMAKE_OSX_ARCHITECTURES={arch}")

        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ['cmake', '--build', '.'] + build_args,
            cwd=self.build_temp
        )


setup(name='pysolnp',
      version=__version__,
      author='Krister Sune Jakobsson',
      author_email='krister.s.jakobsson@gmail.com',
      description='This provides the SOLNP optimization Algorithm.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/KristerSJakobsson/solnp',
      license='Boost Software License',
      ext_modules=[CMakeExtension('pysolnp')],
      cmdclass=dict(build_ext=CMakeBuild),
      zip_safe=False,
      packages=find_packages(),
      include_package_data=True,
      classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Programming Language :: Python :: 3.11",
          "Programming Language :: Python :: 3.12",
          "Programming Language :: Python :: 3.13",
          "Programming Language :: Python :: 3.14",
          "Programming Language :: C++",
          "Operating System :: MacOS :: MacOS X",
          "Operating System :: Microsoft :: Windows",
          "Operating System :: POSIX :: Linux",
          "Development Status :: 6 - Mature"
      ],
      extras_require={
          'dev': requirements_dev,
      },
      data_files=[(
          '.', [
              'requirements-dev.txt',
          ]
      )],
      python_requires='>=3.8'
      )
