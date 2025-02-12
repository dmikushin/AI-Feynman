"""Setup for the AI Feynman software components
 Copyright (C) 2015  Ehsan Azar (dashesy@linux.com)
 Copyright (C) 2019, 2020  Dmitry Mikushin (dmitry@kernelgen.org)

This file basically just uses CMake to compile the given project.

To build the package:
    python3 setup.py build
To build and install:
    python3 setup.py install
To package the wheel (after pip installing twine and wheel):
    python3 setup.py bdist_wheel
To upload the binary wheel to PyPi
    twine upload dist/*.whl
To upload the source distribution to PyPi
    python3 setup.py sdist 
    twine upload <project-name>/<project-name>-*.tar.gz
To exclude certain options in the cmake config use --no:
    for example:
    --no USE_AVX_INSTRUCTIONS: will set -DUSE_AVX_INSTRUCTIONS=no
Additional options:
    --compiler-flags: pass flags onto the compiler, e.g. --compiler-flags "-Os -Wall" passes -Os -Wall onto GCC.
    -G: Set the CMake generator.  E.g. -G "Visual Studio 14 2015"
    --clean: delete any previous build folders and rebuild.  You should do this if you change any build options
             by setting --compiler-flags or --no since the last time you ran a build.  This will
             ensure the changes take effect.
    --set: set arbitrary cmake options e.g. --set CUDA_HOST_COMPILER=/usr/bin/gcc-6.4.0
           passes -DCUDA_HOST_COMPILER=/usr/bin/gcc-6.4.0 to CMake.
"""

project_name = 'feynman'

import os
import re
import sys
import shutil
import platform
import subprocess
import multiprocessing
from distutils import log, dir_util
from math import ceil,floor

from setuptools import setup, Distribution, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install_lib import install_lib
from distutils.version import LooseVersion

from setuptools.command.test import test as TestCommand

def get_extra_cmake_options():
    """read --clean, --no, --set, --compiler-flags, and -G options from the command line and add them as cmake switches.
    """
    _cmake_extra_options = []
    _clean_build_temp = False

    opt_key = None

    has_generator = False

    argv = [arg for arg in sys.argv]  # take a copy
    # parse command line options and consume those we care about
    for arg in argv:
        if opt_key == 'compiler-flags':
            _cmake_extra_options.append('-DCMAKE_CXX_FLAGS={arg}'.format(arg=arg.strip()))
        elif opt_key == 'G':
            has_generator = True
            _cmake_extra_options += ['-G', arg.strip()]
        elif opt_key == 'no':
            _cmake_extra_options.append('-D{arg}=no'.format(arg=arg.strip()))
        elif opt_key == 'set':
            _cmake_extra_options.append('-D{arg}'.format(arg=arg.strip()))

        if opt_key:
            sys.argv.remove(arg)
            opt_key = None
            continue

        if arg == '--clean':
            _clean_build_temp = True
            sys.argv.remove(arg)
            continue

        if arg == '--yes':
            print("The --yes options to setup.py don't do anything since all these options ")
            print("are on by default.  So --yes has been removed.  Do not give it to setup.py.")
            sys.exit(1)
        if arg in ['--no', '--set', '--compiler-flags']:
            opt_key = arg[2:].lower()
            sys.argv.remove(arg)
            continue
        if arg in ['-G']:
            opt_key = arg[1:]
            sys.argv.remove(arg)
            continue

    # If no explicit CMake Generator specification,
    # prefer Ninja on Windows
    if (not has_generator) and (platform.system() == "Windows") and shutil.which("ninja"):
        _cmake_extra_options += ['-G', "Ninja"]

    # Enforce CMake to use the same version of Python, as used for this script
    _cmake_extra_options.append('-DPython3_VERSION={major}.{minor}'.format(major=sys.version_info.major,minor=sys.version_info.minor))

    return _cmake_extra_options, _clean_build_temp

cmake_extra_options,clean_build_temp = get_extra_cmake_options()


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

def rmtree(name):
    """remove a directory and its subdirectories.
    """
    def remove_read_only(func, path, exc):
        excvalue = exc[1]
        if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
            os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            func(path)
        else:
            raise

    if os.path.exists(name):
        log.info('Removing old directory {}'.format(name))
        shutil.rmtree(name, ignore_errors=False, onerror=remove_read_only)


class CMakeBuild(build_ext):

    def get_cmake_version(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("\n*******************************************************************\n" +
                                  " CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions) + 
                               "\n*******************************************************************\n")
        return re.search(r'version\s*([\d.]+)', out.decode()).group(1)

    def run(self):
        cmake_version = self.get_cmake_version()
        if platform.system() == "Windows":
            if LooseVersion(cmake_version) < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        cmake_args = []

        cmake_args += cmake_extra_options 

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() != "Windows":
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            # Do a parallel build
            build_args += ['--', '-j'+str(num_available_cpu_cores(2))]

        build_temp = os.path.abspath(self.build_temp)

        if clean_build_temp:
            rmtree(build_temp)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        cmake_setup = ['cmake', ext.sourcedir] + cmake_args
        cmake_build = ['cmake', '--build', '.'] + build_args

        print("Building extension for Python {}".format(sys.version.split('\n',1)[0]))
        print("Invoking CMake setup: '{}'".format(' '.join(cmake_setup)))
        sys.stdout.flush()
        subprocess.check_call(cmake_setup, cwd=build_temp)
        print("Invoking CMake build: '{}'".format(' '.join(cmake_build)))
        sys.stdout.flush()
        subprocess.check_call(cmake_build, cwd=build_temp)

        # Do in-source install as well, in order to have it visible when
        # package is installed in "develop(ment) mode" (-e --editable)
        install_dir = os.path.join(os.path.abspath(ext.sourcedir), project_name)

        cmake_install_prefix = ['cmake', '-DCMAKE_INSTALL_PREFIX=' + install_dir, '-P', 'cmake_install.cmake' ]

        # Adjust install prefix as shown at LLVM and not widely known:
        # https://llvm.org/docs/CMake.html#id6
        print("Adjusting CMake install prefix: '{}'".format(' '.join(cmake_install_prefix)))
        sys.stdout.flush()
        subprocess.check_call(cmake_install_prefix, cwd=build_temp)

class CMakeInstall(install_lib):

    def install(self):

        build_cmd = self.get_finalized_command('build_ext')
        build_files = build_cmd.get_outputs()
        build_temp = getattr(build_cmd, 'build_temp')

        install_dir = os.path.join(os.path.abspath(self.install_dir), project_name)

        cmake_install_prefix = ['cmake', '-DCMAKE_INSTALL_PREFIX=' + install_dir, '-P', 'cmake_install.cmake' ]

        # Adjust install prefix as shown at LLVM and not widely known:
        # https://llvm.org/docs/CMake.html#id6
        print("Adjusting CMake install prefix: '{}'".format(' '.join(cmake_install_prefix)))
        sys.stdout.flush()
        subprocess.check_call(cmake_install_prefix, cwd=build_temp)

        if (os.path.exists(project_name)):
            dir_util.copy_tree(project_name, install_dir)

def num_available_cpu_cores(ram_per_build_process_in_gb):
    if 'TRAVIS' in os.environ and os.environ['TRAVIS']=='true':
        # When building on travis-ci, just use 2 cores since travis-ci limits
        # you to that regardless of what the hardware might suggest.
        return 2 
    try:
        mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  
        mem_gib = mem_bytes/(1024.**3)
        num_cores = multiprocessing.cpu_count() 
        # make sure we have enough ram for each build process.
        mem_cores = int(floor(mem_gib/float(ram_per_build_process_in_gb)+0.5));
        # We are limited either by RAM or CPU cores.  So pick the limiting amount
        # and return that.
        return max(min(num_cores, mem_cores), 1)
    except ValueError:
        return 2 # just assume 2 if we can't get the os to tell us the right answer.

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = '--ignore docs --ignore ' + os.path.join('ThirdParty', project_name)

    def run_tests(self):
        import shlex
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)

def read_entire_file(fname):
    """Read text out of a file relative to setup.py.
    """
    return open(os.path.join(fname)).read()

setup(
    name=project_name,
    version='2.0.0',
    description='AI Feynman: a Physics-Inspired Method for Symbolic Regression',
    url='https://github.com/SJ001/AI-Feynman',
    author='Silviu-Marian Udrescu',
    author_email='sudrescu@mit.edu',
    license='MIT',
    packages=find_packages(),
    package_data={project_name: ['*.txt']},
    include_package_data=True,
    ext_modules=[CMakeExtension('feynman._feynman')],
    cmdclass=dict(build_ext=CMakeBuild, install_lib=CMakeInstall, test=PyTest),
    python_requires='>3.6',
    install_requires=['matplotlib',
                      'numpy',
                      'seaborn',
                      'sklearn',
                      'sortedcontainers',
                      'sympy >= 1.4',
                      'torch >= 1.4.0',
                      'torchvision',
                      ],
    entry_points={
        'console_scripts': [
            'feynman_sr1 = feynman.symbolic_regress1:main',
            'feynman_sr2 = feynman.symbolic_regress2:main',
            'feynman_sr3 = feynman.symbolic_regress3:main',
        ]},
    distclass=Distribution,
)

