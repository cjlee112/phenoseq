
"""
Phenoseq
========

Phenoseq is an open source package of simulation and
analysis tools for phenotype sequencing experiments.
It can simulate and optimize the design of phenotype
sequencing experiments, and score results from
short-read data from such an experiment.
"""

import warnings
import sys
import subprocess

try:
    from setuptools import setup
    # setuptools can automatically install dependencies for you
    install_requires = ['pygr', 'biopython', 'numpy', 'scipy']
    min_requires = install_requires[:1]
    has_setuptools = True
except ImportError:
    warnings.warn('Setuptools not found, falling back to distutils')
    from distutils.core import setup
    has_setuptools = False

CLASSIFIERS = """
Development Status :: 5 - Production/Stable
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
"""

# split into lines and filter empty ones
CLASSIFIERS = filter(None, CLASSIFIERS.splitlines())

entry_points = {
    'console_scripts': [
        'phenoseq_analyze = phenoseq.analyze:main',
        'phenoseq_cost = phenoseq.simulate:main',
        'phenoseq_pathways = phenoseq.pathways:main',
        'phenoseq_kaks = phenoseq.ka_ks:main',
        'phenoseq_hypergeom = phenoseq.hypergeometric:main',
        ],
    }

def missing_dep(importCmd):
    '''We are forced to check import success in a separate process,
    because successfully installed dependencies STILL fail to import
    in this process (that installed them)!'''
    try:
        subprocess.check_output([sys.executable, '-c', importCmd])
    except subprocess.CalledProcessError:
        return True

def check_deps(ignore_pygr=False):
    'check whether our dependencies are installed'
    deps = []
    if missing_dep('import numpy'):
        deps.append('numpy')
    if missing_dep('from scipy import stats'):
        deps.append('scipy')
    if missing_dep('from Bio import SeqIO'):
        deps.append('biopython')
    if not ignore_pygr:
        if missing_dep('from pygr import cnestedlist'):
            deps.append('pygr')
    return deps

def pip_install(pkg):
    '''try to run pip install in a separate process.  Hmm, this used to
    work great, but now scipy, biopython auto-installs fail if we
    used this method to install numpy (even though numpy imports
    fine, seems to work fine).  But if we run pip install numpy
    *by hand*, everything works great. AAAAAAAAARRRRRRRRRRGGGGGGGGGHHHHHHH!'''
    try:
        subprocess.check_output(['pip', 'install', pkg])
        return True
    except subprocess.CalledProcessError:
        return False

def numpy_install_kluge():
    'setuptools fails to install numpy?!? But pip can install numpy...'
    try:
        import numpy
        return True # already installed
    except ImportError:
        print '''Your environment lacks numpy, required for phenoseq.'''
        try:
            subprocess.check_output(['pip', '--help'])
            print '''Fortunately, it can be installed using pip.
            Type the following command, then re-run the phenoseq installer:

            pip install numpy'''
        except subprocess.CalledProcessError:
            print """Unable to locate pip, so unable to install numpy for you.
            Please install numpy before proceeding.  You may do so using
            Python's pip installer, your system's package manager,
            or by downloading it yourself.  Once numpy is installed,
            re-run the phenoseq installer."""
        sys.exit() # just continuing with install doesn't work; restart.
        ## warnings.warn('Trying to install numpy using pip.  This may take a while...')
        ## if pip_install('numpy'): # see if pip can install it
        # not sure why, but if setup.py runs pip, the subsequent scipy
        # install fails, whereas if the user runs pip, everything is fine!


def try_install(**kwargs):
    'try to install phenoseq using setup()'
    setup(
        name = 'phenoseq',
        version= '0.5',
        description = 'Phenoseq, a simulation and experimental data analysis package for phenotype sequencing',
        long_description = __doc__,
        author = "Christopher Lee",
        author_email='leec@chem.ucla.edu',
        url = 'https://github.com/cjlee112/phenoseq',
        license = 'New BSD License',
        classifiers = CLASSIFIERS,

        packages = ['phenoseq'],
        **kwargs
     )


def main():
    if has_setuptools:
        if not numpy_install_kluge(): # didn't get numpy, so this may fail...
            warnings.warn('''First trying automatic install of
dependencies, even though this usually fails for numpy...\n\n''')
        try:
            try_install(install_requires=install_requires,
                        entry_points=entry_points)
        except:
            warnings.warn('''


Retrying with minimal dependencies. This should work.''')
            try_install(install_requires=min_requires,
                        entry_points=entry_points)
    else:
        try_install()
    
    requires = check_deps(has_setuptools)
    if requires:
        warnings.warn('''It appears you are missing %s,
which setuptools cannot automatically install
(most likely due to numpy auto-installation failing).
Please install the missing package(s) using
pip, easy_install, or some other method, or some
phenoseq functionality will not work.'''
                      % ' and '.join(requires))
    print 'Successfully completed phenoseq install.'
        
if __name__ == '__main__':
    main()
