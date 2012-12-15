
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
        'phenoseq_kaks = phenoseq.ka_ks:run_all',
        'phenoseq_hypergeom = phenoseq.hypergeometric:main',
        ],
    }


def check_deps(ignore_pygr=False):
    'check whether our dependencies are installed'
    deps = []
    try:
        import numpy
    except ImportError:
        deps.append('numpy')
    try:
        from scipy import stats
    except ImportError:
        deps.append('scipy')
    try:
        from Bio import SeqIO
    except ImportError:
        deps.append('biopython')
    if not ignore_pygr:
        try:
            from pygr import cnestedlist
        except ImportError:
            deps.append('pygr')
    return deps

def try_install(**kwargs):
    'try to install phenoseq using setup()'
    setup(
        name = 'phenoseq',
        version= '1.0',
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
