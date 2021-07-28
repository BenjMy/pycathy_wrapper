import sys
from setuptools import setup
# To use:
#	   python setup.py bdist --format=wininst

# trap someone trying to install pyres with something other
#  than python 2 or 3
if not sys.version_info[0] in [2, 3]:
    print('Sorry, pyCATHY not supported in your Python version')
    print('  Supported versions: 2 and 3')
    print('  Your version of Python: {}'.format(sys.version_info[0]))
    sys.exit(1)  # return non-zero value for failure

long_description = ''

try:
    import pypandoc

    long_description = pypandoc.convert('README.md', 'rst')
except:
    pass

setup(name='pyCATHY',
      description='pyCATHY is a Python package to create, run, and post-process CATHY forward and inverse models.',
      long_description=long_description,
      author='B. Mary',
      author_email='benjamin.mary@unipd.it',
      url='https://github.com/?',
      license='New BSD',
      platforms='Windows, Mac OS-X, Unix',
      install_requires=['numpy>=1.7'],
      packages=['pyCATHY'],
      version='1.0.1',
      keywords='hydrogeophysics')
