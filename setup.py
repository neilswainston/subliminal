'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from setuptools import find_packages, setup


setup(name='subliminal-py',
      version='0.1.0',
      description='subliminal-py: Python implementation of modules from the ' +
      'SuBliMinaL Toolbox',
      long_description='subliminal-py: Python implementation of modules ' +
      'from the SuBliMinaL Toolbox',
      url='https://github.com/synbiochem/subliminal-py',
      author='Neil Swainston',
      author_email='neil.swainston@manchester.ac.uk',
      license='MIT',
      classifiers=[
                   'Development Status :: 4 - Beta',
                   'Intended Audience :: Developers',
                   'Topic :: Software Development :: Build Tools',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 2.7'
      ],
      keywords='systems biology',
      packages=find_packages(),
      test_suite='subliminal.utils.test',
      install_requires=['glpk', 'cobra', 'synbiochem-py'])
