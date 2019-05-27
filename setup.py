#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(name='ase_fireworks',
      version='0.1',
      description='A fireworks interface abitrary ASE calculators',
      author='Ben Comer',
      author_email='ben.comer@gatech.edu',
      url='https://github.com/medford-group/ase_fireworks',
      packages=find_packages(),
      install_requires=['spglib', 'numpy','ase','scipy','pymongo','fireworks'],
     )
