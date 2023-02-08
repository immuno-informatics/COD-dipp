#!/usr/bin/env python

'''
setup.py file for Scavager
'''
from setuptools import setup
version = open('VERSION').readline().strip()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name             = 'Scavager',
    version          = version,
    description      = '''Proteomics post-search algorithm''',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    author           = 'Mark Ivanov & Lev Levitsky & Julia Bubis',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'https://bitbucket.org/markmipt/scavager',
    packages         = ['scavager', ],
    install_requires = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Programming Language :: Python :: 3',
                        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    license          = 'License :: OSI Approved :: Apache Software License',
    entry_points     = {'console_scripts': ['scavager = scavager.search:run',
                                            'scav2diffacto = scavager.scav2diffacto:run',
                                            'scav2nsaf = scavager.scav2nsaf:run']}
    )
