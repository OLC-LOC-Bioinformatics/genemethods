#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.util import convert_path
import os
__author__ = 'adamkoziol'

# Find the version
version = dict()
with open(convert_path(os.path.join('genemethods', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)


setup(
    name="genemethods",
    version=version['__version__'],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'SequenceExtractor = genemethods.SequenceExtractor.src.sequenceExtractor:main',
        ],
    },
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/genemethods",
)
