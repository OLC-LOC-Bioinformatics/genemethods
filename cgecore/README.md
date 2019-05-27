# cge_core_module

Core module for the Center for Genomic Epidemiology

This module contains classes and functions needed to run the service wrappers and pipeline scripts

The pypi project can be found here:
https://pypi.org/project/cgecore/

# How to update:
1. Make changes to the modules
2. Bump the version number accordingly in cgecore/__init__.py
3. Install package locally
4. Test the changes locally (for both python2 and python3)
5. Distribute to Pypi

# Install package locally
python2 setup.py install
python3 setup.py install

# Distribute to PyPi
python setup.py sdist upload -r pypi
