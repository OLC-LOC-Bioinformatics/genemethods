#!/usr/bin/env python
import os


def test_imports():
    """
    Test the imports
    """
    for root, dirs, files in os.walk('.'):
        for f in files:
            package = os.path.basename(root)
            module = os.path.splitext(f)[0]
            if f.endswith('.py') and '__' not in f and 'test_' not in f and root != '.' and f != 'setup.py':
                import_statement = 'import {package}.{module}'.format(package=package,
                                                                      module=module)
                if 'blaster' not in import_statement and 'tests' not in import_statement:
                    print(import_statement)
                    exec(import_statement)
                    # Try all the imports from within the file
                    with open(os.path.join(root, f), 'r') as python_file:
                        for line in python_file:
                            if line.startswith('from') or line.startswith('import'):
                                statement = line.rstrip()
                                print(root, f, statement)
                                if len(statement.split(',')) == 1:
                                    print(line.rstrip())
                                    exec(line.rstrip())
                                else:
                                    # from accessoryFunctions.accessoryFunctions import make_path, MetadataObject
                                    if line.startswith('from'):
                                        import_base, import_line = line.split(' import ')
                                        imported_packages = import_line.split(',')
                                        for imported_package in imported_packages:
                                            if imported_package != '/':
                                                print('!!!!', import_base + ' import ', imported_package)
                                    # import sys, os
                                    else:
                                        pass
