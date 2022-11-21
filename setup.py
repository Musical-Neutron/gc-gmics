#!/usr/bin/env python3

from setuptools import setup
from sphinx.setup_command import BuildDoc
# from distutils.extension import Extension

cmdclass = {'build_sphinx': BuildDoc}
ext_modules = []

name = 'gc_paper_plots'
version = '0.1'
release = '0.1.0'

setup(
    name=name,
    version=version,
    description='DESCRIPTION',
    author='Oliver Newton',
    author_email='onewton@cft.edu.pl',
    packages=['gc_paper_plots'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    command_options={
        'build_sphinx': {
            'project': ('setup.py', name.replace('_', ' ').title()),
            'version': ('setup.py', version),
            'release': ('setup.py', release),
            'source_dir': ('setup.py', 'doc')
        }
    },
)
