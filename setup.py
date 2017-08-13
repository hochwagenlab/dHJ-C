#! /usr/bin/env python

import os
import subprocess

from setuptools import setup, find_packages

VERSION = '0.1'

if os.environ.get('DEVDIST'):
    GIT_REV_HEAD = ['git', 'rev-parse', 'HEAD']
    VERSION = VERSION + subprocess.check_output(GIT_REV_HEAD).strip()[:7]

setup(
    name="dhj_c",
    version=VERSION,
    url='https://github.com/hochwagenlab/dHJ-C',
    author='Hochwagen Lab',
    author_email='andi@nyu.edu',
    description="dHJ-C Analysis tools",
    packages=find_packages(),
    tests_require=[
        # 'pytest',
        # 'mock',
    ],
    entry_points='''
        [console_scripts]
        dhj-c=dhj_c_cli:cli
    ''',
)