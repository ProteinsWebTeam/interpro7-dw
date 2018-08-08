#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
from setuptools import setup, find_packages


def get_version():
    filepath = os.path.join(os.path.dirname(__file__), 'i7dw', '__init__.py')

    with open(filepath) as fh:
        text = fh.read()

    m = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', text, re.M)
    return m.group(1)


def get_requirements():
    filepath = os.path.join(os.path.dirname(__file__), 'requirements.txt')

    with open(filepath) as fh:
        requirements = fh.read().splitlines()

    return requirements


setup(
    name='i7dw',
    version=get_version(),
    description='',
    packages=find_packages(),
    zip_safe=False,
    install_requires=get_requirements(),
    entry_points={
        'console_scripts': [
            'buildi7dw = i7dw.cli:cli',
        ]
    }
)
