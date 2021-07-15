#!/usr/bin/env python

from setuptools import setup

setup(
    name='band-dot',
    version='1.0',
    description='Simple band matching for condensed matter and quantum physics',
    author='Michael Lamparski',
    author_email='diagonaldevice@gmail.com',
    url='https://github.com/ExpHP/band-dot/',
    packages=['band_dot'],
    requires=['numpy'],
    scripts=['bin/band-dot'],
)
