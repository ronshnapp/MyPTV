# -*- coding: utf-8 -*-
"""
Created on Sun 20 March 2022

"""

from setuptools import find_packages, setup

setup(
    name='myptv',
    packages=find_packages(include=['myptv', 'myptv.fibers', 'myptv.TsaiModel', 'myptv.extendedZolof', 'myptv.makePlots', 'myptv.sheets']),
    version='1.3.4',
    description='A 3D Particle Tracking Velocimetry library',
    install_requires=['numpy', 'scipy', 'scikit-image','pandas','matplotlib','pyyaml', 'tk', 'Pillow>=9.5.0', 'moviepy==1.0.1'],
    author='Ron Shnapp',
    author_email='ronshnapp@gmail.com',
    license='MIT',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests'
)


