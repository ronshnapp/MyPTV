# -*- coding: utf-8 -*-
"""
Created on Sun 20 March 2022

"""

from setuptools import find_packages, setup

setup(
    name='myptv',
    packages=find_packages(include=['myptv']),
    version='0.6.0',
    description='A 3D Particle Tracking Velocimetry library',
    install_requires=['numpy', 'scipy', 'scikit-image','pandas','matplotlib','pyyaml', 'tk', 'Pillow'],
    author='Ron Shnapp',
    author_email='ronshnapp@gmail.com',
    license='MIT',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests'
)


