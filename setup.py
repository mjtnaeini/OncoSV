#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='OncoSV',
    version='1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'OncoSV=OncoSV.cli:main',
        ],
    },
    author='Marjan Naeini',
    author_email='m.naeini@garvan.org.au',
    description='A comprehensive tool for consensus structural variant calling, somatic variant calling, and complex structural variant network analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'setuptools>=50.0.0',
        'numpy>=1.22.4',
        'pandas>=2.0.0',
        'argparse>=1.4.0',
        'pysam>=0.22.0',
        'networkx>=2.4',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
