from setuptools import setup, find_packages

import nested

setup(
    name='nested',
    version=nested.__version__,
    description='Nested description',
    author=nested.__author__,
    packages=find_packages(),
    install_requires=[
        'bcbio-gff==0.6.4',
        'biopython==1.70',
        'click==6.7',
        'networkx==2.1',
        'PyYAML==3.12'
    ],
    entry_points={
        'console_scripts': [
            'nested-generator=nested.cli.generator:main',
            'nested-nester=nested.cli.nester:main'
        ]
    }
)
