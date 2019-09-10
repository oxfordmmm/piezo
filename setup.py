from setuptools import setup
from piezo import __version__

setup(
    name='piezo',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    packages=['piezo'],
    package_data={'': ['../config/*']},
    install_requires=[
        "numpy >= 1.13",
        "pandas >= 0.23.1",
    ],
    license='University of Oxford',
    scripts=['bin/piezo-vcf-parse.py'],
    long_description=open('README.md').read(),
    zip_safe=False
)
