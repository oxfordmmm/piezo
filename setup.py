from setuptools import setup
from piezo import __version__

setup(
    name='piezo',
    version=__version__,
    author='Philip W Fowler',
    author_email="philip.fowler@ndm.ox.ac.uk",
    description="Predict the effect of a genetic mutation on the effect of an antibiotic using a supplied AMR catalogue.",
    url="https://github.com/philipwfowler/piezo",
    packages=['piezo'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"  ],
    python_requires='>=3.5',
    install_requires=[
        "pandas >= 0.23.1"
    ],
    license='University of Oxford, see LICENSE.md',
    scripts=['bin/piezo-predict.py'],
    zip_safe=False
)
