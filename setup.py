from setuptools import setup, find_packages

setup(
    name='paws',
    version='1.0',
    packages=find_packages(),
    install_requires=[
    "numpy", "astropy", "pathlib", "tqdm", "scipy", "pyfstat", "matplotlib" 
    ],
    author='Damon Cheung',
    author_email='damoncht@umich.edu',
    description='A package to manage: to create condor jobs and analysis the data for directed search of continuous gravitational wave signals using WEAVE',
    url='https://github.com/damondmc/paws',  
)

