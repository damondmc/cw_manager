from setuptools import setup, find_packages

setup(
    name='cw_manager',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        # List your package dependencies here
    ],
    author='Damon Cheung',
    author_email='damoncht@umich.edu',
    description='A package to manage: to create condor jobs and analysis the data for directed search of continuous gravitational wave signalsr',
    url='https://github.com/damondmc/cw_manager',  # Optional
)

