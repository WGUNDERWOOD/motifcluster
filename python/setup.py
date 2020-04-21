from setuptools import setup

with open("../README.md", 'r') as f:
    long_description = f.read()

setup(
   name='motifcluster',
   version='0.0.1',
   description='Motif-Based Spectral Clustering of Weighted Directed Networks',
   license="GPLv3",
   long_description=long_description,
   author='William George Underwood',
   author_email='wgu2@princeton.edu',
   packages=['motifcluster'],
   install_requires=[], #external packages as dependencies
)
