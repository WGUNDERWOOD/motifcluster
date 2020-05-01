import setuptools

with open("README_python.rst", "r") as f:
  long_description = f.read()

setuptools.setup(
  name='motifcluster',
  version='0.0.2',
  description='Motif-Based Spectral Clustering of Weighted Directed Networks',
  license="GPLv3",
  long_description=long_description,
  long_description_content_type="text/x-rst",
  author='William George Underwood, Andrew Elliott',
  author_email='wgu2@princeton.edu',
  packages=['motifcluster'],
  install_requires=[
    "networkx",
    "numpy",
    "scipy",
    "sklearn"
  ],
)
