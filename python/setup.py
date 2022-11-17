import setuptools

with open("README.rst", "r") as f:
  long_description = f.read()

setuptools.setup(
  name='motifcluster',
  version='0.2.3',
  description='Motif-Based Spectral Clustering of Weighted Directed Networks',
  license="GPLv3",
  long_description=long_description,
  long_description_content_type="text/x-rst",
  author='William George Underwood, Andrew Elliott',
  author_email='wgu2@princeton.edu',
  url="https://github.com/WGUNDERWOOD/motifcluster",
  packages=['motifcluster'],
  install_requires=[
    "networkx>=2.4",
    "numpy>=1.18.3",
    "scipy>=1.4.1",
    "scikit-learn>=0.22.2"
  ],
)
