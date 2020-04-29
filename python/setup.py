import setuptools

with open("README_python.rst", "r") as fh:
  long_description = fh.read()

setuptools.setup(
    name='motifcluster',
    version='0.0.1',
    description='Motif-Based Spectral Clustering of Weighted Directed Networks',
    license="GPLv3",
    long_description=long_description,
    author='William George Underwood',
    author_email='wgu2@princeton.edu',
    packages=['motifcluster'],
    install_requires=[
      "networkx",
      "numpy",
      "scipy",
      "sklearn"
    ],
)
