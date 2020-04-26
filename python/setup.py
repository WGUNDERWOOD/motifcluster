from setuptools import setup

setup(
<<<<<<< Updated upstream
    name='motifcluster',
    version='0.0.1',
    description='Motif-Based Spectral Clustering of Weighted Directed Networks',
    license="GPLv3",
    long_description=("Construct motif adjacency matrices for (weighted directed)",
        "networks, and use them for spectral clustering. Also provides random sampling",
        "methods for weighted directed networks."),
=======
<<<<<<< Updated upstream
   name='motifcluster',
   version='0.0.1',
   description='Motif-Based Spectral Clustering of Weighted Directed Networks',
   license="GPLv3",
   long_description="TODO",
   author='William George Underwood',
   author_email='wgu2@princeton.edu',
   packages=['motifcluster'],
   install_requires=[], #external packages as dependencies
=======
    name='motifcluster',
    version='0.0.1',
    description='Motif-Based Spectral Clustering of Weighted Directed Networks',
    license="GPLv3", #TODO licence file?
    long_description=("Construct motif adjacency matrices for (weighted directed)",
        "networks, and use them for spectral clustering. Also provides random sampling",
        "methods for weighted directed networks."), #TODO move this to a md file in python/
>>>>>>> Stashed changes
    author='William George Underwood',
    author_email='wgu2@princeton.edu',
    packages=['motifcluster'],
    install_requires=[
        "networkx",
        "numpy",
        "scipy",
        "sklearn"
    ],
<<<<<<< Updated upstream
=======
>>>>>>> Stashed changes
>>>>>>> Stashed changes
)
