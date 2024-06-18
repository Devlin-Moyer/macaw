from setuptools import find_packages, setup

setup(
    name="macaw",
    description="Metabolic Accuracy Checks And Workflow (MACAW) is a "
    "collection of tests for automatically highlighting "
    "reactions in an arbitrary Genome-Scale Metabolic Model "
    "(GSMM) that are likely to contain or be near an error of "
    "some kind",
    url="https://github.com/Devlin-Moyer/macaw",
    license="MIT",
    version="0.1.0",
    packages=find_packages(),
    install_requires = [
        'cobra==0.29.0',
        'networkx==3.1',
        'pandas==2.0.3',
        'Pebble==5.0.7'
    ]
)
