from setuptools import setup

setup(
   name='pybragg',
   version='0.1',
   description='An implementation of the [Bortfeld bragg curve fit](https://pubmed.ncbi.nlm.nih.gov/9434986/) in python with scipy and numpy.',
   author='Florian Mentzel',
   author_email='florian.mentzel@tu-dortmund.de',
   install_requires=['scipy', 'numpy'],
   packages=['pybragg']
)
