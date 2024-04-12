from setuptools import setup, find_packages

with open("README.md", 'r') as fr:
	description = fr.read()

setup(
    name='FGFinder',
    version='1.0.0',
    url='https://github.com/Borest5543/Fg_Finder',
    license='MIT License',
    author='Túlio Augusto and Jefferson Richard',
    author_email='borest5543@gmail.com',
    keywords='Cheminformatics RDKit Chemistry',
    description='A package to identify functional groups in molecules.',
    long_description = description,
    long_description_content_type = "text/markdown",
    packages=['FGFinder'],
    install_requires=['pandas', 'rdkit'],
	classifiers = [
		'Intended Audience :: Developers',
		'Intended Audience :: End Users/Desktop',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: MIT License',
		'Natural Language :: English',
		'Operating System :: Unix',
		'Operating System :: Microsoft :: Windows',
		'Operating System :: MacOS',
		'Topic :: Scientific/Engineering :: Artificial Intelligence',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 3.8']
)
