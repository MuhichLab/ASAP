from setuptools import setup, find_packages

setup(
    name='ASAP',
    version='1.0.0',
    author='Steven Wilson',
    author_email='sawilso6@asu.edu',
    url='https://github.com/MuhichLab/ASAP',
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[ 
        "sympy",
        "pymatgen>=2022.0.17"
    ],
   )
