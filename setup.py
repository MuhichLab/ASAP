from setuptools import setup, find_packages

setup(
    name='ASAPy',
    version='1.1.0',
    author='Steven Wilson',
    author_email='sawilso6@asu.edu',
    url='https://github.com/MuhichLab/ASAP',
    packages=["ASAP",],
   )


'''
add this later
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[ 
        "sympy",
        "pymatgen>=2022.0.17"
    ],
'''
