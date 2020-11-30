from setuptools import setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='demask',
    author='Daniel Munro',
    author_email='dan@dmun.ro',
    version='1.0',
    description='Amino acid substitution impact prediction.',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/Singh-Lab/DeMaSk',
    packages=['demask'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'configargparse',
        'numpy',
    ],
    extras_require={
        'tests': [
        ],
        'docs': [
            'sphinx',
            'sphinx_rtd_theme',
            'sphinx_autodoc_typehints',
        ]
    }
)
