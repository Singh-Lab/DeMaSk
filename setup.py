from setuptools import setup

setup(
    name='demask',
    author='Daniel Munro',
    version='1.0',
    packages=['demask'],
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
