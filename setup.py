from setuptools import find_packages, setup

setup(
    name='fishsim',
    version='1.0.0',
    author='Jenkin Tsui, William Choi, Luna Liu',
    author_email='jenkin.tsui@aya.yale.edu',
    description='A simulator for barcoding images',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'scipy',
        'scikit-image',
        'click',
        'pyyaml',
        'black',
        'pytest',
        'pathlib',
    ],
    entry_points={
        'console_scripts': [
            'fishsim = fishsim.cli:main',
        ],
    },
)
