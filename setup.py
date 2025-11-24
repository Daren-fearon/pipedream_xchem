"""
Setup script for Pipedream XChem Integration
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name='pipedream-xchem',
    version='1.0.2',
    author='Daren Fearon',
    author_email='xchem@diamond.ac.uk',
    description='Automated Pipedream refinement pipeline for XChem crystallographic datasets',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/xchem/pipedream_xchem',
    license='GPL-3.0',
    py_modules=[
        'pipedream_xchem',
        'collate_pipedream_results',
        'export_pipedream',
        'pipedream_post_process',
        'pipedream_plots',
        'pipedream_thresholds',
    ],
    python_requires='>=3.10',
    install_requires=[
        'gemmi',
        'jinja2',
        'matplotlib',
        'numpy',
        'pandas',
        'paramiko',
        'plotly',
        'pyyaml',
        'rdkit',
        'scikit-learn',
        'tqdm',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Operating System :: POSIX :: Linux',
    ],
    entry_points={
        'console_scripts': [
            'pipedream-xchem=pipedream_xchem:main',
            'pipedream-collate=collate_pipedream_results:main',
            'pipedream-export=export_pipedream:main',
            'pipedream-postprocess=pipedream_post_process:main',
        ],
    },
    include_package_data=True,
    package_data={
        '': ['templates/*.html', 'Data/*.csv', '*.yaml'],
    },
    keywords='crystallography xray-diffraction protein-structure refinement xray pipedream buster',
)
