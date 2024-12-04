from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from adecoi import __version__, _program


setup(name='adecoi',
      version=__version__,
      packages=find_packages(),
      scripts=[
            "adecoi/scripts/Snakefile",
            "adecoi/scripts/cns_runner.smk"
            ],
      package_data={"adecoi":["data/*"]},
      install_requires=[
            "biopython>=1.70",
            "mako>=1.1",
            "matplotlib",
            "seaborn"
        ],
      description='ADE COI',
      url='https://github.com/aineniamh/adecoi',
      author='Aine OToole, Rachel Colquhoun & Rambaut Group',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = adecoi.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
