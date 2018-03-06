#!/usr/bin/env python2.7

from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup(name = 'PhenotypeSeeker',
      version ='0.0.1',
      description = 'Identifies phenotype-specific k-mers, creates phenotype \
          prediction models and using this model predicts phenotype from \
          sequenced bacterial strains',
      author = 'Erki Aun',
      author_email = 'erki.aun@ut.ee',
      url = 'https://www.bioinfo.ut.ee/PhenotypeSeeker',
      py_modules = [
          'PhenotypeSeeker/modeling', 'PhenotypeSeeker/prediction',
          'PhenotypeSeeker/__init__'
          ],
      scripts = ['scripts/phenotypeseeker'],
      cmdclass={'build_ext':build_ext},
      setup_requires=['numpy'],
      install_requires = [
          'scipy', 'Biopython', 'cogent', 'scikit-learn'
          ],
      )


