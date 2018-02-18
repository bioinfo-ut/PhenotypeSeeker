#!/usr/bin/env python2.7

from setuptools import setup

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
      install_requires = [
          'Biopython', 'cogent', 'numpy', 'scipy', 'scikit-learn'
          ],
      )


