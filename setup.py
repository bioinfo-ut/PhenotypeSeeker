#!/usr/bin/env python3

from setuptools import setup

setup(name = 'PhenotypeSeeker',
      version ='0.5.0',
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
          'scipy', 'Biopython', 'scikit-learn', 'xgboost', 'multiprocess',
          'pandas'
          ],
      )


