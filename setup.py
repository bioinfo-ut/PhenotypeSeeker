#!/usr/bin/env python3

from setuptools import setup

setup(name = 'PhenotypeSeeker',
      version ='0.7.3',
      description = 'Identifies phenotype-specific k-mers, creates phenotype \
          prediction models and using this model predicts phenotype from \
          sequenced bacterial strains',
      author = 'Erki Aun',
      author_email = 'erki.aun@ut.ee',
      url = 'https://www.bioinfo.ut.ee/PhenotypeSeeker',
      py_modules = [
          'PhenotypeSeeker/modeling', 'PhenotypeSeeker/prediction',
          'PhenotypeSeeker/annotation', 'PhenotypeSeeker/__init__'
          ],
      scripts = ['scripts/phenotypeseeker'],
      install_requires = [
          'numpy==1.18.1', 'scipy==1.4.1', 'Biopython==1.76', 'scikit-learn==0.22.1', 'xgboost==1.0.1', 'multiprocess==0.70.9',
          'pandas==1.0.1', 'ete3==3.1.1', 'matplotlib==2.1.0'
          ],
      )


