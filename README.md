# PhenotypeSeeker
Identify phenotype-specific k-mers and predict phenotype using sequenced bacterial strains

## Introduction

PhenotypeSeeker identifies phenotype-specific k-mers, generates phenotype prediction model and predicts the phenotype from sequencing data. 

PhenotypeSeeker consist of two subprograms: 'PhenotypeSeeker modeling' and 'PhenotypeSeeker prediction', which both take either assembled contigs or raw-read data as an input. 


PhenotypeSeeker uses statistical model that can be trained automatically on isolates with known phenotype. The prediction of phenotypes takes less than a second per isolate if assembled genomes are used and less than a minute per isolate if raw sequencing data are used. Therefore, PhenotypeSeeker is well suited for predicting phenotypes from large sequencing datasets.

The method is implemented in Python programming language and can be run on low-end Linux server and/or on laptop computers.

## Installation

PhenotypeSeeker supports Linux operating systems.

To install PhenotypeSeeker, open the command-line and type in the following commands:
```
sudo apt-get install git
git clone http://github.com/bioinfo-ut/PhenotypeSeeker.git
cd PhenotypeSeeker
sh install.sh

phenotypeseeker --version
```

### Install locally without sudo privileges

Python, python-pip, python-dev and git must be installed to install PhenotypeSeeker locally without sudo rights!

To install PhenotypeSeeker locally, open the command-line and type in the following commands:
```
git clone http://github.com/bioinfo-ut/PhenotypeSeeker.git
cd PhenotypeSeeker
sh local_install.sh

~/.local/bin/phenotypeseeker --version
```

That's it with PhenotypeSeeker!

## Usage

We provide the example analysis scripts for 'PhenotypeSeeker modeling' and 'PhenotypeSeeker prediction'. For detailed information about example analysis, please refer to the PhenotypeSeeker/example/README.md

To run the example analysis of 'PhenotypeSeeker modeling' and 'PhenotypeSeeker prediction' execute the example scripts in PhenotypeSeeker/example directory:
```
cd example
./test_PS_modeling.sh
./test_PS_prediction.sh
```


## Contact
If you require further assistance please contact us by sending email to
erki.aun@ut.ee
