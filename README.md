# PhenotypeSeeker
Identify phenotype-specific k-mers and predict phenotype using sequenced bacterial strains
## Introduction
## Installation
PhenotypeSeeker currently supports Linux operating systems.

With python2.7 and pip2.7 installed, easiest way to install PhenotypeSeeker is to type
```
pip2.7 install PhenotypeSeeker
```
on the command-line (may need sudo command).

Or download the latest version of PhenotypeSeeker from the GitHub repository:
```
git clone https://github.com/bioinfo-ut/PhenotypeSeeker.git
```
Then, change to PhenotypeSeeker main directory and run the installer (may also need sudo command):
```
cd PhenotypeSeeker
pip2.7 install .
```
This will install PhenotypeSeeker and the required Python packages. PhenotypeSeeker executable is added to $PATH.

Test the installation, by typing:
```
phenotypeseeker --version
```
That's it with PhenotypeSeeker! 
You also need to install prerequisites.

### Prerequisites

To install and run PhenotypeSeeker, python2.7 and pip2.7 need to be installed on your system.

To install pip2.7 with python2.7 and curl installed, open the command-line and type in:
```
curl https://bootstrap.pypa.io/ez_setup.py -o - | sudo python2.7
sudo easy_install pip2.7
```

You also need to install yourself:

Genometester4
```
git clone https://github.com/bioinfo-ut/GenomeTester4.git
cd src
make
```
mash

You can download the pre-compiled mash binary from 
```
https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar
```
and then put it directory specified in $PATH variable (e.g. /usr/bin).

Or you can clone the source repository from Github
```
git clone https://github.com/marbl/Mash.git
```
And install it following the instructions specified in
```
https://github.com/marbl/Mash/blob/master/INSTALL.txt
```
## Example

To create the phenotype prediction model with PhenotypeSeeker, open the command-line and type in:
```
phenotypeseeker modeling data.pheno 
```
Where "data.pheno" is an input text file containing tab separated lists of sampleID's, sample FastA/FastQ file addresses and sample phenotype values (one or more column).

Example of "data.pheno" inputfile with binary resistance phenotypes for 5 antibiotics:
```

SampleID        Address Ciprofloxacin   Imipenem        Meropenem       Tobramycin Colistin         
PA2_D4  /storage8/erkia/data/pseudomonas_genomes/PA2_D4.fasta   0       NA	0 	NA	0
PA2_F7  /storage8/erkia/pseudomonas_genomes/PA2_F7.fasta 	1	1	1	0	1
1D4     /storage8/erkia/pseudomonas_genomes/1D4.fasta           1	1 	1 	NA	1
6G2     /storage8/erkia/pseudomonas_genomes/6G2.fasta           1	1	1       1       1
PA4_A1  /storage8/erkia/pseudomonas_genomes/PA4_A1.fasta 	0	0	0       0       0
6F4     /storage8/erkia/pseudomonas_genomes/6F4.fasta           0       NA	0 	NA	0
5H2     /storage8/erkia/pseudomonas_genomes/5H2.fasta           1	1 	1 	NA	1
PA2_C6  /storage8/erkia/pseudomonas_genomes/PA2_C6.fasta 	1	0	1	1	0
PA3_A4  /storage8/erkia/pseudomonas_genomes/PA3_A4.fasta 	0	0	0	0       0
PA3_A10 /storage8/erkia/pseudomonas_genomes/PA3_A10.fasta       0	1	0	NA	0
```
Example of "data.pheno" inputfile with continuous resistance phenotypes for 5 antibiotics:
```
SampleID        Address Ciprofloxacin   Imipenem        Meropenem	Tobramycin Colistin
PA3_E8  /storage8/erkia/data/pseudomonas_genomes/PA3_E8.fasta   0.1250  2.0000  3.0000  4.0000  2.0000
7D2     /storage8/erkia/pseudomonas_genomes/7D2.fasta           2.0000  1.5000  3.0000  12.0000 0.5000
PA3_H11 /storage8/erkia/pseudomonas_genomes/PA3_H11.fasta       0.0940  4.0000  4.0000  8.0000  1.5000
PA2_H12 /storage8/erkia/pseudomonas_genomes/PA2_H12.fasta       0.0940  32.0000 3.0000  0.5000  2.0000
5E7     /storage8/erkia/pseudomonas_genomes/5E7.fasta           0.1250  3.0000  3.0000  1.0000  256.0000
PA3_H4  /storage8/erkia/pseudomonas_genomes/PA3_H4.fasta        0.1250  6.0000  0.0940  1.0000  1.0000
PA2_C4  /storage8/erkia/pseudomonas_genomes/PA2_C4.fasta        32.0000 32.0000 0.1250  0.2500  2.0000
PA2_A3  /storage8/erkia/pseudomonas_genomes/PA2_A3.fasta        16.0000 2.0000  0.2500  0.2500  NA
PA3_A9  /storage8/erkia/pseudomonas_genomes/PA3_A9.fasta        32.0000 4.0000  6.0000  0.5000  NA
6F7     /storage8/erkia/pseudomonas_genomes/6F7.fasta           32.0000 8.0000  32.0000 2.0000  1.0000
```
To limit the analysis on a selection of phenotypes present in "data.pheno", specify the columns with "--mpheno" option. For example, to analyse only the phenotypes of ciprofloxacin and meropenem from example "data.pheno" file, type in:
```
phenotypeseeker modeling data.pheno --mpheno 1 3
```
## Contact
