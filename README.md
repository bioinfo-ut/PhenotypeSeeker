# PhenotypeSeeker
Identify phenotype-specific k-mers and predict phenotype using sequenced bacterial strains
## Introduction
## Installation
PhenotypeSeeker currently supports Linux operating systems and needs python2.7 and pip2.7 to be installed, 
#### Install with pip
Easiest way to install PhenotypeSeeker is to type
```
pip2.7 install PhenotypeSeeker
```
on the command-line (may need sudo command).
#### Install from source
Download the PhenotypeSeeker from the GitHub repository:
```
git clone https://github.com/bioinfo-ut/PhenotypeSeeker.git
```
Then, change to PhenotypeSeeker main directory and run the installer (may also need sudo command):
```
cd PhenotypeSeeker
pip2.7 install .
```
#### Test the installation
This will install PhenotypeSeeker and the required Python packages. PhenotypeSeeker executable is added to $PATH.

Test the installation, by typing:
```
phenotypeseeker --version
```
That's it with PhenotypeSeeker! 
You also need to install prerequisites.

### Prerequisites

To install and run PhenotypeSeeker, python2.7, pip2.7, Genometester4 and mash need to be installed on your system.

#### Pip2.7
To install pip2.7 with python2.7 and curl installed, open the command-line and type in:
```
curl https://bootstrap.pypa.io/ez_setup.py -o - | sudo python2.7
sudo easy_install pip2.7
```

You also need to install yourself:

#### Genometester4
```
git clone https://github.com/bioinfo-ut/GenomeTester4.git
cd src
make
```
#### mash

You can download the pre-compiled mash binary from 
```
https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar
```
and then put it into the directory specified in $PATH variable (e.g. /usr/bin).

Or you can clone the source repository from Github
```
git clone https://github.com/marbl/Mash.git
```
and install it following the instructions specified in
```
https://github.com/marbl/Mash/blob/master/INSTALL.txt
```
## Examples

### Launching the PhenotypeSeeker

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
To limit the analysis on a selection of phenotypes present in "data.pheno", specify the phenotype columns with "--mpheno" option. For example, to analyse only the phenotypes of ciprofloxacin and meropenem from example "data.pheno" file, type in:
```
phenotypeseeker modeling data.pheno --mpheno 1 3
```
### Outputfiles of PhenotypeSeeker
After launcing the PhenotypeSeeker it starts counting the k-mers from input samples. 

To save the temporary k-mer lists, the temporary repository **"K-mer_lists"** is created.

The process of k-mer counting results in the **"k-mer_matrix.txt"** file, where the number of each counted k-mer (rows) in each sample (columns) is recorded.

Then, if "-w" (weights) option is choosed it estimates the distances of samples with program mash, which output is saved in **"mash_distances.mat"**. 

Those distances are used to calculate GSC (Gerstein, Sonnhammer and Chothia) weights, which are considered in forthcoming analysis to correct for clonal population structure.

Then the **"chi-squared_test_results.txt"** or **"t-test_results.txt"** file is created depending on which type of phenotype (binary or continuous accordingly) were used. 

The columns in "chi-squared_test_results.txt" represent (1) the tested k-mer sequences, (2) the chi-squared statistic values, (3) the pvalues of chi-squared statistics, (4) the numbers of samples with the specific k-mers and (5) the names of the samples with the specific k-mers.

Example of **"chi-squared_test_results.txt"** file:
```
ACCTCTGGGTGGCGAA	0.18	9.13E-01	3	| PA2_F7 PA2_G7 PA3_F7
AGCGCTACTCTTGATC	0.07	9.67E-01	8	| PA4_H4 PA3_G6 PA3_G4 PA3_A5 PA4_D2 PA3_G3 PA3_G9 PA3_G8
AAACTTTCGACACAGC	0.18	9.13E-01	3	| PA2_A2 PA2_A1 PA3_D6
GACCTGTCCCTCTCCA	0.18	9.13E-01	3	| PA2_B10 PA2_H6 PA2_C10
GCTCAATCGCTAAAGA	2.42	2.98E-01	2	| PA3_F10 PA3_C1
AAGTCGCTGGATTTCG	4.73	9.39E-02	7	| 1C9 PA3_A8 6E7 PA3_F7 PA2_F12 PA3_C7 PA3_A2
GGCACCCCGTTGGCCA	3.96	1.38E-01	9	| 6E5 1P4 PA3_B8 PA3_C9 PA2_B8 5I6 PA3_C4 PA2_F9 PA2_F4
AGCGCTGCTGCGCGAT	0.02	9.92E-01	2	| PA4_B3 PA4_H2
```

The columns in "t-test_results.txt" represent (1) the tested k-mer sequences, (2) the Welch's t-test statistic values, (3) the pvalues of t-test statistics, (4) the mean phenotype values of samples missing the specific k-mers, (5) the mean phenotype values of samples having the specific k-mers, (6) the numbers of samples with the specific k-mers and (7) the names of the samples with the specific k-mers.

Example of "t-test_results.txt" file:
```
CGAGCTCCCGGCA	24.86	3.14E-62	5.0	-0.54	4	| PA3_F5 PA2_F2 PA2_E2 PA2_G1
ATTGTCAAGGCGC	-10.89	4.06E-11	-3.27	-0.38	3	| PA4_F4 PA4_B2 PA2_H6
TTGCGAAGATAAA	24.47	1.71E-61	5.0	-0.48	2	| PA2_G6 5I3
AGCAAGGTGAGTC	-11.48	1.20E-23	-3.0	-0.39	3	| PA2_C8 PA3_C8 PA3_E1
CCATTTTTGTTAC	-9.02	5.81E-10	-2.85	-0.38	4	| PA2_C2 PA2_A2 PA3_B5 PA2_A1
CGAACTATAACCC	24.47	1.69E-61	5.0	-0.48	2	| PA3_G8 PA3_G7
GTCGCCCTCTGGA	-8.2	3.86E-06	-3.03	-0.37	4	| PA4_F4 PA3_H10 PA4_F1 PA3_H11
AATGCCTACGAAG	-11.46	1.30E-23	-3.0	-0.4	2	| PA4_A4 PA4_A5
CTGTTGGATCCCA	6.88	2.79E-06	3.24	-0.68	13	| PA3_E8 PA2_D10 PA3_D1 PA3_C10 PA3_C3 PA3_D7 PA3_H7 PA3_D3 PA3_F9 PA2_E9 PA3_C5 PA2_H3 PA4_G1
CACAACCCCAACA	-11.48	1.20E-23	-3.0	-0.39	3	| PA2_C8 PA2_C9 PA2_C4
```

The following file **"k-mers_filtered_by_pvalue.txt"** contains the subset of k-mers having the p-value lower than the p-value cut-off specified by "--pvalue" (with -B or -FDR) option(s).

Follows the regression analysis with top "--n_kmers" lowest p-valued k-mers from **"k-mers_filtered_by_pvalue.txt"** file.

The regression model is outputted in **"log_reg_model.pkl"** or **"lin_reg_model.pkl"** file depending on which type of phenotype (binary or continuous accordingly) were used. 

The **"summary_of_(log\lin)_reg_analysis.txt"** contains the information about conducted regression analysis. For example the regularisation parameter choosen, the actual vs predicted phenotypes of test set samples and the model-evaluation metrics.

The columns in the following file **"k-mers_and_coefficients_in_(log/lin)_reg_model.txt"** represent (1) the k-mers used in regression model as parameter, (2) the regression model coefficients of k-mers, (3) the numbers of samples with the specific k-mers and (4) the names of the samples with the specific k-mers.

Example of **"k-mers_and_coefficients_in_(log/lin)_reg_model.txt"** file:
```
K-mer	coef._in_lin_reg_model	No._of_samples_with_k-mer	Samples_with_k-mer
CAAGGCACTTTCC	-2.68E-01	143	| PA2_F7 1D4 6G2 PA4_A1 6F4 5H2 PA2_C6 PA3_A4 PA3_A10 1P3 PA4_F2 PA4_F4 PA4_B2 PA2_C3 PA2_B4 6F2 6B5 PA3_B10 PA3_H10 PA2_E4 PA2_G7 PA4_A4 7P3 PA2_C8 PA3_C8 PA3_E1 PA2_D9 PA2_C2 PA2_B2 PA2_A2 PA2_D1 PA4_D4 PA2_B1 PA3_D10 PA3_F1 PA2_G8 PA2_A4 6G5 PA2_H12 PA2_H1 PA2_C9 PA4_G2 PA2_F8 PA4_A5 PA4_E2 PA4_F1 PA2_G12 PA3_G2 PA3_F6 PA3_E5 PA3_E10 PA3_B6 PA2_F6 PA2_B10 PA2_A10 PA2_A3 6G1 PA2_H6 PA3_F4 PA4_E4 PA4_B3 PA4_H4 PA4_H1 PA3_H11 PA4_E1 PA2_E1 PA2_A9 1D3b 6E9 6F3 PA3_H2 PA3_G6 PA3_G4 5C6 PA2_A1 PA2_A7 PA2_B3 PA2_D3 PA4_D1 PA2_F1 PA4_C1 PA2_G4 PA3_B3 PA2_E6 PA2_F3 1P4 6F7 PA3_F3 PA2_E8 PA3_C2 PA3_A6 PA3_A5 PA2_H4 PA2_C7 PA4_D2 PA3_F10 PA3_B8 PA4_H2 5D7 PA3_B4 PA2_D6 PA3_E8 PA3_E7 PA3_H9 PA2_D7 PA3_C9 PA2_B8 5F7 PA3_G10 5D4 PA4_C2 PA3_H5 PA3_D6 PA2_E7 5H6 PA3_A8 PA3_A11 PA3_A7 6E7 PA3_F8 6B9 PA3_H4 PA3_H3 PA3_G1 PA3_E4 PA3_G8 6A1 PA3_B1 PA3_A1 7D2 7D1 1D3a 1D6 5E7 PA2_C10 PA3_D5 PA2_B7 PA3_F7 PA2_F12 PA3_C7 PA3_A2 PA2_F4 PA4_G4
CTGCCGGGAACGA	0.00E+00	24	| PA2_B2 PA3_D8 5D4 6F9 PA3_D1 PA3_C3 PA3_D7 PA3_H7 6E7 PA3_G9 6B1 PA3_D3 PA2_G6 PA3_G8 PA3_C1 6A1 PA3_D2 PA3_C5 PA2_H3 PA2_F9 PA3_F7 PA2_F12 PA3_C7 PA3_A2
AACGCTTTCGGAC	0.00E+00	36	| 5B8 PA2_G3 7P1 PA2_F3 PA2_A8 PA2_E3 PA3_D8 PA3_B2 5C3 PA3_A8 6B7 6F9 PA3_D1 PA3_C10 PA3_C3 PA3_D7 PA3_B7 PA3_H7 6E7 6B9 PA3_G9 6B1 PA3_D3 PA2_G6 6A1 PA3_F9 PA3_D2 PA2_E9 PA3_G7 PA3_C5 PA2_H3 PA3_D5 PA3_F7 PA2_F12 PA3_C7 PA3_A2
CTCGGGGCGGAAA	0.00E+00	19	| PA3_D1 PA3_C3 PA3_D7 PA3_H7 6E7 PA3_G9 6B1 PA3_D3 PA2_F2 PA2_E2 PA3_D2 PA3_C5 PA2_H3 PA2_F9 PA3_F7 PA2_F12 PA3_C7 PA3_A2 PA2_G1
GCCTGTAGGTCGC	0.00E+00	69	| 1D4 6F4 PA3_A4 1P3 PA4_F4 PA4_B2 PA2_E4 PA2_B9 7P3 PA3_D10 5B8 PA2_G3 PA3_F1 PA2_F8 PA4_E2 PA2_H6 PA4_E4 PA2_F1 PA2_G4 PA2_A8 PA2_E3 PA3_B4 PA3_D8 PA4_C3 PA4_G3 PA4_D3 PA4_E3 PA4_F3 PA3_B2 5I6 5C3 PA2_D10 PA3_A8 6B7 PA3_D1 PA3_C10 PA3_C3 PA3_A7 PA3_D7 PA3_B7 PA3_H7 PA3_H8 6E7 PA3_G9 6B1 PA3_D3 PA3_F5 PA2_F2 PA2_E2 PA3_H3 PA3_G1 PA3_F9 PA3_D2 PA2_E9 PA3_G7 PA3_C5 PA2_H3 PA3_B1 PA3_A1 7D2 7D1 5I3 PA3_D5 PA3_F7 PA2_F12 PA3_C7 PA3_A2 PA2_G1 PA4_G1
TCGCGGTCTACGA	0.00E+00	42	| PA2_C9 6E5 PA3_B4 5D4 5I6 PA2_D10 6B7 PA3_D1 PA3_C10 PA3_C4 PA3_A11 PA3_C3 PA3_A7 PA3_D7 PA3_H7 6E7 PA3_F8 6B9 PA3_G9 6B1 PA3_D3 PA3_F5 PA2_F2 PA2_E2 PA2_G6 PA3_H3 PA3_G1 PA3_G8 PA3_C1 6A1 PA3_F9 PA3_D2 PA2_E9 PA3_C5 PA2_H3 5I3 PA2_B7 PA3_F7 PA3_C7 PA3_A2 PA2_F4 PA2_G1
CCTTGACCGAACG	7.17E-01	35	| 1D4 PA2_B9 PA2_B2 PA4_A5 1C9 PA4_C3 PA4_G3 PA4_D3 PA4_E3 PA4_F3 5D4 5I6 PA3_D1 PA3_C3 PA3_D7 PA3_H7 PA3_G9 6B1 PA3_D3 PA3_H4 PA2_G6 PA3_G8 6A1 PA3_D2 PA3_C5 PA2_H3 PA3_B1 PA3_A1 7D2 7D1 5I3 PA2_F9 PA3_F7 PA3_C7 PA3_A2
```


## Contact
