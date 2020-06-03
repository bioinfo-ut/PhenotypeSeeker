# PhenotypeSeeker
Identify phenotype-specific k-mers and predict phenotype using sequenced bacterial strains

## Introduction

PhenotypeSeeker identifies phenotype-specific k-mers, generates phenotype prediction model and predicts the phenotype from sequencing data. 

PhenotypeSeeker consist of two subprograms: 'PhenotypeSeeker modeling' and 'PhenotypeSeeker prediction', which both take either assembled contigs or raw-read data as an input. 


PhenotypeSeeker uses statistical model that can be trained automatically on isolates with known phenotype. The prediction of phenotypes takes less than a second per isolate if assembled genomes are used and less than a minute per isolate if raw sequencing data are used. Therefore, PhenotypeSeeker is well suited for predicting phenotypes from large sequencing datasets.

The method is implemented in Python programming language and can be run on low-end Linux server and/or on laptop computers.

## Installation

PhenotypeSeeker supports Linux operating systems.

To install the PhenotypeSeeker, open the command-line and type in the following commands:
```
sudo apt-get install git
git clone http://github.com/bioinfo-ut/PhenotypeSeeker.git
cd PhenotypeSeeker
sh install.sh

phenotypeseeker --version
```

### Install locally without sudo privileges

Python3, python3-pip, python3-dev and git must be installed to install PhenotypeSeeker locally without sudo rights!

To install the PhenotypeSeeker locally, open the command-line and type in the following commands:
```
git clone http://github.com/bioinfo-ut/PhenotypeSeeker.git
cd PhenotypeSeeker
sh local_install.sh

~/.local/bin/phenotypeseeker --version
```

That's it with PhenotypeSeeker!

## Usage

The following examples and instructions for PhenotypeSeeker are based on the analysis of example dataset. For more detailed examples and instructions, please refer to the user_manual.md

### Automated example analysis scripts

We provide the automated example analysis scripts for 'PhenotypeSeeker modeling' and 'PhenotypeSeeker prediction'. For detailed information about example analysis, please refer to the PhenotypeSeeker/example/README.md

To run the automated example analysis of 'PhenotypeSeeker modeling' and 'PhenotypeSeeker prediction' execute the example scripts in PhenotypeSeeker/example directory:
```
cd example
./test_PS_modeling.sh
./test_PS_prediction.sh
```

The C.difficile genomes used in the example analysis originate from the European Nucleotide Archive [EMBL:PRJEB11776 ((http://www.ebi.ac.uk/ena/data/view/PRJEB11776)] and the binary phenotypes of azithromycin resistance for these genomes are adapted from Drouin et al. (2016. Predictive computational phenotyping and biomarker discovery using reference-free genome comparisons. BMC Genomics [Internet]. 17(1):754. Available from: http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2889-6).

### Launching the "PhenotypeSeeker modeling" (step by step guide)

First download and unpack the folder containing C. difficile FASTA files and inputfile for 'PhenotypeSeeker modeling' example analysis:

```
wget http://bioinfo.ut.ee/PhenotypeSeeker/PS_modeling_example_files.tar.gz
tar -zxvf PS_modeling_example_files.tar.gz
```

To create the phenotype prediction model with PhenotypeSeeker, type in:

```
phenotypeseeker modeling PS_modeling_example_files/data.pheno 
```
Where "data.pheno" is an input text file containing tab separated lists of (1) sampleID's, (2) sample FastA/FastQ file addresses and (3) sample phenotype values (one or more column). The first line of the input file is header specifying the names of the columns. The names of the phenotype columns are used in the names of the output files of the corresponding phenotype analyses.

Head of "data.pheno" inputfile with binary resistance phenotype for azithromycin:
```
ID	Addresses	Azithromycin
VL_0464	PS_modeling_example_files/VL_0464.fasta	0
VL_0456	PS_modeling_example_files/VL_0456.fasta	1
VL_0453	PS_modeling_example_files/VL_0453.fasta	1
VL_0442	PS_modeling_example_files/VL_0442.fasta	0
```

#### Outputfiles of PhenotypeSeeker modeling

The created **"chi-squared_test_results_Azithromycin.txt"** contains the statistical test results for every k-mer. More specifically, the columns in the file represent (1) the tested k-mer sequences, (2) the chi-squared statistic values, (3) the pvalues of chi-squared statistics, (4) the numbers of samples with the specific k-mers and (5) the names of the samples with the specific k-mers.

Head of **"chi-squared_test_results_Azithromycin.txt"** file:
```
GATAGAACTATTAAAC	6.05	4.86E-02	5	| VL_0198 VL_0456 VL_0088 VL_0137 VL_0302
AAATAAACTTACCTAT	19.46	5.96E-05	18	| VL_0064 VL_0040 VL_0453 VL_0065 VL_0091 VL_0303 VL_0233 VL_0073 VL_0464 VL_0098 VL_0159 VL_0442 VL_0377 VL_0288 VL_0047 VL_0346 VL_0252 VL_0296
ACTTGAGTATGCTATA	10.83	4.44E-03	8	| VL_0004 VL_0378 VL_0040 VL_0198 VL_0456 VL_0088 VL_0137 VL_0302
ATAATGGGTCCATTTA	22.43	1.35E-05	13	| VL_0004 VL_0216 VL_0368 VL_0378 VL_0040 VL_0112 VL_0276 VL_0369 VL_0198 VL_0456 VL_0088 VL_0137 VL_0302
```
The file **"k-mers_filtered_by_pvalue_Azithromycin.txt"** contains the subset of k-mers having the p-value of chi-squared test lower than 0.05 (default of "--pvalue" option).

The regression model is outputted in **"log_reg_model_Azithromycin.pkl"**

The **"summary_of_log_reg_analysis_Azithromycin.txt"** contains the information about conducted regression analysis. For example the regularisation parameter choosen, the actual vs predicted phenotypes of test set samples and the model-evaluation metrics.

The created **"k-mers_and_coefficients_in_log_reg_model_Azithromycin.txt"** contains the regression model coefficients for every k-mer used in the model as parameter. More specifically, the columns in the file represent (1) the k-mers used in regression model as parameter, (2) the regression model coefficients of k-mers, (3) the numbers of samples with the specific k-mers and (4) the names of the samples with the specific k-mers.

### Launching the "PhenotypeSeeker prediction" (step by step guide)

"PhenotypeSeeker prediction" predicts the phenotypes of input samples using the model previously created with "PhenotypeSeeker modeling". 

Therefore the "PhenotypeSeeker modeling" outputfiles "log_reg_model_Azithromycin.pkl" and "k-mers_and_coefficients_in_log_reg_model_Azithromycin.txt" are essential to run the "PhenotypeSeeker prediction".

Download and unpack the folder containing C. difficile FASTA files and inputfiles for 'PhenotypeSeeker prediction' example analysis:
```
wget http://bioinfo.ut.ee/PhenotypeSeeker/PS_prediction_example_files.tar.gz
tar -zxvf PS_prediction_example_files.tar.gz
```

To predict the phenotypes of samples under study, type in:
```
phenotypeseeker prediction PS_prediction_example_files/inputfile1 PS_prediction_example_files/inputfile2
```
Where: 

"inputfile1" is a text file containing tab separated lists of (1) sampleID's and (2) sample FastA/FastQ file addresses.

Head of "inputfile1":
```
VL_0068	PS_prediction_example_files/VL_0068.fasta
VL_0052	PS_prediction_example_files/VL_0052.fasta
VL_0155	PS_prediction_example_files/VL_0155.fasta
VL_0145	PS_prediction_example_files/VL_0145.fasta
```

"inputfile2" is a text file containing the tab separated list of (1) the name of the phenotype to predict, (2) corresponding model ("log_reg_model_Azithromycin.pkl") address and (3) corresponding k-mer list ("k-mers_and_coefficients_in_log_reg_model_Azithromycin.txt") address.
```
Azithromycin	PS_prediction_example_files/log_reg_model_Azithromycin.pkl	PS_prediction_example_files/k-mers_and_coefficients_in_log_reg_model_Azithromycin.txt
```

#### Outputfiles of PhenotypeSeeker prediction

After launcing the "PhenotypeSeeker prediction" it starts counting the k-mers from input samples, followed by the detection of presence or absence of the model specific k-mers in each sample.

The data of model specific k-mers presence or absence in each sample is used to predict the phenotypes for samples.

The results of predictions are saved into **"predictions_Azithromycin.txt"**.

Head of **"predictions_Azithromycin.txt"**
```
Sample_ID	predicted_phenotype	probability_for_predicted_class
VL_0068	1	1.0
VL_0052	0	0.53
VL_0155	0	0.53
VL_0145	1	1.0
```
## PhenotypeSeeker web tool

PhenotypeSeeker web tool enabels to use phenotype predictive models pre-trained on large datasets.

[Click here](http://bioinfo.ut.ee/PhenotypeSeeker) to visit PhenotypeSeeker webpage and use our web-tool for phenotype prediction from bacterial genome data.

## Contact
If you require further assistance please contact us by sending email to
erki.aun@ut.ee
