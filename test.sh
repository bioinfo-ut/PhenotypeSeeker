#!/bin/bash

# PhenotypeSeeker test commands
mkdir PS_testing
cd PS_testing

# Download C. difficile FASTA files and inputfile for PhenotypeSeeker
echo Downloading the folder with example files ...
echo
wget http://bioinfo.ut.ee/PhenotypeSeeker/PS_modeling_test_files.tar.gz

#Unpack the downloaded folder
echo Unpacking the folder ...
echo
tar -zxvf PS_modeling_example_files.tar.gz

#Launch the "PhenotypeSeeker modeling"
echo
echo "Launching the phenotypeseeker modeling:"
echo "phenotypeseeker modeling PS_modeling_example_files/data.pheno"
echo

phenotypeseeker modeling PS_modeling_example_files/data.pheno 