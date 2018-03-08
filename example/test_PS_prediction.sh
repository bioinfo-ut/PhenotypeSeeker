#!/bin/bash
# Script for testing phenotype prediction with PhenotypeSeeker, using
# the previously generated model on Clostridium difficile genomes 
# with unknown phenotype.

mkdir PS_prediction_example_analysis
cd PS_prediction_example_analysis

# Download C. difficile FASTA files and inputfiles for PhenotypeSeeker prediction
echo Downloading the folder with example files ...
echo
wget http://bioinfo.ut.ee/PhenotypeSeeker/PS_prediction_example_files.tar.gz

#Unpack the downloaded folder
echo Unpacking the folder ...
echo
tar -zxvf PS_prediction_example_files.tar.gz

#Launch the "PhenotypeSeeker prediction"
echo
echo "Launching the phenotypeseeker prediction:"
echo "phenotypeseeker prediction PS_prediction_example_files/inputfile1 PS_prediction_example_files/inputfile2"
echo

phenotypeseeker prediction PS_prediction_example_files/inputfile1 PS_prediction_example_files/inputfile2

echo "Finished!"
