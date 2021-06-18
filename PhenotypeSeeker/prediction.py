#!/usr/bin/env python3

__author__ = "Erki Aun"
__version__ = "0.7.3"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

from subprocess import call
import math
import warnings
warnings.showwarning = lambda *args, **kwargs: None

import pkg_resources
pkg_resources.require("numpy==1.18.1", "scikit-learn==0.22.1")

import numpy as np
from sklearn.externals import joblib

def get_kmers(phenotypes_to_predict):
    call(["mkdir", "-p", "K-mer_lists"])
    for phenotype in phenotypes_to_predict:
        k_mer_list_2_allocate_from = phenotypes_to_predict[phenotype][1]
        outname = "K-mer_lists/k-mers_" + phenotype + ".txt"
        with open(outname, "w+") as f1:
            call(
            	["cut -f1 "  + str(k_mer_list_2_allocate_from) 
            	+ " | tail -n +2"], 
            	shell=True, stdout=f1
            	)

def format_kmer_db(phenotypes_to_predict):
    for phenotype in phenotypes_to_predict:
        with open("K-mer_lists/k-mers_" + phenotype + ".txt") as f1:
            with open(
            	    "K-mer_lists/k-mer_db_" + phenotype + ".txt", "w+"
            	    ) as f2:
                counter = 1 
                for line in f1:
                    line = str(counter) + "\t1\t" + line
                    f2.write(line)
                    counter += 1

def kmer_filtering_by_freq_cutoff_in_sample(
        samples_info, min_freq, phenotypes_to_predict
        ):
    for phenotype in phenotypes_to_predict:
        for ID in samples_info:
            with open(
                    "K-mer_lists/" + ID + "_k-mer_counts_"
                    + phenotype  + ".txt"
                    ) as f1:
                with open(
                        "K-mer_lists/" + ID + "_k-mer_counts_filtered_"
                        + phenotype + ".txt", "w+"
                        ) as f2:
                    for line in f1:
                        if "#TextDatabase" in line:
                            continue
                        list1 = line.strip().split()
                        if len(list1) == 2 or list1[2] >= min_freq:
                            f2.write(line)
                        else:
                            list1[2] = "0"
                            f2.write("\t".join(list1)+"\n")

def map_samples_prediction(samples_info, phenotypes_to_predict):
    # Takes k-mers of model as feature space and maps input samples 
    # k-mer lists to that feature space. A vector of k-mers frequency 
    # information is created for every sample.
    for phenotype in phenotypes_to_predict:
        for ID in samples_info:
            call(
            	["gmer_counter -db K-mer_lists/k-mer_db_" + phenotype 
            	+ ".txt " +  samples_info[ID][0] + " > K-mer_lists/" + ID 
            	+ "_k-mer_counts_" + phenotype  + ".txt"], 
            	shell=True
            	)

def parse_prediction_input_file1(inputfilename):
    # Parses info from tabulated input file into    samples directory.
    # Stores the order of samples in "samples_order" list.
    samples = {}
    samples_order = []
    n_o_s = 0
    with open(inputfilename) as f1:
        for line in f1:
            if line == "\n":
                break
            line = line.strip()
            list1 = line.split()
            samples[list1[0]] = list1[1:]
            samples_order.append(list1[0])
            n_o_s += 1
    return(samples, samples_order, n_o_s)

def parse_prediction_input_file2(inputfilename):
    # Parses info from tabulated input file into phenotypes_to_predict 
    # directory.
    phenotypes_to_predict = {}
    with open(inputfilename) as f1:
        for line in f1:
            if line == "\n":
                break
            line = line.strip()
            list1 = line.split()
            phenotypes_to_predict[list1[0]] = list1[1:]
    return(phenotypes_to_predict)

def vectors_to_matrix_prediction(samples_order, phenotypes_to_predict):
    # Takes all vectors with k-mer frequency information and inserts 
    # them into matrix of dimensions "number of samples" x "number of 
    # k-mers (features).
    for phenotype in phenotypes_to_predict:
        kmer_matrix = open("K-mer_lists/k-mer_matrix_" + phenotype  + ".txt", "w")
        kmer_list_files = [
            "K-mer_lists/" + item + "_k-mer_counts_filtered_" + phenotype
            + ".txt" for item in samples_order
            ]
        for line in zip(*[open(item) for item in kmer_list_files]):
            kmer_matrix.write('\t'.join([j.split()[2].strip() for j in line]) + "\n")

def predict(samples_order, phenotypes_to_predict):
    # Generating a binary k-mer presence matrix
    for phenotype in phenotypes_to_predict:
        kmers_presence_matrix = []
        with open("K-mer_lists/k-mer_matrix_" + phenotype  + ".txt") as f1:
            for line in f1:
                kmers_presence_matrix.append([0 if int(x) == 0 else 1 for x in line.split()])
        kmers_presence_matrix = np.array(kmers_presence_matrix).transpose()
        
        #Loading regression model
        model = joblib.load(phenotypes_to_predict[phenotype][0])
        predictions = model.predict(kmers_presence_matrix)
        
        with open("predictions_" + phenotype + ".txt", "w+") as f1:
            model_name_short = phenotypes_to_predict[phenotype][0].split("_model")[0].split("/")[-1]
            if model_name_short in ("log_reg", "NB", "RF", "SVM", "XGBC"):
                predict_proba = model.predict_proba(kmers_presence_matrix)
                f1.write("Sample_ID\tpredicted_phenotype\t" \
                    "probability_for_predicted_class\n")
                for ID, prediction, proba in zip(
                        samples_order, predictions, predict_proba
                        ): 
                    f1.write(ID + "\t" + str(prediction) 
                        + "\t" + str(round(proba[1], 2))  + "\n")
            else:
                f1.write("Sample_ID\tpredicted_phenotype\n")
                for ID, prediction in zip(samples_order, predictions):
                    f1.write(ID + "\t" + str(prediction) + "\n")

def prediction(args):
    samples, samples_order, n_o_s = parse_prediction_input_file1(
        args.inputfile1
        )
    phenotypes_to_predict = parse_prediction_input_file2(args.inputfile2)
    
    get_kmers(phenotypes_to_predict)
    format_kmer_db(phenotypes_to_predict)
    map_samples_prediction(samples, phenotypes_to_predict)
    kmer_filtering_by_freq_cutoff_in_sample(
        samples, args.c, phenotypes_to_predict)
    vectors_to_matrix_prediction(samples_order, phenotypes_to_predict
        )
    predict(samples_order, phenotypes_to_predict)            
