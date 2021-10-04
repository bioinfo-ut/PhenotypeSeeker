#!/usr/bin/env python3

__author__ = "Erki Aun"
__version__ = "0.7.3"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

from functools import partial
from subprocess import call
from collections import OrderedDict
import math

import joblib
import pkg_resources
import numpy as np
import sys
import warnings

pkg_resources.require("numpy==1.18.1", "scikit-learn==0.22.1")
warnings.showwarning = lambda *args, **kwargs: None

class Input():
    
    samples = OrderedDict()
    phenos = OrderedDict()

    @classmethod
    def get_samples(cls, inputfile):
        # Parses info from tabulated input file into    samples directory.
        # Stores the order of samples in "samples_order" list.
        with open(inputfile) as inp:
            for line in inp:
                if line.strip():
                    sample_name = line.split()[0]
                    cls.samples[sample_name] = Samples.from_inputfile(line)

    @classmethod
    def get_phenos(cls, inputfile):
        # Parses info from tabulated input file into    samples directory.
        # Stores the order of samples in "samples_order" list.
        n_o_s = 0
        with open(inputfile) as inp:
            for line in inp:
                if line.strip():
                    pheno = line.split()[0]
                    cls.phenos[pheno] = Phenotypes.from_inputfile(line)

class Samples():

    no_samples = 0

    def __init__(self, name, address):
        self.name = name
        self.address = address

        Samples.no_samples += 1

    @classmethod
    def from_inputfile(cls, line):
        name, address = \
            line.split()[0], line.split()[1]
        return cls(name, address)

    def map_samples(self, pheno):
        # Takes k-mers of model as feature space and maps input samples 
        # k-mer lists to that feature space. A vector of k-mers frequency 
        # information is created for every sample.
        call(
            ["gmer_counter -db K-mer_lists/k-mer_db_" + pheno
            + ".txt " + self.address + " > K-mer_lists/" + self.name 
            + "_k-mer_counts_" + pheno  + ".txt"], shell=True
            )

    def kmer_counts(self, pheno):
        with open(
                "K-mer_lists/" + self.name +
                "_k-mer_counts_"+ pheno  + ".txt"
            ) as counts:
            with open(
                    "K-mer_lists/" + self.name + "_k-mer_counts_filtered_"
                    + pheno + ".txt", "w+"
                    ) as counts_filtered:
                counts.readline()
                for line in counts:
                    splitted = line.split()
                    if splitted[2] >= pheno.cutoff:
                        splitted[2] = "1"
                        counts_filtered.write("\t".join(splitted)+"\n")
                    else:
                        splitted[2] = "0"
                        counts_filtered.write("\t".join(splitted)+"\n")

class Phenotypes():

    cutoff = 1
    no_phenotypes = 0

    def __init__(self, name, model, kmers, pca, pca_model, scaler):
        self.name = name
        self.model = model
        self.kmers = kmers
        self.pca = pca
        self.pca_model = pca_model
        self.scaler = scaler
        self.kmer_matrix = None

        Phenotypes.no_phenotypes += 1

    @classmethod
    def from_inputfile(cls, line):
        name, model_adre = line.split()[0], line.split()[1] 
        model_pkg = joblib.load(model_adre)
        model = model_pkg['model']
        kmers = model_pkg['kmers']

        pca = False
        pca_model = None
        scaler = None
        if model_pkg['pca']:
            pca_model = model_pkg['pca_model']
            scaler = model_pkg['scaler']
            pca = True
        return cls(name, model, kmers, pca, pca_model, scaler)

    def set_kmer_db(self):
        with open("K-mer_lists/k-mer_db_" + self.name + ".txt", "w+") as db:
            for kmer in self.kmers:
                db.write(f"{kmer}\t1\t{kmer}")

    def get_inp_matrix(self):
        # Takes all vectors with k-mer frequency information and inserts 
        # them into matrix of dimensions "number of samples" x "number of 
        # k-mers (features).
        kmer_counts = [
            "K-mer_lists/" + sample + "_k-mer_counts_filtered_" + self.name
            + ".txt" for sample in Input.samples.keys()
            ]
        for line in zip(*[open(counts) for counts in kmer_counts]):
            self.kmer_matrix = np.array([j.split()[2].strip() for j in line]).transpose()
        if self.pca:
            km

    def predict(self):

        #Loading regression model
        model = joblib.load(self.model)
        predictions = model.predict(self.kmer_matrix)
        
        with open("predictions_" + phenotype + ".txt", "w+") as out:
            if self.pred_type == "binary":
                predict_proba = model.predict_proba(self.kmer_matrix)
                f1.write("Sample_ID\tpredicted_phenotype\t" \
                    "probability_for_predicted_class\n")
                for sample, prediction, proba in zip(
                        Input.samples.keys(), predictions, predict_proba
                        ): 
                    f1.write(f"{sample}\t{str(prediction)}\t{str(round(proba[1], 2))}\n")
            else:
                f1.write("Sample_ID\tpredicted_phenotype\n")
                for sample, prediction in zip(Input.samples.keys(), predictions):
                    f1.write(sample + "\t" + str(prediction) + "\n")

def prediction(args):
    sys.stderr.write("\x1b[1;1;101m######                   PhenotypeSeeker                   ######\x1b[0m\n")
    sys.stderr.write("\x1b[1;1;101m######                      prediction                       ######\x1b[0m\n\n")

    Phenotypes.cufoff = args.c
    call(["mkdir", "-p", "K-mer_lists"])
    Input.get_samples(args.inputfile1)
    Input.get_phenos(args.inputfile2)

    
    for pheno in Input.phenos.values():
        sys.stderr.write(f"\x1b[1;32mPredicting the phenotypes for {pheno.name}.\x1b[0m\n")
        pheno.set_kmer_db()
        map(partial(
            lambda x: x.map_samples(), pheno.name),
            Input.samples.values())
        map(partial(
            lambda x: x.kmer_counts(), pheno.name),
            Input.samples.values())
        pheno.get_inp_matrix()
        pheno.predict()

    sys.stderr.write("\n\x1b[1;1;101m######          PhenotypeSeeker prediction finished          ######\x1b[0m\n")         
