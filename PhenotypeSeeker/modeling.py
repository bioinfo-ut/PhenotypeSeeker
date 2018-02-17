#!/usr/bin/python2.7

__author__ = "Erki Aun"
__version__ = "1.0"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

from itertools import permutations
from subprocess import call, Popen, PIPE
import math
import sys
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from cogent import LoadTree
from cogent.align.weights.methods import GSC
from scipy import stats
from sklearn.externals import joblib
from sklearn.linear_model import Lasso, LogisticRegression, Ridge
from sklearn.metrics import classification_report, r2_score, mean_squared_error
from sklearn.model_selection import GridSearchCV, train_test_split
import Bio
import sklearn.datasets
import numpy as np

class Printer():
    # Print things to stdout on one line dynamically.
    def __init__(self,data):
        sys.stderr.write("\r\x1b[K"+data.__str__())
        sys.stderr.flush()


def write_to_stderr_if(
        previousPercent, currentKmerNum, totalKmers, text, phenotype=""
        ):
    currentPercent = currentKmerNum/totalKmers*100    
    if int(currentPercent) > previousPercent:
        output = "\t" + phenotype + "%d%% of %d " % (
            currentPercent,totalKmers
            ) + text
        Printer(output)
        previousPercent = currentPercent
    currentKmerNum += 1
    return(previousPercent, currentKmerNum)

def parse_modeling_input_file(inputfilename):
    # Parses info from tabulated input file into samples directory.
    # Stores the order of samples in "samples_order" list.
    # Counts the number of samples and phenotypes and stores those
    # values in n_o_s and n_o_p variables, respectively.
    samples = {}
    samples_order = []
    n_o_s = 0
    headerline = False
    phenotype_scale = "binary"
    phenotypes = []
    with open(inputfilename) as f1:
        for line in f1:
            if line == "\n":
                break
            line = line.strip()
            list1 = line.split()
            if list1[0] == "ID":
                phenotypes = list1[2:]
                headerline = True
            else:
                for item in list1[2:]:
                    if item != "1" and item != "0" and item != "NA":
                        phenotype_scale = "continuous"
                samples[list1[0]] = list1[1:]
                samples_order.append(list1[0])
                n_o_s += 1
    n_o_p = len(list1[2:])
    return(
    	samples, samples_order, n_o_s, n_o_p, 
    	phenotype_scale, headerline, phenotypes
    	)

def kmer_list_generator(samples_info, kmer_length, freq, input_samples):
    # Makes "K-mer_lists" directory where all lists are stored.
    # Generates k-mer lists for every sample in sample_names variable 
    # (list or dict).
    sys.stderr.write("Generating the k-mer lists:\n")
    totalFiles = len(input_samples)
    currentSampleNum = 1
    call(["mkdir", "-p", "K-mer_lists"])
    for item in input_samples:
        out_name = "K-mer_lists/" + item + "_output.txt"
        call(
        	["glistmaker " + str(samples_info[item][0]) + " -o K-mer_lists/" 
        	+ item + " -w " + kmer_length + " -c " + freq], 
        	shell=True
        	)
        with open(out_name, "w+") as f2:
            call(
            	["glistquery", "K-mer_lists/" + item + "_" 
            	+ kmer_length + ".list"], 
            	stdout=f2
            	)
        output = "\t%d of %d lists generated." % (currentSampleNum,totalFiles)
        Printer(output)
        currentSampleNum += 1

def kmer_frequencies(samples):
    # Counts the k-mers presence frequencies in samples
    sys.stderr.write(
        "\nCounting the k-mers presence frequencies in samples:\n"
        )
    totalFiles = len(samples)
    currentSampleNum = 1
    dict_of_frequencies = {}
    for item in samples:
        with open("K-mer_lists/" + item + "_output.txt") as f1:
            for line in f1:
                if line.split()[0] not in dict_of_frequencies:
                    dict_of_frequencies[line.split()[0]] = 1
                else:
                    dict_of_frequencies[line.split()[0]] = dict_of_frequencies[
                    line.split()[0]
                    ] + 1
        output = "\t%d of %d samples counted for k-mers presence." % (
            currentSampleNum, totalFiles
            )
        Printer(output)
        currentSampleNum += 1
    return(dict_of_frequencies)

def kmer_filtering_by_frequency(dict_of_frequencies, min_freq, max_freq):
    #Filters k-mers by their frequency in samples.
    sys.stderr.write(
        "\nFiltering the k-mers based on their " \
        "presence frequencies in samples:\n"
        )
    kmers_passed = 0
    totalKmers = len(dict_of_frequencies)
    currentKmerNum = 1.0
    previousPercent = 0.0
    f1 = open("K-mer_lists/k-mers_filtered_by_freq.txt", "a")
    for key, value in dict_of_frequencies.iteritems():
        if int(value) >= int(min_freq) and int(value) <= int(max_freq):
            f1.write(key + "\n")
            kmers_passed += 1
        previousPercent, currentKmerNum = write_to_stderr_if(
            previousPercent, currentKmerNum, totalKmers, "k-mers filtered."
            )
    sys.stderr.write("\nBuilding the feature vector from filtered k-mers.\n")
    f1.close()
    return(kmers_passed) 

def map_samples_modeling(sample_names, kmer_length):
    # Takes k-mers, which passed frequency filtering as feature space 
    # and maps samples k-mer lists to that feature space. A vector of 
    # k-mers frequency information is created for every sample.
    sys.stderr.write("Mapping samples to the feature vector space:\n")
    totalFiles = len(sample_names)
    currentSampleNum = 1
    for item in sample_names:
        out_name = "K-mer_lists/"+ item + "_output2.txt"
        with open(out_name, "w+") as f1:
            call(
            	["glistquery", "K-mer_lists/" + item + "_" + kmer_length 
            	+ ".list", "-f", "K-mer_lists/k-mers_filtered_by_freq.txt"],
            	stdout=f1
            	)
        output = "\t%d of %d samples mapped." % (currentSampleNum,totalFiles)
        Printer(output)
        currentSampleNum += 1

def vectors_to_matrix_modeling(samples_order):
    # Takes all vectors with k-mer frequency information and inserts 
    # them into matrix of dimensions "number of samples" x "number of 
    # k-mers (features).
    sys.stderr.write("\nConverting the set of vectors into data matrix.\n")
    totalFiles = len(samples_order)
    currentSampleNum = 1
    call(["mv", "K-mer_lists/k-mers_filtered_by_freq.txt", "k-mer_matrix.txt"])
    for item in samples_order:
        with open("K-mer_lists/tmp1.txt", "w+") as f1:
            with open("K-mer_lists/tmp2.txt", "w+") as f2:
                call(
                	["cut", "-f", "2", "K-mer_lists/" + item + "_output2.txt"],
                	stdout=f2
                	)
                call(
                	["paste", "k-mer_matrix.txt", "K-mer_lists/tmp2.txt"], 
                	stdout=f1
                	)
                call(["mv", "K-mer_lists/tmp1.txt", "k-mer_matrix.txt"])
        output = "\t%d of %d vectors inserted into matrix." % (
            currentSampleNum,totalFiles
            )
        Printer(output)
        currentSampleNum += 1

def mash_caller(samples_info, freq):
    #Estimating phylogenetic distances between samples using mash
    sys.stderr.write("\nEstimating the Mash distances between samples...\n")
    mash_args = ["mash", "sketch", "-o", "reference", "-m", freq]
    for item in samples_info:
        mash_args.append(samples_info[item][0])
    process = Popen(mash_args, stderr=PIPE)
    for line in iter(process.stderr.readline, ''):
        Printer(line.strip())
    Printer("")
    with open("mash_distances.mat", "w+") as f1:
        call(["mash", "dist", "reference.msh", "reference.msh"], stdout=f1)

def mash_output_to_distance_matrix(samples_order, mash_distances):
    with open(mash_distances) as f1:
        with open("distances.mat", "w+") as f2:
            counter = 0
            f2.write(samples_order[counter])
            for line in f1:
                distance = line.split()[2]
                f2.write("\t" + distance)
                counter += 1
                if counter%len(samples_order) == 0:
                    if counter != len(samples_order)**2:
                        f2.write(
                        	"\n" + samples_order[counter/len(samples_order)]
                        	)

def distance_matrix_modifier(distance_matrix):
    # Modifies distance matrix to be suitable argument 
    # for Bio.Phylo.TreeConstruction._DistanceMatrix function
    distancematrix = []
    with open(distance_matrix) as f1:
        counter = 2
        for line in f1:
            line = line.strip()
            list1 = line.split()
            distancematrix.append(list1[1:counter])
            counter += 1
    for i in range(len(distancematrix)):
        for j in range(len(distancematrix[i])):
            distancematrix[i][j] = float(distancematrix[i][j])
    return(distancematrix)

def distance_matrix_to_phyloxml(names_order_in_dist_mat, distance_matrix):
    #Converting distance matrix to phyloxml
    dm = _DistanceMatrix(names_order_in_dist_mat, distance_matrix)
    tree_xml = DistanceTreeConstructor().nj(dm)
    with open("tree_xml.txt", "w+") as f1:
        Bio.Phylo.write(tree_xml, f1, "phyloxml")

def phyloxml_to_newick(phyloxml):
    #Converting phyloxml to newick
    with open("tree_newick.txt", "w+") as f1:
        Bio.Phylo.convert(phyloxml, "phyloxml", f1, "newick")

def newick_to_GSC_weights(newick_tree):
    # Calculating Gerstein Sonnhammer Coathia weights from Newick 
    # string. Returns dictionary where sample names are keys and GSC 
    # weights are values.
    sys.stderr.write("Calculating the Gerstein Sonnhammer Coathia " \
        "weights from mash distance matrix...")
    tree=LoadTree(newick_tree)
    weights=GSC(tree)
    for item in weights:
        weights[item] = 1 - weights[item]
    return(weights)

def weighted_t_test(
	    kmer_matrix, samples, samples_order, weights, number_of_phenotypes, 
	    phenotypes, kmers_to_analyse, phenotypes_to_analyze=False, FDR=False, 
        headerline=False
	    ):
    # Calculates weighted Welch t-tests results for every k-mer
    sys.stderr.write(
        "\nConducting the k-mer specific weighted Welch t-tests:\n"
        )
    pvalues_all_phenotypes = []
    nr_of_kmers_tested_all_phenotypes = []
    outputfiles = []
    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)
    for k in phenotypes_to_analyze:
        currentKmerNum = 1.0
        previousPercent = 0.0
        pvalues = []
        counter = 0
        NA = False
        with open(kmer_matrix) as f1:
            if headerline:
                outputfile = "t-test_results_" + phenotypes[k-1] + ".txt"
                f2 = open(outputfile, "w+")
                phenotype = phenotypes[k-1] + ": " 
            elif number_of_phenotypes > 1:
                outputfile = "t-test_results_" + str(k) + ".txt"
                f2 = open(outputfile, "w+")
                phenotype = "phenotype " + str(k) + ": "
            else:
                outputfile = "t-test_results.txt"
                f2 = open(outputfile, "w+")
                phenotype = ""
            outputfiles.append(outputfile)
            for line in f1:
                samp_w_pheno_specified = 0
                samples_x = []
                x = []
                y = []
                x_weights = []
                y_weights = []
                line=line.strip()
                kmer=line.split()[0]
                list1=line.split()[1:]
                for j in range(len(list1)):
                    if samples[samples_order[j]][k] != "NA":
                        samp_w_pheno_specified += 1
                        if list1[j] == "0":
                            y.append(float(samples[samples_order[j]][k]))
                            y_weights.append(weights[samples_order[j]])
                        else:
                            x.append(float(samples[samples_order[j]][k]))
                            x_weights.append(weights[samples_order[j]])
                            samples_x.append(samples_order[j])
                    else:
                        NA = True
                if NA == True:
                    if len(x) < 2 or len(y) < 2:
                        continue
                    elif (len(x) >= samp_w_pheno_specified - 1 
                    	    or len(y) >= samp_w_pheno_specified -1):
                        continue
                counter += 1

                #Parametes for group containig the k-mer
                wtd_mean_y = np.average(y, weights=y_weights)
                sumofweightsy = sum(y_weights)
                sumofweightsy2 = sum(i**2 for i in y_weights)
                vary = (
                	(sumofweightsy / (sumofweightsy**2 - sumofweightsy2))
                	* sum(y_weights * (y - wtd_mean_y)**2)
                	)
    
                #Parameters for group not containig the k-mer
                wtd_mean_x = np.average(x, weights=x_weights)
                sumofweightsx = sum(x_weights)
                sumofweightsx2 = sum(i**2 for i in x_weights)
                varx = (
                	(sumofweightsx / (sumofweightsx**2 - sumofweightsx2))
                	* sum(x_weights * (x - wtd_mean_x)**2)
                	)

                #Calculating the weighted Welch's t-test results
                dif = wtd_mean_x-wtd_mean_y
                sxy = math.sqrt((varx/sumofweightsx)+(vary/sumofweightsy))
                df = (
                	(((varx/sumofweightsx)+(vary/sumofweightsy))**2)
                	/ (
                		(((varx/sumofweightsx)**2) / (sumofweightsx - 1))
                		+ ((vary/sumofweightsy)**2 / (sumofweightsy - 1))
                		)
                	)
                t = dif/sxy
                pvalue = stats.t.sf(abs(t), df)*2
                pvalues.append(pvalue)

                f2.write(
                	kmer + "\t" + str(round(t, 2)) + "\t" + "%.2E" % pvalue
                	+ "\t" + str(round(wtd_mean_x, 2)) + "\t"
                	+ str(round(wtd_mean_y,2)) + "\t" + str(len(samples_x))
                	+ "\t| " + " ".join(samples_x) + "\n"
                	)
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "tests conducted.", phenotype
                    )
        pvalues_all_phenotypes.append(pvalues)
        nr_of_kmers_tested_all_phenotypes.append(counter)
        if len(phenotypes_to_analyze) > 1 and k != phenotypes_to_analyze[-1]:
            sys.stderr.write("\n")
    f1.close()
    f2.close()
    return(
    	pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes, outputfiles
    	)

def t_test(
	    kmer_matrix, samples, samples_order, number_of_phenotypes, phenotypes,
	    kmers_to_analyse, phenotypes_to_analyze=False, FDR=False, 
        headerline=False
	    ):
    # Calculates Welch t-test results for every k-mer
    sys.stderr.write("\nConducting the k-mer specific Welch t-tests:\n")
    pvalues_all_phenotypes = []
    nr_of_kmers_tested_all_phenotypes = []
    outputfiles = []
    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)
    for k in phenotypes_to_analyze:
        currentKmerNum = 1.0
        previousPercent = 0.0
        pvalues = []
        counter = 0
        NA = False
        with open(kmer_matrix) as f1:
            if headerline:
                outputfile = "t-test_results_" + phenotypes[k-1] + ".txt"
                f2 = open(outputfile, "w+")
                phenotype = phenotypes[k-1] + ": "
            elif number_of_phenotypes > 1:
                outputfile = "t-test_results_" + str(k) + ".txt"
                f2 = open(outputfile, "w+")
                phenotype = "phenotype " + str(k) + ": "
            else:
                outputfile = "t-test_results.txt"
                f2 = open(outputfile, "w+")
                phenotype = ""
            outputfiles.append(outputfile)
            for line in f1:
                samp_w_pheno_specified = 0
                samples_x = []
                x = []
                y = []
                line=line.strip()
                kmer=line.split()[0]
                list1=line.split()[1:]
                for j in range(len(list1)):
                    if samples[samples_order[j]][k] != "NA":
                        samp_w_pheno_specified += 1
                        if list1[j] == "0":
                            y.append(float(samples[samples_order[j]][k]))
                        else:
                            x.append(float(samples[samples_order[j]][k]))
                            samples_x.append(samples_order[j])
                    else:
                        NA = True 
                if NA == True:
                    if len(x) < 2 or len(y) < 2:
                        continue
                    elif (len(x) >= samp_w_pheno_specified - 1 
                    	    or len(y) >= samp_w_pheno_specified -1):
                        continue
                counter += 1
 
                #Calculating the Welch's t-test results using scipy.stats
                meanx = round((sum(x)/len(x)), 2)
                meany = round((sum(y)/len(y)), 2)
                ttest = stats.ttest_ind(x, y, equal_var=False)

                pvalues.append(ttest[1])

                f2.write(
                	kmer + "\t%.2f\t%.2E\t" % ttest + str(meanx) + "\t"
                	+ str(meany) + "\t" + str(len(samples_x))  +"\t| "
                	+ " ".join(samples_x) + "\n"
                	)
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "tests conducted.", phenotype
                    )
        pvalues_all_phenotypes.append(pvalues)
        nr_of_kmers_tested_all_phenotypes.append(counter)
        if len(phenotypes_to_analyze) > 1 and k != phenotypes_to_analyze[-1]:
            sys.stderr.write("\n")
    f1.close()
    f2.close()
    return(
    	pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes, outputfiles
    	)

def weighted_chi_squared(
	    kmer_matrix, samples, samples_order, weights, number_of_phenotypes,
	    phenotypes, phenotypes_to_analyze, kmers_to_analyse, FDR=False, 
        headerline=False
	    ):
    # Calculates Chi-squared tests for every k-mer
    sys.stderr.write(
        "\nConducting the k-mer specific weighted chi-square tests:\n"
        )
    nr_of_kmers_tested_all_phenotypes = []
    pvalues_all_phenotypes = []
    outputfiles = []
    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)
    for k in phenotypes_to_analyze:
        currentKmerNum = 1.0
        previousPercent = 0.0
        counter = 0
        pvalues = []
        with open(kmer_matrix) as f1:
            if headerline:
                outputfile = ("chi-squared_test_results_" 
                	          + phenotypes[k-1] + ".txt")
                f2 = open(outputfile, "w+")
                phenotype = phenotypes[k-1] + ": "
            elif number_of_phenotypes > 1:
                outputfile = "chi-squared_test_results_" + str(k) + ".txt"
                f2 = open(outputfile, "w+")
                phenotype = "phenotype " + str(k) + ": "
            else:
                outputfile = "chi-squared_test_results.txt"
                f2 = open(outputfile, "w+")
                phenotype = ""
            outputfiles.append(outputfile)
            for line in f1:
                samples_x = []

                line=line.strip()
                kmer=line.split()[0]
                list1=line.split()[1:]

                weights_of_res_w_kmer = 0
                weights_of_res_wo_kmer = 0
                weights_of_sens_w_kmer = 0
                weights_of_sens_wo_kmer = 0

                for j in range(len(list1)):
                    if samples[samples_order[j]][k] != "NA":
                        if (list1[j] != "0" 
                        	    and samples[samples_order[j]][k] == "1"):
                            weights_of_res_w_kmer += weights[samples_order[j]]
                            samples_x.append(samples_order[j])
                        if (list1[j] == "0" 
                        	    and samples[samples_order[j]][k] == "1"):
                            weights_of_res_wo_kmer += weights[samples_order[j]]
                        if (list1[j] != "0" 
                        	    and samples[samples_order[j]][k] == "0"):
                            weights_of_sens_w_kmer += weights[samples_order[j]]
                            samples_x.append(samples_order[j])
                        if (list1[j] == "0" 
                        	    and samples[samples_order[j]][k] == "0"):
                            weights_of_sens_wo_kmer += weights[
                                samples_order[j]
                                ]

                weights_of_res_samples = (weights_of_res_w_kmer
                	                     + weights_of_res_wo_kmer)
                weights_of_sens_samples = (weights_of_sens_w_kmer
                	                      + weights_of_sens_wo_kmer)
                weights_of_samples_w_kmer = (weights_of_res_w_kmer
                	                        + weights_of_sens_w_kmer)
                weights_of_samples_wo_kmer = (weights_of_res_wo_kmer
                	                         + weights_of_sens_wo_kmer)
                weights_of_samples_total = (weights_of_res_samples
                	                       + weights_of_sens_samples)

                weights_of_res_w_kmer_exp = (
                	(weights_of_res_samples*weights_of_samples_w_kmer) 
                	/ float(weights_of_samples_total)
                	)
                weights_of_res_wo_kmer_exp = (
                	(weights_of_res_samples*weights_of_samples_wo_kmer)
                	/ float(weights_of_samples_total)
                	)
                weights_of_sens_w_kmer_exp = (
                	(weights_of_sens_samples*weights_of_samples_w_kmer)
                	/ float(weights_of_samples_total)
                	)
                weights_of_sens_wo_kmer_exp = (
                	(weights_of_sens_samples*weights_of_samples_wo_kmer)
                	/ float(weights_of_samples_total)
                	) 

                chisquare_results = stats.chisquare(
                	[
                	weights_of_res_w_kmer, weights_of_res_wo_kmer,
                	weights_of_sens_w_kmer, weights_of_sens_wo_kmer
                	],
                	[
                	weights_of_res_w_kmer_exp, weights_of_res_wo_kmer_exp,
                	weights_of_sens_w_kmer_exp, weights_of_sens_wo_kmer_exp
                	],
                	1
                	)
                
                pvalues.append(chisquare_results[1])

                f2.write(
                	kmer + "\t%.2f\t%.2E\t" % chisquare_results 
                	+ str(len(samples_x)) +"\t| " + " ".join(samples_x) + "\n"
                	)
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "tests conducted.", phenotype
                    )
        pvalues_all_phenotypes.append(pvalues)
        nr_of_kmers_tested_all_phenotypes.append(counter)
        if len(phenotypes_to_analyze) > 1 and k != phenotypes_to_analyze[-1]:
            sys.stderr.write("\n")
    f1.close()
    f2.close()
    return(
    	pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes, outputfiles
    	)

def chi_squared(
	    kmer_matrix, samples, samples_order, number_of_phenotypes, phenotypes,
	    phenotypes_to_analyze, kmers_to_analyse, FDR=False, headerline=False
	    ):
    # Calculates Chi-squared tests for every k-mer
    sys.stderr.write("\nConducting the k-mer specific chi-square tests:\n")
    nr_of_kmers_tested_all_phenotypes = []
    pvalues_all_phenotypes = []
    outputfiles = []
    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)
    for k in phenotypes_to_analyze:
        currentKmerNum = 1.0
        previousPercent = 0.0
        counter = 0
        pvalues = []
        with open(kmer_matrix) as f1:
            if headerline:
                outputfile = ("chi-squared_test_results_" 
                	+ phenotypes[k-1] + ".txt")
                f2 = open(outputfile, "w+")
                phenotype = phenotypes[k-1] + ": "
            elif number_of_phenotypes > 1:
                outputfile = "chi-squared_test_results_" + str(k) + ".txt"
                f2 = open(outputfile, "w+")
                phenotype = "phenotype " + str(k) + ": "
            else:
                outputfile = "chi-squared_test_results.txt"
                f2 = open(outputfile, "w+")
                phenotype = ""
            outputfiles.append(outputfile)
            for line in f1:
                samples_x = []                
                counter += 1

                line=line.strip()
                kmer=line.split()[0]
                list1=line.split()[1:]

                res_w_kmer = 0
                res_wo_kmer = 0
                sens_w_kmer = 0
                sens_wo_kmer = 0

                for j in range(len(list1)):
                    if samples[samples_order[j]][k] != "NA":
                        if (list1[j] != "0" 
                        	    and samples[samples_order[j]][k] == "1"):
                            res_w_kmer += 1
                            samples_x.append(samples_order[j])
                        if (list1[j] == "0" 
                        	    and samples[samples_order[j]][k] == "1"):
                            res_wo_kmer += 1
                        if (list1[j] != "0" 
                        	    and samples[samples_order[j]][k] == "0"):
                            sens_w_kmer += 1
                            samples_x.append(samples_order[j])
                        if (list1[j] == "0" 
                        	    and samples[samples_order[j]][k] == "0"):
                            sens_wo_kmer += 1

                res_samples = (res_w_kmer + res_wo_kmer)
                sens_samples = (sens_w_kmer + sens_wo_kmer)
                samples_w_kmer = (res_w_kmer + sens_w_kmer)
                samples_wo_kmer = (res_wo_kmer + sens_wo_kmer)
                samples_total = res_samples+sens_samples

                res_w_kmer_exp = ((res_samples * samples_w_kmer)
                	             / float(samples_total))
                res_wo_kmer_exp = ((res_samples * samples_wo_kmer) 
                	              / float(samples_total))
                sens_w_kmer_exp = ((sens_samples * samples_w_kmer)
                	              / float(samples_total))
                sens_wo_kmer_exp = ((sens_samples * samples_wo_kmer)
                	               / float(samples_total))

                chisquare_results = stats.chisquare(
                	[res_w_kmer, res_wo_kmer, sens_w_kmer, sens_wo_kmer],
                	[
                	res_w_kmer_exp, res_wo_kmer_exp, 
                	sens_w_kmer_exp, sens_wo_kmer_exp
                	],
                	1
                	)
                
                pvalues.append(chisquare_results[1])

                f2.write(
                	kmer + "\t%.2f\t%.2E\t" % chisquare_results 
                	+ str(len(samples_x))  +"\t| " + " ".join(samples_x) + "\n"
                	)
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "tests conducted.", phenotype
                    )
        pvalues_all_phenotypes.append(pvalues)
        nr_of_kmers_tested_all_phenotypes.append(counter)
        if len(phenotypes_to_analyze) > 1 and k != phenotypes_to_analyze[-1]:
            sys.stderr.write("\n")
    return(
    	pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes, outputfiles)

def kmer_filtering_by_pvalue(
	    test_results, pvalue, number_of_phenotypes, phenotype_scale,
	    nr_of_kmers_tested_all_phenotypes, pvalues_all_phenotypes, phenotypes, 
        kmer_limit, kmers_to_analyse, phenotypes_to_analyze=False, FDR=False, 
        B=False, headerline=False
	    ):
    # Filters the k-mers by their p-value achieved in statistical 
    # testing.
    sys.stderr.write("\nFiltering the k-mers by p-value:\n")
    kmers_passed_all_phenotypes = []
    counter = 0
    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)
    for j, k in enumerate(phenotypes_to_analyze):
        currentKmerNum = 1.0
        previousPercent = 0.0
        kmers_passed = []
        f1 = open(test_results[j])
        if headerline:
            f2 = open(
            	"k-mers_filtered_by_pvalue_" + phenotypes[k-1] + ".txt", "w+")
            phenotype = phenotypes[k-1] + ": "
        elif number_of_phenotypes > 1:
            f2 = open("k-mers_filtered_by_pvalue_" + str(k) + ".txt", "w+")
            phenotype = "phenotype " + str(k) + ": "
        else:
            f2 = open("k-mers_filtered_by_pvalue.txt", "w+")
            phenotype = ""
        if phenotype_scale == "continuous":
            f2.write(
            	"K-mer\tWelch's_t-statistic\tp-value\t+_group_mean\
            	\t-_group_mean\tNo._of_samples_with_k-mer\
            	\tSamples_with_k-mer\n"
            	)
        else:
            f2.write(
            	"K-mer\tChi-square_statistic\tp-value\
            	\tNo._of_samples_with_k-mer\tSamples_with_k-mer\n"
            	)
        number_of_kmers = 0
        pvalues_ascending = sorted(pvalues_all_phenotypes[counter])
        max_pvalue_by_limit = float('%.2E' % pvalues_ascending[kmer_limit-1])
        if B:
            for line in f1:
                max_pvalue_by_B = (
                    pvalue/nr_of_kmers_tested_all_phenotypes[counter]
                    )
                list1 = line.split()
                if float(list1[2]) < (max_pvalue_by_B):
                    f2.write(line)
                    if (number_of_kmers < kmer_limit 
                            and float(list1[2]) <= max_pvalue_by_limit):
                        kmers_passed.append(list1[0])
                        number_of_kmers += 1
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "k-mers filtered.", phenotype
                    )
        elif FDR:
            max_pvalue_by_FDR = 0
            for i, item in enumerate(pvalues_ascending):
                if  (item  < (
                	    (float(i+1) 
                	    / nr_of_kmers_tested_all_phenotypes[counter]) * pvalue
                	    )):
                    highest_sign_pvalue = item
                elif item > pvalue:
                    break
            for line in f1:
                list1 = line.split()
                if float(list1[2]) <= max_pvalue_by_FDR:
                    f2.write(line)
                    if (number_of_kmers < kmer_limit 
                            and float(list1[2]) <= max_pvalue_by_limit):
                        kmers_passed.append(list1[0])
                        number_of_kmers += 1
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "k-mers filtered.", phenotype
                    )
        else:
            for line in f1:
                list1 = line.split()
                if  float(list1[2]) < pvalue:
                    f2.write(line)
                    if (number_of_kmers < kmer_limit 
                            and float(list1[2]) <= max_pvalue_by_limit):
                        kmers_passed.append(list1[0])
                        number_of_kmers += 1
                previousPercent, currentKmerNum = write_to_stderr_if(
                    previousPercent, currentKmerNum, 
                    kmers_to_analyse, "k-mers filtered.", phenotype
                    )
        if len(phenotypes_to_analyze) > 1 and k != phenotypes_to_analyze[-1]:
            sys.stderr.write("\n")
        if len(kmers_passed) == 0:
            f2.write("\nNo k-mers passed the filtration by p-value.\n")
        counter += 1
        f1.close()
        f2.close()
        kmers_passed_all_phenotypes.append(kmers_passed)
    return(kmers_passed_all_phenotypes)

def linear_regression(
	    kmer_matrix, samples, samples_order, alphas, number_of_phenotypes,
	    kmers_passed_all_phenotypes, penalty, n_splits, weights, testset_size,
	    phenotypes, use_of_weights, phenotypes_to_analyze=False, 
        headerline=False
	    ):
    # Applies linear regression with Lasso regularization on k-mers
    # that passed the filtering by p-value of statistical test. K-mers
    # presence/absence (0/1) in samples are used as independent
    # parameters, resistance value (continuous) is used as dependent
    # parameter.

    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)

    if len(phenotypes_to_analyze) > 1:
        sys.stderr.write("\nConducting the linear regression analysis:\n")
    elif headerline:
        sys.stderr.write("\nConducting the linear regression analysis of " 
            +  phenotypes[0] + " data...\n")
    else:
        sys.stderr.write("\nConducting the linear regression analysis...\n")

    for j, k in enumerate(phenotypes_to_analyze):
        #Open files to write results of	linear regression
        if headerline:
            f1 = open("summary_of_lin_reg_analysis" 
            	     + phenotypes[k-1] + ".txt", "w+")
            f2 = open("k-mers_and_coefficients_in_lin_reg_model_" 
            	     + phenotypes[k-1] + ".txt", "w+")
            model_filename = "lin_reg_model_" + phenotypes[k-1] + ".pkl"
            if len(phenotypes_to_analyze) > 1:
                sys.stderr.write("\tregression analysis of " 
                    +  phenotypes[k-1] + " data...\n")
        elif number_of_phenotypes > 1:
            f1 = open("summary_of_lin_reg_analysis" + str(k) + ".txt", "w+")
            f2 = open("k-mers_and_coefficients_in_lin_reg_model_" 
            	     + str(k) + ".txt", "w+")
            model_filename = "lin_reg_model_" +	phenotypes[k-1] + ".pkl"
            sys.stderr.write("\tregression analysis of phenotype " 
                +  str(k) + " data...\n")
        else:
            f1 = open("summary_of_lin_reg_analysis.txt", "w+")
       	    f2 = open("k-mers_and_coefficients_in_lin_reg_model.txt", "w+")
            model_filename = "lin_reg_model.txt"

        if len(kmers_passed_all_phenotypes[j]) == 0:
       	    f1.write("No k-mers passed the step of k-mer selection for " \
       	    	"regression analysis.\n")
       	    continue

        # Generating a binary k-mer presence/absence matrix and a list
        # of k-mer names based on information in k-mer_matrix.txt
        kmers_presence_matrix = []
        features = []
        Phenotypes = [samples[item][k] for item in samples_order]
        with open(kmer_matrix) as f3:
            for line in f3:
                if line.split()[0] in kmers_passed_all_phenotypes[j]:
                    features.append(line.split()[0])
                    kmers_presence_matrix.append(
                    	map(
                    		lambda x: 0 if x == 0 else 1, 
                    		map(int, line.split()[1:])
                    		)
                    	)
        f3.close()

        # Converting data into Python array formats suitable to use in
        # sklearn modelling. Also, deleting information associated with
        # stains missing the phenotype data
        features = np.array(features)
        Phenotypes = np.array(Phenotypes)
        kmers_presence_matrix = np.array(kmers_presence_matrix).transpose()
        samples_in_analyze = np.array(samples_order)
        to_del = []
        for i, item in enumerate(Phenotypes):
            if item == "NA":
                to_del.append(i)
        kmers_presence_matrix = np.delete(kmers_presence_matrix, to_del, 0)
        Phenotypes = map(float, np.delete(Phenotypes, to_del, 0))            
        samples_in_analyze = np.delete(samples_in_analyze, to_del, 0)

        # Insert data into linear regression dataset 
        dataset = sklearn.datasets.base.Bunch(
        	data=kmers_presence_matrix, target=Phenotypes,
        	target_names=np.array(["resistant", "sensitive"]),
        	feature_names=features
        	)  
        f1.write("Dataset:\n%s\n\n" % dataset)

        # Defining linear regression parameters    
        if penalty == 'L1':
            lin_reg = Lasso()        
        if penalty == 'L2':
            lin_reg = Ridge() 
        
        # Generate grid search classifier where parameters
        # (like regularization strength) are optimized by
        # cross-validated grid-search over a parameter grid.
        parameters = {'alpha': alphas}
        clf = GridSearchCV(lin_reg, parameters, cv=n_splits)

        # Fitting the linear regression model to dataset
        # (with or without considering the weights). Writing results
        # into corresponding files.
        if testset_size != 0.0:
            if penalty == 'L2' and use_of_weights == "+":
                array_weights = np.array(
                	[weights[item] for item in samples_in_analyze]
                	)
                (
                X_train, X_test, sample_weights_train, sample_weights_test,
                y_train, y_test, samples_train, 
                samples_test
                ) = train_test_split(
                dataset.data, array_weights, dataset.target,
                samples_in_analyze, test_size=testset_size, random_state=55
                )
                model = clf.fit(
                	X_train, y_train,
                	sample_weight=sample_weights_train
                	)
            else:
                (
                X_train, X_test, y_train, y_test, samples_train, samples_test
                ) = train_test_split(
                dataset.data, dataset.target, samples_in_analyze,
                test_size=testset_size, random_state=55
                )
                model = clf.fit(X_train, y_train)
            f1.write('Parameters:\n%s\n\n' % model)
            f1.write("Grid scores (R2 score) on development set: \n")
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score']
            for mean, std, params in zip(
            	    means, stds, clf.cv_results_['params']
            	    ):
                f1.write(
                	"%0.3f (+/-%0.03f) for %r \n" % (mean, std * 2, params)
                	)
            f1.write("\nBest parameters found on development set: \n")
            for key, value in clf.best_params_.iteritems():
                f1.write(key + " : " + str(value) + "\n")      
            f1.write("\nModel predictions on test set:\nSample_ID " \
            	    "Acutal_phenotype Predicted_phenotype\n")
            for u in range(len(samples_test)):
                f1.write('%s %s %s\n' % (
                	samples_test[u], y_test[u], clf.predict(X_test)[u]
                	))
            train_y_prediction = clf.predict(X_train)
            f1.write('\nTraining set:\n')
            f1.write('Mean squared error: %s\n' % \
            	     mean_squared_error(y_train, train_y_prediction))
            f1.write("The coefficient of determination:"
                + " %s\n" % clf.score(X_train, y_train))
            f1.write("The Spearman correlation coefficient and p-value:" \
                " %s, %s \n\n" % stats.spearmanr(y_train, train_y_prediction))
            test_y_prediction = clf.predict(X_test)
            f1.write('Test set:\n')
            f1.write('Mean squared error: %s\n' 
                % mean_squared_error(y_test, test_y_prediction))
            f1.write('The coefficient of determination:'
                + ' %s\n' % clf.score(X_test, y_test)) 
            f1.write("The Spearman correlation coefficient and p-value:" \
                + " %s, %s \n" % stats.spearmanr(y_test, test_y_prediction))
        else:
            if penalty == 'L2' and use_of_weights == "+":
                array_weights = np.array(
                	[weights[item] for item in samples_in_analyze]
                	)
                model = clf.fit(
                	dataset.data, dataset.target, sample_weight=array_weights
                	)
            else:
                model = clf.fit(dataset.data, dataset.target)
            f1.write('Parameters:\n%s\n\n' % model)
            f1.write("Grid scores (R2 score) on development set:")
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score'] 
            for mean, std, params in zip(
            	    means, stds, clf.cv_results_['params']
            	    ):
                f1.write("%0.3f (+/-%0.03f) for %r \n" % (
                	mean, std * 2, params
                	))
            f1.write("\nBest parameters found on development set: \n")
            for key, value in clf.best_params_.iteritems():
                f1.write(key + " : " + str(value) + "\n")
            y_prediction = clf.predict(dataset.data) 
            f1.write('\nMean squared error on the dataset: %s\n' % \
            	    mean_squared_error(dataset.target, y_prediction))
            f1.write("The coefficient of determination of the dataset: " \
            	    "%s\n" % clf.score(X_train, y_train))
            f1.write('The Spearman correlation coefficient and p-value of ' \
                'the dataset: %s, %s \n' % stats.spearmanr(
                    dataset.target, y_prediction
                    )
                )

        joblib.dump(model, model_filename)
        kmers_presence_matrix = np.array(kmers_presence_matrix).transpose()
        f2.write("K-mer\tcoef._in_lin_reg_model\tNo._of_samples_with_k-mer\
        	    \tSamples_with_k-mer\n")
        for x in range(len(clf.best_estimator_.coef_)):
            samples_with_kmer = [i for i,j in zip(
            	samples_in_analyze, kmers_presence_matrix[x]
            	) if j != 0]
            f2.write("%s\t%s\t%s\t| %s\n" % (
            	features[x], clf.best_estimator_.coef_[x],
            	len(samples_with_kmer), " ".join(samples_with_kmer)
            	))        
        f1.close()
        f2.close()

def logistic_regression(
	    kmer_matrix, samples, samples_order, alphas, number_of_phenotypes, 
	    kmers_passed_all_phenotypes, penalty, n_splits, weights, testset_size,
	    phenotypes, use_of_weights, phenotypes_to_analyze=False, 
	    headerline=False
	    ):
    # Applies logistic regression with Lasso regularization on k-mers
    # that passed the filtering by p-value of statistical test. K-mers
    # presence/absence (0/1) in samples are used as independent
    # parameters, resistance value (0/1) is used as dependent 
    # parameter.
    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)
    
    if len(phenotypes_to_analyze) > 1:
        sys.stderr.write("\nConducting the logistic regression analysis:\n")
    elif headerline:
        sys.stderr.write("\nConducting the logistic regression analysis of " 
            +  phenotypes[0] + " data...\n")
    else:
        sys.stderr.write("\nConducting the logistic regression analysis...\n")

    for j, k in enumerate(phenotypes_to_analyze):
        #Open files to write results of	logistic regression
        if headerline:
            f1 = open(
            	"summary_of_log_reg_analysis_" + phenotypes[k-1] + ".txt", "w+"
            	)
            f2 = open("k-mers_and_coefficients_in_log_reg_model_" 
            	     + phenotypes[k-1] + ".txt", "w+")
            model_filename = "log_reg_model_" + phenotypes[k-1] + ".pkl"
            if len(phenotypes_to_analyze) > 1:
                sys.stderr.write("\tregression analysis of " 
                    +  phenotypes[k-1] + " data...\n")
        elif number_of_phenotypes > 1:
            f1 = open("summary_of_log_reg_analysis_" + str(k) + ".txt", "w+")
            f2 = open("k-mers_and_coefficients_in_log_reg_model_" 
            	     + str(k) + ".txt", "w+")
            model_filename = "log_reg_model_" +	str(k) + ".pkl"
            sys.stderr.write("\tregression analysis of phenotype " 
                +  str(k) + " data...\n")
        else:
            f1 = open("summary_of_log_reg_analysis.txt", "w+")
       	    f2 = open("k-mers_and_coefficients_in_log_reg_model.txt", "w+")
            model_filename = "log_reg_model.pkl"
        
        if len(kmers_passed_all_phenotypes[j]) == 0:
            f1.write("No k-mers passed the step of k-mer selection for \
            	regression analysis.\n")
            continue

        # Generating a binary k-mer presence/absence matrix and a list
        # of k-mer names based on information in k-mer_matrix.txt 
        kmers_presence_matrix = []
        features = []
        Phenotypes = [samples[item][k] for item in samples_order]
        with open(kmer_matrix) as f3:
            for line in f3:
                if line.split()[0] in kmers_passed_all_phenotypes[j]:
                    features.append(line.split()[0])
                    kmers_presence_matrix.append(map(
                    	lambda x: 0 if x == 0 else 1,
                    	map(int, line.split()[1:])
                    	))
        f3.close()

        # Converting data into Python array formats suitable to use in
        # sklearn modelling. Also, deleting information associated with
        # stains missing the phenotype data
        features = np.array(features)
        Phenotypes = np.array(Phenotypes)
        kmers_presence_matrix = np.array(kmers_presence_matrix).transpose()
        samples_in_analyze = np.array(samples_order)
        to_del = []
        for i, item in enumerate(Phenotypes):
            if item == "NA":
                to_del.append(i)
        kmers_presence_matrix = np.delete(kmers_presence_matrix, to_del, 0)
        Phenotypes = map(int, np.delete(Phenotypes, to_del, 0))            
        samples_in_analyze = np.delete(samples_in_analyze, to_del, 0)

        #Insert data into logistic regression dataset  
        dataset = sklearn.datasets.base.Bunch(
        	data=kmers_presence_matrix, target=Phenotypes,
        	target_names=np.array(["resistant", "sensitive"]),
        	feature_names=features
        	) 
        f1.write("Dataset:\n%s\n\n" % dataset)

        #Defining logistic regression parameters    
        if penalty == "L1":
            log_reg = LogisticRegression(penalty='l1', solver='saga')        
        if penalty == "L2":
            log_reg = LogisticRegression(penalty='l2', solver='saga')
       
        # Generate grid search classifier where parameters
        # (like regularization strength) are optimized by
        # cross-validated grid-search over a parameter grid. 
        Cs = list(map(lambda x: 1/x, alphas))
        parameters = {'C':Cs}
        clf = GridSearchCV(log_reg, parameters, cv=n_splits)

        # Fitting the logistic regression model to dataset 
        # (with or without considering the weights). Writing logistic
        # regression results into corresponding files.
        if testset_size != 0.0:
            if use_of_weights == "+":
                array_weights = np.array(
                	[weights[item] for item in samples_in_analyze]
                	)
                (
                    X_train, X_test, sample_weights_train, sample_weights_test,
                    y_train, y_test, samples_train, samples_test
                    ) = train_test_split(
                    dataset.data, array_weights, dataset.target,
                    samples_in_analyze, test_size=testset_size,
                    stratify=dataset.target, random_state=55
                    )
                model = clf.fit(
                	X_train, y_train, sample_weight=sample_weights_train
                	)
            else:
                (
                    X_train, X_test, y_train, y_test, samples_train,
                    samples_test
                    ) = train_test_split(
                    dataset.data, dataset.target, samples_in_analyze,
                    stratify=dataset.target, test_size=testset_size, 
                    random_state=55
                    )
                model = clf.fit(X_train, y_train)
            f1.write('Parameters:\n%s\n\n' % model)
            f1.write("Grid scores (mean accuracy) on development set:\n")
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score']
            for mean, std, params in zip(
            	    means, stds, clf.cv_results_['params']
            	    ):
                f1.write("%0.3f (+/-%0.03f) for %r \n" % (
                	mean, std * 2, params
                	))
            f1.write("\nBest parameters found on development set: \n")
            for key, value in clf.best_params_.iteritems():
                f1.write(key + " : " + str(value) + "\n")
            f1.write("\n\nModel predictions on test set:\nSample_ID \
            	Acutal_phenotype Predicted_phenotype\n")
            for u in range(len(samples_test)):
                f1.write('%s %s %s\n' % (
                	samples_test[u], y_test[u], clf.predict(X_test)[u]
                	))
            f1.write("\nTraining set: \n")
            f1.write("Mean accuracy: %s\n" % \
            	    clf.score(X_train, y_train))
            f1.write("\nTest set: \n")
            f1.write('Mean accuracy: %s\n\n' % clf.score(X_test, y_test))
            f1.write('Classification report:\n %s\n\n' % classification_report(
            	y_test, clf.predict(X_test), 
            	target_names=["sensitive", "resistant"]
            	))  
        else:
            if use_of_weights == "+":
                array_weights = np.array(
                	[weights[item] for item in samples_in_analyze])
                model = clf.fit(
                	dataset.data, dataset.target, 
                	fit_params={'sample_weight': array_weights}
                	)
            else:
                model = clf.fit(dataset.data, dataset.target)
            f1.write('Parameters:\n%s\n\n' % model)
            f1.write("Grid scores on development set:\n")
            means = clf.cv_results_['mean_test_score']
            stds = clf.cv_results_['std_test_score']
            for mean, std, params in zip(
            	    means, stds, clf.cv_results_['params']
            	    ):
                f1.write("%0.3f (+/-%0.03f) for %r \n" % (
                	mean, std * 2, params
                	))
            f1.write("\nBest parameters found on development set: \n")
            for key, value in clf.best_params_.iteritems():
                f1.write(key + " : " + str(value) + "\n")
            f1.write("\nMean accuracy on the dataset: %s\n" % clf.score(
            	dataset.data, dataset.target
            	))
            f1.write('Classification report:\n %s\n\n' % classification_report(
            	dataset.target, model.predict(dataset.data), 
            	target_names=["sensitive", "resistant"]
            	))
        
        joblib.dump(model, model_filename)
        kmers_presence_matrix = np.array(kmers_presence_matrix).transpose()
        f2.write("K-mer\tcoef._in_log_reg_model\tNo._of_samples_with_k-mer\
        	    \tSamples_with_k-mer\n")
        for x in range(len(clf.best_estimator_.coef_[0])):
            samples_with_kmer = [i for i,j in zip(
            	samples_in_analyze, kmers_presence_matrix[x]
            	) if j != 0]
            f2.write("%s\t%s\t%s\t| %s\n" % (
            	features[x], clf.best_estimator_.coef_[0][x],
            	len(samples_with_kmer), " ".join(samples_with_kmer)
            	))
        f1.close()
        f2.close()

def ReverseComplement(kmer):
    # Returns the reverse complement of kmer
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return("".join([seq_dict[base] for base in reversed(kmer)]))

def string_set(string_list):
    # Removes subsequences from kmer_list
    return set(i for i in string_list
               if not any(i in s for s in string_list if i != s))

def overlap(a, b, min_length=3):
    # Returns the overlap of kmer_a and kmer_b if overlap equals or 
    # exceeds the min_length. Otherwise returns 0.
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def pick_overlaps(reads, min_olap):
    # Takes kmer_list as an input. Generates pairwise permutations of 
    # the kmers in kmer list. Finds the overlap of each pair. Returns 
    # the lists of kmers and overlap lengths of the pairs which overlap
    # by min_olap or more nucleotides.
    reada, readb, olap_lens = [], [], []
    for a, b in permutations(reads, 2):
        olap_len = overlap(a, b, min_length=min_olap)
        if olap_len > 0:
            reada.append(a)
            readb.append(b)
            olap_lens.append(olap_len)
    return reada, readb, olap_lens

def kmer_assembler(kmer_list, min_olap=None):
    # Assembles the k-mers in kmer_list which overlap by at least 
    # min_olap nucleotides.

    kmer_length = len(kmer_list[0])
    if min_olap == None:
        min_olap = kmer_length-1
    assembled_kmers = []

    # Adding the reverse-complement of each k-mer
    kmer_list = kmer_list + map(ReverseComplement, kmer_list)

    # Find the overlaping k-mers
    kmers_a, kmers_b, olap_lens = pick_overlaps(kmer_list, min_olap)

    while olap_lens != []:
        set_a = set(kmers_a)
        set_b = set(kmers_b)

        # Picking out the assembled k-mers which have no sufficient
        # overlaps anymore.
        for item in kmer_list:
            if (item not in set_a and item not in set_b
                    and ReverseComplement(item) not in assembled_kmers):
                assembled_kmers.append(item)

        # Generating new kmer_list, where overlaping elements from previous
        # kmer_list are assembled.
        kmer_list = []
        for i, olap in enumerate(olap_lens):
            kmer_list.append(kmers_a[i] + kmers_b[i][olap:])

        # Removing substrings of other elements from kmer_list.
        kmer_list = list(string_set(kmer_list))

        # Find the overlaping elements in new generated kmer_list.
        kmers_a, kmers_b, olap_lens = pick_overlaps(kmer_list, min_olap)

    for item in kmer_list:
        # Picking out the assembled k-mers to assembled_kmers set.
        if (ReverseComplement(item) not in assembled_kmers):
            assembled_kmers.append(item)
    return(assembled_kmers)

def assembling(
        kmers_passed_all_phenotypes, phenotypes, number_of_phenotypes, 
        phenotypes_to_analyze=False, headerline=False
        ):
    # Assembles the input k-mers and writes assembled sequences
    # into "assembled_kmers.txt" file in FastA format.

    if not phenotypes_to_analyze:
        phenotypes_to_analyze = range(1,number_of_phenotypes+1)

    if len(phenotypes_to_analyze) > 1:
        sys.stderr.write(
            "Assembling the k-mers used in regression model of:\n"
            )
    elif headerline:
        sys.stderr.write("Assembling the k-mers used in regression model of " 
            +  phenotypes[0] + " data...\n")
    else:
        sys.stderr.write("Assembling the k-mers used in regression model...\n")

    for j, k in enumerate(phenotypes_to_analyze):
        #Open files to write the results of k-mer assembling
        if headerline:
            f1 = open("assembled_kmers_" + phenotypes[k-1] + ".fasta", "w+")
            if len(phenotypes_to_analyze) > 1:
                sys.stderr.write("\t" + phenotypes[k-1] + "...\n")
        elif number_of_phenotypes > 1:
            f1 = open("assembled_kmers_" + str(k) + ".fasta", "w+")
            sys.stderr.write("\tphenotype " + str(k) + "...\n")
        else:
            f1 = open("assembled_kmers.fasta", "w+")
        
        kmers_to_assemble = kmers_passed_all_phenotypes[j]
        assembled_kmers = sorted(
            kmer_assembler(kmers_to_assemble), key = len
            )[::-1]
        for i, item in enumerate(assembled_kmers):
            f1.write(">seq_" + str(i+1) + "_length_" 
                + str(len(item)) + "\n" + item + "\n")
    f1.close()

def modeling(args):

    if args.alphas == None:
        alphas = np.logspace(
            math.log10(args.alpha_min),
            math.log10(args.alpha_max), num=args.n_alphas)
    else: 
        alphas = np.array(args.alphas)

    (
        samples, samples_order, n_o_s, n_o_p, phenotype_scale, headerline, phenotypes
        ) = parse_modeling_input_file(args.inputfile)

    if args.min == "0":
        args.min = 2
    if args.max == "0":
        args.max = n_o_s - 2
    
    kmer_list_generator(samples, args.length, args.cutoff, samples_order)
    dict_of_frequencies = kmer_frequencies(samples_order)
    kmers_to_analyse = kmer_filtering_by_frequency(
        dict_of_frequencies , args.min, args.max
        )
    map_samples_modeling(samples, args.length)
    vectors_to_matrix_modeling(samples_order)
    
    call(["rm -r K-mer_lists/"], shell = True)
    
    weights = []
    if args.weights == "+":   
        mash_caller(samples, args.cutoff)
        mash_output_to_distance_matrix(samples_order, "mash_distances.mat")
        dist_mat = distance_matrix_modifier("distances.mat")
        distance_matrix_to_phyloxml(samples_order, dist_mat)   
        phyloxml_to_newick("tree_xml.txt")
        weights = newick_to_GSC_weights("tree_newick.txt")
   
    if phenotype_scale == "continuous":
        if args.weights == "+":
            (
                pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes,
                test_result_files
                ) = weighted_t_test(
                "k-mer_matrix.txt", samples, samples_order, weights, n_o_p,
                phenotypes, kmers_to_analyse, args.mpheno, args.FDR, headerline
                )
        else:
            (
            pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes,
            test_result_files
            ) = t_test(
            "k-mer_matrix.txt", samples, samples_order, n_o_p, phenotypes,
            kmers_to_analyse, args.mpheno, args.FDR, headerline
            )
    elif phenotype_scale == "binary":
        if args.weights == "+":
            (
            pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes,
            test_result_files
            ) = weighted_chi_squared(
            "k-mer_matrix.txt", samples, samples_order, weights, n_o_p,
            phenotypes, args.mpheno, kmers_to_analyse, args.FDR, headerline
            )
        else:
            (
            pvalues_all_phenotypes, nr_of_kmers_tested_all_phenotypes,
            test_result_files
            ) = chi_squared(
            "k-mer_matrix.txt", samples, samples_order, n_o_p, phenotypes,
            args.mpheno, kmers_to_analyse, args.FDR, headerline
            )   

    kmers_passed_all_phenotypes = kmer_filtering_by_pvalue(
        test_result_files, args.pvalue, n_o_p, phenotype_scale, 
        nr_of_kmers_tested_all_phenotypes, pvalues_all_phenotypes, phenotypes,
        args.n_kmers, kmers_to_analyse, args.mpheno, args.FDR, args.Bonferroni, 
        headerline
        )
    
    if phenotype_scale == "continuous":
        linear_regression(
            "k-mer_matrix.txt", samples, samples_order, alphas, n_o_p,
            kmers_passed_all_phenotypes, args.regularization, args.n_splits,
            weights, args.testset_size, phenotypes, args.weights,
            args.mpheno, headerline
            )
    elif phenotype_scale == "binary":
        logistic_regression(
            "k-mer_matrix.txt", samples, samples_order, alphas, n_o_p,
            kmers_passed_all_phenotypes, args.regularization, args.n_splits,
            weights, args.testset_size, phenotypes, args.weights,
            args.mpheno, headerline
            )
    if args.assembly == "+":
        assembling(
            kmers_passed_all_phenotypes, phenotypes, n_o_p, args.mpheno, 
            headerline
            )
