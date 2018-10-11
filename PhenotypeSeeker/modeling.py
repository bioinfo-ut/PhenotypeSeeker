#!/usr/bin/python2.7

from __future__ import print_function

__author__ = "Erki Aun"
__version__ = "0.3.1"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

from itertools import chain, izip, izip_longest, permutations
from subprocess import call, Popen, PIPE, check_output
import math
import sys
import warnings
warnings.showwarning = lambda *args, **kwargs: None

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from cogent import LoadTree
from cogent.align.weights.methods import GSC
from collections import Counter, OrderedDict
from multiprocess import Manager, Pool, Value
from scipy import stats
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import (Lasso, LogisticRegression, Ridge, ElasticNet,
    SGDClassifier)
from sklearn.naive_bayes import BernoulliNB, GaussianNB
from sklearn.svm import SVC
from sklearn.metrics import (
    classification_report, r2_score, mean_squared_error, recall_score,
    roc_auc_score, average_precision_score, matthews_corrcoef, cohen_kappa_score,
    confusion_matrix, accuracy_score
    )
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, train_test_split
from functools import partial
import xgboost as xgb
import Bio
import numpy as np
import pandas as pd
import sklearn.datasets

class Input():

    samples = OrderedDict()
    phenotypes_to_analyse = OrderedDict()
    pool = None
    lock = None
    
    @classmethod
    def get_input_data(cls, inputfilename, take_logs):
        # Read the data from inputfile into "samples" directory
        samples = OrderedDict()
        Samples.take_logs = take_logs
        with open(inputfilename) as inputfile:
            for i, line in enumerate(inputfile):
                if i == 0:
                    firstline = line.split()
                    Samples.no_phenotypes = len(firstline)-2
                    if firstline[0] == "SampleID":
                        Samples.phenotypes = firstline[2:]
                        Samples.headerline = True
                        continue
                    else:
                        for j in xrange(1, Samples.no_phenotypes + 1):
                            Samples.phenotypes.append("phenotype_%s" %j)
                sample_name = line.split()[0]
                cls.samples[sample_name] = (
                    Samples.from_inputfile(line)
                    )

    # ---------------------------------------------------------
    # Set parameters for multithreading
    @classmethod
    def get_multithreading_parameters(cls):
        cls.lock = Manager().Lock()
        cls.pool = Pool(Samples.num_threads)

    # ---------------------------------------------------------
    # Functions for processing the command line input arguments

    @classmethod
    def Input_args(
            cls, alphas, alpha_min, alpha_max, n_alphas,
            gammas, gamma_min, gamma_max, n_gammas, 
            min_samples, max_samples, mpheno, kmer_length,
            cutoff, num_threads, pvalue_cutoff, kmer_limit,
            FDR, B, binary_classifier, regressor, penalty, max_iter,
            tol, l1_ratio, testset_size, kernel, n_iter,
            n_splits
            ):
        cls._get_phenotypes_to_analyse(mpheno)
        phenotypes.alphas = cls._get_alphas(
            alphas, alpha_min, alpha_max, n_alphas
            )
        phenotypes.gammas = cls._get_gammas(
            gammas, gamma_min, gamma_max, n_gammas
            )
        Samples.min_samples, Samples.max_samples = cls._get_min_max(
            min_samples, max_samples
            )
        Samples.kmer_length = kmer_length
        Samples.cutoff = cutoff
        Samples.num_threads = num_threads
        phenotypes.pvalue_cutoff = pvalue_cutoff
        phenotypes.kmer_limit = kmer_limit
        phenotypes.FDR = FDR
        phenotypes.B = B
        phenotypes.penalty = penalty.upper()
        phenotypes.max_iter = max_iter
        phenotypes.tol = tol
        phenotypes.l1_ratio = l1_ratio
        phenotypes.testset_size = testset_size
        phenotypes.kernel = kernel
        phenotypes.n_iter = n_iter
        phenotypes.n_splits = n_splits

        cls.get_model_name(regressor, binary_classifier)

    @staticmethod
    def get_model_name(regressor, binary_classifier):
        if phenotypes.scale == "continuous":
            if regressor == "lin":
                phenotypes.model_name_long = "linear regression"
                phenotypes.model_name_short = "lin_reg"
            elif regressor == "XGBR":
                phenotypes.model_name_long = "XGBRegressor"
                phenotypes.model_name_short = "XGBR"
        elif phenotypes.scale == "binary":
            if binary_classifier == "log":
                phenotypes.model_name_long = "logistic regression"
                phenotypes.model_name_short = "log_reg"
            elif binary_classifier == "SVM":
                phenotypes.model_name_long = "support vector machine"
                phenotypes.model_name_short = "SVM"
            elif binary_classifier == "RF":
                phenotypes.model_name_long = "random forest"
                phenotypes.model_name_short = "RF"
            elif binary_classifier == "NB":
                phenotypes.model_name_long = "Naive Bayes"
                phenotypes.model_name_short = "NB"
            elif binary_classifier == "XGBC":
                phenotypes.model_name_long = "XGBClassifier"
                phenotypes.model_name_short = "XGBC"
        
    @staticmethod
    def _get_alphas(alphas, alpha_min, alpha_max, n_alphas):       
        # Generating the vector of alphas (hyperparameters in regression analysis)
        # based on the given command line arguments.
        if alphas == None:
            alphas = np.logspace(
                math.log10(alpha_min),
                math.log10(alpha_max), num=n_alphas)
        else: 
            alphas = np.array(alphas)
        return alphas

    @staticmethod
    def _get_gammas(gammas, gamma_min, gamma_max, n_gammas):
        # Generating the vector of gammas 
        # (hyperparameters in SVM kernel analysis)
        # based on the given command line arguments.
        if gammas == None:
            gammas = np.logspace(
                math.log10(gamma_min),
                math.log10(gamma_max), num=n_gammas)
        else: 
            gammas = np.array(gammas)
        return gammas

    @staticmethod
    def _get_min_max(min_samples, max_samples):
        # Set the min and max arguments to default values
        min_samples = int(min_samples)
        if min_samples == 0:
            min_samples = 2
        max_samples = int(max_samples)
        if max_samples == 0:
            max_samples = Samples.no_samples - 2
        return min_samples, max_samples

    @classmethod
    def _get_phenotypes_to_analyse(cls, mpheno):
        if not mpheno:
            phenotypes_to_analyze = xrange(Samples.no_phenotypes)
        else: 
            phenotypes_to_analyze = map(lambda x: x-1, mpheno)
        for item in phenotypes_to_analyze:
            cls.phenotypes_to_analyse[Samples.phenotypes[item]] = \
                phenotypes(Samples.phenotypes[item])

class Samples():

    no_samples = 0
    no_phenoypes = 0
    phenotypes = []
    take_logs = None
    headerline = None

    kmer_length = None
    cutoff = None
    min_samples = None
    max_samples = None
    num_threads = None

    mash_distances_args = []

    def __init__(self, name, address, phenotypes, weight=1):
        self.name = name
        self.address = address
        self.phenotypes = phenotypes
        self.weight = weight
    

        Samples.no_samples += 1

    def get_kmer_lists(self):
        # Makes "K-mer_lists" directory where all lists are stored.
        # Generates k-mer lists for every sample in names_of_samples variable 
        # (list or dict).
        call(["mkdir", "-p", "K-mer_lists"])
        call(
            ["glistmaker " + self.address + " -o K-mer_lists/" 
            + self.name + " -w " + self.kmer_length + " -c " + self.cutoff], 
            shell=True
            )
        Input.lock.acquire()
        stderr_print.currentSampleNum.value += 1
        Input.lock.release()
        stderr_print.print_progress("lists generated.")

    def map_samples(self):
        # Takes k-mers, which passed frequency filtering as 
        # feature space and maps samples k-mer list to that 
        # feature space. A vector of k-mers frequency information 
        # is created for every sample.
        outputfile = "K-mer_lists/" + self.name + "_mapped.txt"
        with open(outputfile, "w+") as outputfile:
            call(
                [
                "glistquery", "K-mer_lists/" + self.name + "_" + self.kmer_length +
                ".list", "-l", "K-mer_lists/feature_vector_" + self.kmer_length +
                ".list"
                ]
                , stdout=outputfile)
        Input.lock.acquire()
        stderr_print.currentSampleNum.value += 1
        Input.lock.release()
        stderr_print.print_progress("samples mapped.")

    @classmethod
    def from_inputfile(cls, line):
        sample_phenotypes = {}
        name, address, phenotype_list = \
            line.split()[0], line.split()[1], line.split()[2:]
        if not all(x == "0" or x == "1" or x == "NA" for x in phenotype_list):
            phenotypes.scale = "continuous"
        if cls.take_logs:
            phenotype_list = map(lambda x: math.log(x, 2), phenotype_list)
        for i,j in izip(cls.phenotypes, phenotype_list):
            sample_phenotypes[i] = j
        return cls(name, address, sample_phenotypes)

    @classmethod
    def get_feature_vector(cls):
        glistmaker_args = ["glistmaker"] + \
            [sample.address for sample in Input.samples.values()] + \
            [
            '-c', cls.cutoff, '-w', cls.kmer_length, '-o', 'K-mer_lists/feature_vector'
            ]
        call(glistmaker_args)


    # -------------------------------------------------------------------
    # Functions for calculating the mash distances and GSC weights for
    # input samples.
    
    def get_mash_sketches(self):
    	mash_args = "cat " + self.address + "| mash sketch - -o K-mer_lists/" + self.name
    	process = Popen(mash_args, shell=True, stderr=PIPE)
        for line in iter(process.stderr.readline, ''):
            stderr_print(line.strip())

    @classmethod
    def get_weights(cls):
        cls.get_mash_distances()
        cls._mash_output_to_distance_matrix(Input.samples.keys(), "mash_distances.mat")
        dist_mat = cls._distance_matrix_modifier("distances.mat")
        cls._distance_matrix_to_phyloxml(Input.samples.keys(), dist_mat)   
        cls._phyloxml_to_newick("tree_xml.txt")
        stderr_print("Calculating the Gerstein Sonnhammer Coathia " \
            "weights from mash distance matrix...")
        weights = cls._newick_to_GSC_weights("tree_newick.txt")
        for key, value in weights.iteritems():
            Input.samples[key].weight = value

    @classmethod
    def get_mash_distances(cls):
        mash_args = "mash paste reference.msh K-mer_lists/*.msh"
        process = Popen(mash_args, shell=True, stderr=PIPE)
        for line in iter(process.stderr.readline, ''):
            stderr_print(line.strip())
        with open("mash_distances.mat", "w+") as f1:
            call(["mash", "dist", "reference.msh", "reference.msh"], stdout=f1)

    @classmethod
    def _mash_output_to_distance_matrix(cls, names_of_samples, mash_distances):
        with open(mash_distances) as f1:
            with open("distances.mat", "w+") as f2:
                counter = 0
                f2.write(names_of_samples[counter])
                for line in f1:
                    distance = line.split()[2]
                    f2.write("\t" + distance)
                    counter += 1
                    if counter%cls.no_samples == 0:
                        if counter != cls.no_samples**2:
                            f2.write(
                                "\n" + names_of_samples[counter/cls.no_samples]
                                )

    @staticmethod
    def _distance_matrix_modifier(distance_matrix):
        # Modifies distance matrix to be suitable argument 
        # for Bio.Phylo.TreeConstruction._DistanceMatrix function
        distancematrix = []
        with open(distance_matrix) as f1:
            counter = 2
            for line in f1:
                line = line.strip().split()
                distancematrix.append(line[1:counter])
                counter += 1
        for i in xrange(len(distancematrix)):
            for j in xrange(len(distancematrix[i])):
                distancematrix[i][j] = float(distancematrix[i][j])
        return(distancematrix)

    @staticmethod
    def _distance_matrix_to_phyloxml(samples_order, distance_matrix):
        #Converting distance matrix to phyloxml
        dm = _DistanceMatrix(samples_order, distance_matrix)
        tree_xml = DistanceTreeConstructor().nj(dm)
        with open("tree_xml.txt", "w+") as f1:
            Bio.Phylo.write(tree_xml, f1, "phyloxml")

    @staticmethod
    def _phyloxml_to_newick(phyloxml):
        #Converting phyloxml to newick
        with open("tree_newick.txt", "w+") as f1:
            Bio.Phylo.convert(phyloxml, "phyloxml", f1, "newick")

    @staticmethod
    def _newick_to_GSC_weights(newick_tree):
        # Calculating Gerstein Sonnhammer Coathia weights from Newick 
        # string. Returns dictionary where sample names are keys and GSC 
        # weights are values.
        tree=LoadTree(newick_tree)
        weights=GSC(tree)
        for item in weights:
            weights[item] = 1 - weights[item]
        return(weights)

class stderr_print():
    # --------------------------------------------------------
    # Functions and variables necessarry to show the progress 
    # information in standard error.

    currentSampleNum = Value("i", 0)
    currentKmerNum = Value("i", 0)
    previousPercent = Value("i", 0)

    def __init__(self,data):
        sys.stderr.write("\r\x1b[K"+data.__str__())
        sys.stderr.flush()

    @classmethod
    def check_progress(cls, totalKmers, text, phenotype=""):
        currentPercent = (cls.currentKmerNum.value/float(totalKmers))*100
        if int(currentPercent) > cls.previousPercent.value:
            output = "\t" + phenotype + "%d%% of " % (
                currentPercent
                ) + text
            cls.previousPercent.value = int(currentPercent)
            cls(output)

    @classmethod
    def print_progress(cls, txt):
        output = "\t%d of %d %s" % (
            cls.currentSampleNum.value, Samples.no_samples,
            txt
            )
        cls(output)

class phenotypes():

    scale = "binary"

    model_name_long = None
    model_name_short = None

    # Multithreading parameters
    vectors_as_multiple_input = []
    progress_checkpoint = Value("i", 0)
    no_kmers_to_analyse = Value("i", 0)
    test_outputfiles = dict()

    # Filtering parameters
    pvalue_cutoff = None
    kmer_limit = None
    FDR = None
    B = None

    # Machine learning parameters
    model = None
    best_model = None
    model_fitted = None
    clf = None
    penalty = None
    max_iter = None
    tol = None
    l1_ratio = None
    hyper_parameters = None
    alphas = None
    gammes = None
    testset_size = 0.0
    kernel = None
    n_iter = None
    n_splits = None
    xgb_train = None
    xgb_test = None

    # ML output file holders
    summary_file = None
    coeff_file = None
    model_file = None

    def __init__(self, name):
        self.name = name
        self.pvalues = None
        self.kmers_for_ML = set()
        self.skl_dataset = None
        self.ML_df = pd.DataFrame()
        self.ML_df_train = None
        self.ML_df_test = None
        self.X_train = None
        self.y_train = None
        self.X_test = None
        self.y_test = None
        self.train_weights = None
        self.test_weights = None
        self.X_dataset = None
        self.y_dataset = None
        self.weights_dataset = None

    # -------------------------------------------------------------------
    # Functions for calculating the association test results for kmers.
    @classmethod
    def preparations_for_kmer_testing(cls):
        if phenotypes.scale == "continuous":
            sys.stderr.write("\nConducting the k-mer specific Welch t-tests:\n")
        else:
            sys.stderr.write("\nConducting the k-mer specific chi-square tests:\n")
        cls.get_params_for_kmers_testing()

    def test_kmers_association_with_phenotype(self):
        stderr_print.currentKmerNum.value = 0
        stderr_print.previousPercent.value = 0
        pvalues_from_all_threads = Input.pool.map(
            partial(
                self.get_kmers_tested, Input.samples.values()
                ), self.vectors_as_multiple_input
            )
        self.pvalues = \
            sorted(list(chain(*pvalues_from_all_threads)))
        sys.stderr.write("\n")
        self.concatenate_test_files(self.name)

    @classmethod
    def get_params_for_kmers_testing(cls):
        cls._split_sample_vectors_for_multithreading()
        cls._splitted_vectors_to_multiple_input()
        cls.no_kmers_to_analyse.value = int(
            check_output(
                ['wc', '-l', "K-mer_lists/" + Input.samples.keys()[0] + "_mapped.txt"]
                ).split()[0]
            )
        cls.progress_checkpoint.value = int(
            math.ceil(cls.no_kmers_to_analyse.value/(100*Samples.num_threads))
            )

    @staticmethod
    def _split_sample_vectors_for_multithreading():
        for sample in Input.samples:
            call(
                [
                "split -a 5 -d -n r/" + str(Samples.num_threads) + \
                " K-mer_lists/" + sample + "_mapped.txt " + \
                "K-mer_lists/" + sample + "_mapped_"
                ],
                shell=True
                )

    @classmethod
    def _splitted_vectors_to_multiple_input(cls):
        vectors_as_multiple_input = []
        for i in xrange(Samples.num_threads):
            cls.vectors_as_multiple_input.append(
                [
                "K-mer_lists/" + sample + "_mapped_%05d" %i \
                for sample in Input.samples
                ]
                )
        

    def get_kmers_tested(self, samples, split_of_kmer_lists):
        
        pvalues = []
        counter = 0

        multithreading_code = split_of_kmer_lists[0][-5:]
        test_results_file = open(self.test_result_output(
            multithreading_code
            ), "w")
        text1_4_stderr = self.get_text1_4_stderr()
        text2_4_stderr = "tests conducted."
        for line in izip_longest(
                *[open(item) for item in split_of_kmer_lists], fillvalue = ''
            ):
            counter += 1
            if counter%self.progress_checkpoint.value == 0:
                Input.lock.acquire()
                stderr_print.currentKmerNum.value += self.progress_checkpoint.value
                Input.lock.release()
                stderr_print.check_progress(
                    self.no_kmers_to_analyse.value, text2_4_stderr, text1_4_stderr
                )
            kmer = line[0].split()[0]
            kmer_presence_vector = [j.split()[1].strip() for j in line]

            if phenotypes.scale == "binary":
                pvalue = self.conduct_chi_squared_test(
                    kmer, kmer_presence_vector,
                    test_results_file, samples
                    )
            elif phenotypes.scale == "continuous":
                pvalue = self.conduct_t_test(
                    kmer, kmer_presence_vector,
                    test_results_file, samples
                    )
            if pvalue:
                pvalues.append(pvalue)
        Input.lock.acquire()
        stderr_print.currentKmerNum.value += counter%self.progress_checkpoint.value
        Input.lock.release()
        stderr_print.check_progress(
            self.no_kmers_to_analyse.value, text2_4_stderr, text1_4_stderr
        )
        test_results_file.close()
        return(pvalues)

    def test_result_output(self, code):
        if phenotypes.scale == "continuous":
            beginning_text = "t-test_results_"
        else:
            beginning_text = "chi-squared_test_results_"

        if Samples.headerline:
            outputfile = beginning_text + \
                self.name + "_" + code + ".txt"
        elif len(Input.phenotypes_to_analyse) > 1:
            outputfile = beginning_text + \
                self.name + "_" + code + ".txt"
        else:
            outputfile = beginning_text + code + ".txt"
        return outputfile

    def get_text1_4_stderr(self):
        if Samples.headerline:
            text1_4_stderr = self.name + ": "
        elif len(Input.phenotypes_to_analyse) > 1:
            text1_4_stderr = self.name + ": "
        else:
            text1_4_stderr = ""
        return text1_4_stderr

    def conduct_t_test(
        self, kmer, kmer_presence_vector,
        test_results_file, samples
        ):
        samples_w_kmer = []
        x = []
        y = []
        x_weights = []
        y_weights = []
        
        self.get_samples_distribution_for_ttest(
            x, y, x_weights, y_weights, kmer_presence_vector,
            samples_w_kmer, samples
            )

        if len(x) < Samples.min_samples or len(y) < 2 or len(x) > Samples.max_samples:
            return None

        t_statistic, pvalue, mean_x, mean_y = self.t_test(
            x, y, x_weights, y_weights
            )

        test_results_file.write(
            kmer + "\t" + str(round(t_statistic, 2)) + "\t" + \
            "%.2E" % pvalue + "\t" + str(round(mean_x, 2)) + "\t" + \
            str(round(mean_y,2)) + "\t" + str(len(samples_w_kmer)) + "\t| " + \
            " ".join(samples_w_kmer) + "\n"
            )
        return pvalue

    def get_samples_distribution_for_ttest(
            self, x, y, x_weights, y_weights,
            kmer_presence_vector, samples_w_kmer,
            samples
            ):
        for index, sample in enumerate(samples):
            sample_phenotype = sample.phenotypes[self.name]
            if sample_phenotype != "NA":
                if kmer_presence_vector[index] == "0":
                    y.append(float(sample_phenotype))
                    y_weights.append(sample.weight)
                else:
                    x.append(float(sample_phenotype))
                    x_weights.append(sample.weight)
                    samples_w_kmer.append(sample.name)

    @staticmethod
    def t_test(x, y, x_weights, y_weights):
        #Parametes for group containig the k-mer
        wtd_mean_y = np.average(y, weights=y_weights)
        sumofweightsy = sum(y_weights)
        ybar = np.float64(sum([i*j for i,j in izip(y, y_weights)])/sumofweightsy)
        vary = sum([i*j for i,j in izip(y_weights, (y - ybar)**2)])/(sumofweightsy-1)
        
        #Parameters for group not containig the k-mer
        wtd_mean_x = np.average(x, weights=x_weights)
        sumofweightsx = sum(x_weights)
        xbar = np.float64(sum([i*j for i,j in izip(x, x_weights)])/sumofweightsx)
        varx = sum([i*j for i,j in izip(x_weights, (x - xbar)**2)])/(sumofweightsx-1)

        #Calculating the weighted Welch's t-test results
        dif = wtd_mean_x-wtd_mean_y
        sxy = math.sqrt((varx/sumofweightsx)+(vary/sumofweightsy))
        df = (((varx/sumofweightsx)+(vary/sumofweightsy))**2) / \
            ((((varx/sumofweightsx)**2)/(sumofweightsx-1)) + \
                ((vary/sumofweightsy)**2/(sumofweightsy-1)))
        t= dif/sxy
        pvalue = stats.t.sf(abs(t), df)*2

        return t, pvalue, wtd_mean_x, wtd_mean_y

    def conduct_chi_squared_test(
        self, kmer, kmer_presence, test_results_file,
        samples
        ):
        samples_w_kmer = []
        (
        w_pheno_w_kmer, w_pheno_wo_kmer, wo_pheno_w_kmer, wo_pheno_wo_kmer,
        no_samples_wo_kmer
        ) = self.get_samples_distribution_for_chisquared(
            kmer_presence, samples_w_kmer, samples
            )
        no_samples_w_kmer = len(samples_w_kmer)
        if no_samples_w_kmer < Samples.min_samples or no_samples_wo_kmer < 2 \
            or no_samples_w_kmer > Samples.max_samples:
            return None
        (w_pheno, wo_pheno, w_kmer, wo_kmer, total) = self.get_totals_in_classes(
            w_pheno_w_kmer, w_pheno_wo_kmer, wo_pheno_w_kmer, wo_pheno_wo_kmer
            )

        (
        w_pheno_w_kmer_expected, w_pheno_wo_kmer_expected,
        wo_pheno_w_kmer_expected, wo_pheno_wo_kmer_expected
        ) = self.get_expected_distribution(
            w_pheno, wo_pheno, w_kmer, wo_kmer, total)

        chisquare_results = stats.chisquare(
            [
            w_pheno_w_kmer, w_pheno_wo_kmer,
            wo_pheno_w_kmer, wo_pheno_wo_kmer
            ],
            [
            w_pheno_w_kmer_expected, w_pheno_wo_kmer_expected, 
            wo_pheno_w_kmer_expected, wo_pheno_wo_kmer_expected
            ],
            1
            )
        test_results_file.write(
            kmer + "\t%.2f\t%.2E\t" % chisquare_results 
            + str(no_samples_w_kmer)  +"\t| " + " ".join(samples_w_kmer) + "\n"
            )
        pvalue = chisquare_results[1]
        return pvalue

    def get_samples_distribution_for_chisquared(
            self, kmers_presence_vector, samples_w_kmer,
            samples
            ):
        no_samples_wo_kmer = 0
        with_pheno_with_kmer = 0
        with_pheno_without_kmer = 0
        without_pheno_with_kmer = 0
        without_pheno_without_kmer = 0
        for index, sample in enumerate(samples):
            if sample.phenotypes[self.name] == "1":
                if (kmers_presence_vector[index] != "0"):
                    with_pheno_with_kmer += sample.weight 
                    samples_w_kmer.append(sample.name)
                else:
                    with_pheno_without_kmer += sample.weight
                    no_samples_wo_kmer += 1
            elif sample.phenotypes[self.name] == "0":
                if (kmers_presence_vector[index] != "0"):
                    without_pheno_with_kmer += sample.weight
                    samples_w_kmer.append(sample.name)
                else:
                    without_pheno_without_kmer += sample.weight
                    no_samples_wo_kmer += 1
        return(
            with_pheno_with_kmer, with_pheno_without_kmer,
            without_pheno_with_kmer, without_pheno_without_kmer,
            no_samples_wo_kmer
            )

    @staticmethod
    def get_totals_in_classes(
            w_pheno_w_kmer, w_pheno_wo_kmer, wo_pheno_w_kmer, wo_pheno_wo_kmer
            ):
        w_pheno = (w_pheno_w_kmer + w_pheno_wo_kmer)
        wo_pheno = (
            wo_pheno_w_kmer + wo_pheno_wo_kmer
            )
        w_kmer = (w_pheno_w_kmer + wo_pheno_w_kmer)
        wo_kmer = (
            w_pheno_wo_kmer + wo_pheno_wo_kmer
            )
        total = w_pheno + wo_pheno
        return w_pheno, wo_pheno, w_kmer, wo_kmer, total

    @staticmethod
    def get_expected_distribution(w_pheno, wo_pheno, w_kmer, wo_kmer, total):
        w_pheno_w_kmer_expected = ((w_pheno * w_kmer)
                         / float(total))
        w_pheno_wo_kmer_expected = ((w_pheno * wo_kmer) 
                          / float(total))
        wo_pheno_w_kmer_expected  = ((wo_pheno * w_kmer)
                          / float(total))
        wo_pheno_wo_kmer_expected = ((wo_pheno * wo_kmer)
                           / float(total))
        return(
            w_pheno_w_kmer_expected, w_pheno_wo_kmer_expected,
            wo_pheno_w_kmer_expected, wo_pheno_wo_kmer_expected
            )

    def concatenate_test_files(self, phenotype):
        if phenotypes.scale == "continuous":
            beginning_text = "t-test_results_"
        else:
            beginning_text = "chi-squared_test_results_"
        if Samples.no_phenotypes > 1:
            outputfile = beginning_text + phenotype + ".txt"
        else:
            outputfile = beginning_text[:-1] + ".txt"
        self.test_outputfiles[phenotype] = outputfile
        if Samples.headerline:
            call(
                [
                "cat " + beginning_text + phenotype + "_* > " +
                outputfile
                ],
                shell=True
                )
            for l in xrange(Samples.num_threads):
                call(
                    [
                    "rm " + beginning_text + phenotype +
                    "_%05d.txt" % l
                    ],
                    shell=True
                    )
        elif Samples.no_phenotypes > 1:
            call(
                [
                "cat " + beginning_text + phenotype + "_* > " +
                outputfile
                ],
                shell=True
                )
            for l in xrange(Samples.num_threads):
                call(
                    ["rm " + beginning_text + phenotype + "_%05d.txt" %l],
                    shell=True
                    )     
        else:
            call(
                [
                "cat " + beginning_text + "* > " + outputfile +
                " && rm " + beginning_text + "*"
                ],
                shell=True
                )

    # -------------------------------------------------------------------
    # Functions for filtering the k-mers based on the p-values of
    # conducted tests.
    def get_kmers_filtered(self):
        # Filters the k-mers by their p-value achieved in statistical 
        # testing.
        pvalues = self.pvalues
        phenotype = self.name
        nr_of_kmers_tested = float(len(pvalues))
        self.get_pvalue_cutoff(pvalues, nr_of_kmers_tested)
        max_pvalue_by_limit = float('%.2E' % pvalues[self.kmer_limit-1])

        stderr_print.currentKmerNum.value = 0
        stderr_print.previousPercent.value = 0
        text1_4_stderr = self.get_text1_4_stderr()
        text2_4_stderr = "k-mers filtered."
        checkpoint = int(math.ceil(nr_of_kmers_tested/100))
        inputfile = open(self.test_outputfiles[phenotype])
        outputfile = open(self.kmers_filtered_output(phenotype), "w")
        self.write_headerline(outputfile)

        counter = 0
        for line in inputfile:
            counter += 1
            line_to_list = line.split()
            if float(line_to_list[2]) < self.pvalue_cutoff:
                outputfile.write(line)
                if float(line_to_list[2]) <= max_pvalue_by_limit:
                        self.kmers_for_ML.add(line_to_list[0])
            if counter%checkpoint == 0:
                stderr_print.currentKmerNum.value += checkpoint
                stderr_print.check_progress(
                    nr_of_kmers_tested, text2_4_stderr, text1_4_stderr
                )

        stderr_print.currentKmerNum.value += counter%checkpoint
        stderr_print.check_progress(
            nr_of_kmers_tested, text2_4_stderr, text1_4_stderr
            )
        sys.stderr.write("\n")
        if len(self.kmers_for_ML) == 0:
            outputfile.write("\nNo k-mers passed the filtration by p-value.\n")
        inputfile.close()
        outputfile.close()

    def get_pvalue_cutoff(self, pvalues, nr_of_kmers_tested):
        if self.B:
            self.pvalue_cutoff = (cls.pvalue_cutoff/nr_of_kmers_tested)
        elif self.FDR:
            pvalue_cutoff_by_FDR = 0
            for index, pvalue in enumerate(pvalues):
                if  (pvalue  < (
                        (index+1) 
                        / nr_of_kmers_tested) * self.pvalue_cutoff
                        ):
                    pvalue_cutoff_by_FDR = pvalue
                elif item > pvalue:
                    break
            self.pvalue_cutoff = pvalue_cutoff_by_FDR

    def kmers_filtered_output(self, phenotype):
        if Samples.headerline:
            outputfile = "k-mers_filtered_by_pvalue_" + self.name + ".txt"
        elif Samples.no_phenotypes > 1:
            outputfile = "k-mers_filtered_by_pvalue_" + self.name + ".txt"
        else:
            outputfile = "k-mers_filtered_by_pvalue.txt"
        return outputfile

    @staticmethod
    def write_headerline(outputfile):
        if phenotypes.scale == "continuous":
            outputfile.write(
                "K-mer\tWelch's_t-statistic\tp-value\t+_group_mean\
                \t-_group_mean\tNo._of_samples_with_k-mer\
                \tSamples_with_k-mer\n"
                )
        elif phenotypes.scale == "binary":
            outputfile.write(
                "K-mer\tChi-square_statistic\tp-value\
                \tNo._of_samples_with_k-mer\tSamples_with_k-mer\n"
                )

    # -------------------------------------------------------------------
    # Functions for generating the machine learning model.  
    @classmethod
    def preparations_for_modeling(cls):
        if len(Input.phenotypes_to_analyse) > 1:
            sys.stderr.write("Generating the " + cls.model_name_long + " model:\n")
        elif Samples.headerline:
            sys.stderr.write("Generating the " + cls.model_name_long + " model of " 
                +  Samples.phenotypes[0] + " data...\n")
        else:
            sys.stderr.write("Generating the " + cls.model_name_long + " model...\n")
        cls.set_model()
        cls.set_hyperparameters()
        cls.get_best_model()


    @classmethod
    def set_model(cls):
        if cls.scale == "continuous":
            if cls.model_name_short == "lin_reg":
                # Defining linear regression parameters    
                if cls.penalty == 'L1':
                    cls.model = Lasso(max_iter=cls.max_iter, tol=cls.tol)        
                if cls.penalty == 'L2':
                    cls.model = Ridge(max_iter=cls.max_iter, tol=cls.tol)
                if cls.penalty == 'elasticnet' or "L1+L2":
                    cls.model = ElasticNet(
                        l1_ratio=cls.l1_ratio, max_iter=cls.max_iter, tol=cls.tol
                        )
            elif cls.model_name_short == "XGBR":
                cls.model = xgb.XGBRegressor()
        elif cls.scale == "binary":
            if cls.model_name_long == "logistic regression":
                #Defining logistic regression parameters
                if cls.penalty == "L1":
                    cls.model = LogisticRegression(
                        penalty='l1', solver='liblinear',
                        max_iter=cls.max_iter, tol=cls.tol
                        )        
                elif cls.penalty == "L2":
                    cls.model = LogisticRegression(
                        penalty='l2', solver='saga',
                        max_iter=cls.max_iter, tol=cls.tol
                        )
                elif cls.penalty == "elasticnet" or "L1+L2":
                    cls.model = SGDClassifier(
                        penalty='elasticnet', l1_ratio=cls.l1_ratio,
                        max_iter=cls.max_iter, tol=cls.tol, loss='log'
                        )
            elif cls.model_name_long == "support vector machine":
                cls.model = SVC(
                    kernel=cls.kernel, probability=True,
                    max_iter=cls.max_iter, tol=cls.tol
                    ) 
            elif cls.model_name_long == "random forest":
                cls.model = RandomForestClassifier(n_estimators=100)
            elif cls.model_name_long == "Naive Bayes":
                cls.model = BernoulliNB()
            elif cls.model_name_short == "XGBC":
                cls.model = xgb.XGBClassifier()

    @classmethod
    def set_hyperparameters(cls):
        if cls.scale == "continuous":
            if cls.model_name_short == "lin_reg":
                # Defining linear regression parameters    
                cls.hyper_parameters = {'alpha': cls.alphas}
        elif cls.scale == "binary":
            if cls.model_name_long == "logistic regression":
                #Defining logistic regression parameters
                if cls.penalty == "L1" or "L2":
                    Cs = list(map(lambda x: 1/x, cls.alphas))
                    cls.hyper_parameters = {'C':Cs}
                elif penalty == "elasticnet":
                    cls.hyper_parameters = {'alpha': cls.alphas}
            elif cls.model_name_long == "support vector machine":
                Cs = list(map(lambda x: 1/x, cls.alphas))
                Gammas = list(map(lambda x: 1/x, cls.gammas))
                if cls.kernel == "linear":
                    cls.hyper_parameters = {'C':Cs}
                if cls.kernel == "rbf":
                    cls.hyper_parameters = {'C':Cs, 'gamma':Gammas}

    @classmethod
    def get_best_model(cls):
        if cls.scale == "continuous":
            if cls.model_name_short == "lin_reg":
                cls.best_model = GridSearchCV(
                    cls.model, cls.hyper_parameters, cv=cls.n_splits
                    )
            elif cls.model_name_short == "XGBR":
                cls.best_model = cls.model
        elif cls.scale == "binary":
            if cls.model_name_long == "logistic regression":
                cls.best_model = GridSearchCV(
                    cls.model, cls.hyper_parameters, cv=cls.n_splits
                    )
            elif cls.model_name_long == "support vector machine":
                if cls.kernel == "linear":
                    cls.best_model = GridSearchCV(
                        cls.model, cls.hyper_parameters, cv=cls.n_splits
                        )
                if cls.kernel == "rbf":
                    cls.best_model = RandomizedSearchCV(
                        cls.model, cls.hyper_parameters,
                        n_iter=cls.n_iter, cv=cls.n_splits
                        )
            elif cls.model_name_short in ("RF", "NB", "XGBC"):
                cls.best_model = cls.model

    def machine_learning_modelling(self):
        if len(Input.phenotypes_to_analyse) > 1:
            sys.stderr.write("\tof " 
                +  self.name + " data...\n")
        self.get_outputfile_names()
        if len(self.kmers_for_ML) == 0:
            self.summary_file.write("No k-mers passed the step of k-mer filtering for " \
                "machine learning modelling.\n")
            return
        self.get_dataframe_for_machine_learning()
        self.fit_model()
        self.cross_validation_results()

        self.summary_file.write('\nTraining set:\n')
        self.predict(self.X_train, self.y_train)
        if self.testset_size != 0.0:
            self.summary_file.write('\nTest set:\n')
            self.predict(self.X_test, self.y_test)

        joblib.dump(self.model, self.model_file)
        self.write_model_coefficients_to_file()

        self.summary_file.close()
        self.coeff_file.close()
        self.model_file.close()


    def get_outputfile_names(self):
        if Samples.headerline:
            summary_file = "summary_of_" + self.model_name_short + "_analysis_" \
                + self.name + ".txt"
            coeff_file = "k-mers_and_coefficients_in_" + self.model_name_short \
                + "_model_" + self.name + ".txt"
            model_file = self.model_name_short + "_model_" + self.name + ".pkl"
        elif len(Input.phenotypes_to_analyse) > 1:
            summary_file = "summary_of_" + self.model_name_short + "_analysis_" \
                + self.name + ".txt"
            coeff_file = "k-mers_and_coefficients_in_" + self.model_name_short \
                + "_model_" + self.name + ".txt"
            model_file = self.model_name_short +"_model_" + self.name + ".pkl"
        else:
            summary_file = "summary_of_" + self.model_name_short + "_analysis.txt"
            coeff_file = "k-mers_and_coefficients_in_" + self.model_name_short \
                + "_model.txt"
            model_file = self.model_name_short + "_model.txt"
        
        self.summary_file = open(summary_file, "w")
        self.coeff_file = open(coeff_file, "w")
        self.model_file = open(model_file, "w")



    def get_dataframe_for_machine_learning(self):
        kmer_lists = ["K-mer_lists/" + sample + "_mapped.txt" for sample in Input.samples]
        for line in izip_longest(*[open(item) for item in kmer_lists], fillvalue = ''):
            if line[0].split()[0] in self.kmers_for_ML:
                self.ML_df[line[0].split()[0]] = [int(j.split()[1].strip()) for j in line]
        self.ML_df = self.ML_df.astype(bool).astype(int)
        self.ML_df['phenotype'] = [
            sample.phenotypes[self.name] for sample in Input.samples.values()
            ]
        self.ML_df['weight'] = [
            sample.weight for sample in Input.samples.values()
            ]
        self.ML_df.index = Input.samples.keys()
        self.ML_df = self.ML_df.loc[self.ML_df.phenotype != 'NA']
        self.ML_df = self.ML_df.T.drop_duplicates().T
        self.skl_dataset = sklearn.datasets.base.Bunch(
            data=self.ML_df.iloc[:,0:-2].values, target=self.ML_df['phenotype'].values,
            target_names=np.array(["resistant", "sensitive"]),
            feature_names=self.ML_df.iloc[:,0:-2].columns.values
            )
        if self.testset_size != 0.0:
            if phenotypes.scale == "continuous":
                stratify = None
            elif phenotypes.scale == "binary":
                stratify = self.skl_dataset.target
            self.ML_df_train, self.ML_df_test = train_test_split(
                self.ML_df, test_size=self.testset_size,
                stratify=stratify, random_state=55
                )
            self.X_train = self.ML_df_train.iloc[:,0:-2]
            self.y_train = self.ML_df_train.iloc[:,-2:-1]
            self.weights_train = self.ML_df_train.iloc[:,-1:]
            self.X_test = self.ML_df_test.iloc[:,0:-2]
            self.y_test = self.ML_df_test.iloc[:,-2:-1]
            self.weights_test = self.ML_df_test.iloc[:,-1:]
        else:
            self.X_train = self.ML_df.iloc[:,0:-2]
            self.y_train = self.ML_df.iloc[:,-2:-1]
            self.weights_train = self.ML_df.iloc[:,-1:]

        if phenotypes.scale == "continuous":
            self.y_train = self.y_train.astype(float)
            if self.testset_size != 0.0:
                self.y_test = self.y_test.astype(float)    
        elif phenotypes.scale == "binary":
            self.y_train = self.y_train.astype(int)
            if self.testset_size != 0.0:
                self.y_test = self.y_test.astype(int)

        self.summary_file.write("Dataset:\n%s\n\n" % self.skl_dataset)  


    def fit_model(self):
        if self.scale == "continuous":
            if self.model_name_short == "lin_reg":
                if self.penalty in ("L1", "elasticnet"):
                    self.model_fitted = self.best_model.fit(self.X_train, self.y_train.values.flatten())
                elif self.penalty == L2:
                    self.model_fitted = self.best_model.fit(
                        self.X_train, self.y_train.values.flatten(),
                        sample_weight=self.weights_train.values.flatten()
                        )
            elif self.model_name_short == "XGBR":
                self.model_fitted = self.best_model.fit(self.X_train.values, self.y_train.values.flatten())
        elif self.scale == "binary":
            if self.model_name_short == "XGBC":
                self.model_fitted = self.best_model.fit(self.X_train.values, self.y_train.values.flatten())
            else:
                self.model_fitted = self.best_model.fit(
                    self.X_train, self.y_train.values.flatten(),
                    sample_weight=self.weights_train.values.flatten()
                    )


    def cross_validation_results(self):
        if self.model_name_short not in ("RF", "NB", "XGBC", "XGBR"):
            self.summary_file.write('Parameters:\n%s\n\n' % self.model)
            if self.phenotype == "continuous":
                self.summary_file.write("Grid scores (R2 score) on development set: \n")
            elif self.phenotype == "binary":
                self.summary_file.write("Grid scores (mean accuracy) on development set: \n")
            means = self.best_model.cv_results_['mean_test_score']
            stds = self.best_model.cv_results_['std_test_score']
            params = self.best_model.cv_results_['params']
            for mean, std, param in izip(
                    means, stds, params
                    ):
                self.summary_file.write(
                    "%0.3f (+/-%0.03f) for %r \n" % (mean, std * 2, params)
                    )
            self.summary_file.write("\nBest parameters found on development set: \n")
            for key, value in self.best_model.best_params_.iteritems():
                self.summary_file.write(key + " : " + str(value) + "\n")

    def predict(self, dataset, labels):
        predictions = self.best_model.predict(dataset.values)
        self.summary_file.write("\nModel predictions on samples:\nSample_ID " \
            "Acutal_phenotype Predicted_phenotype\n")
        for index, row in dataset.iterrows():
                self.summary_file.write('%s %s %s\n' % (
                    index, labels.loc[index].values[0],
                    self.best_model.predict(row.reshape(1, -1))[0]
                    ))
        if self.scale == "continuous":
            self.model_performance_regressor(dataset, labels.values.flatten(), predictions)
        elif self.scale == "binary":
            self.model_performance_classifier(dataset, labels.values.flatten(), predictions)

    def model_performance_regressor(self, dataset, labels, predictions):
        self.summary_file.write('\nMean squared error: %s\n' % \
                 mean_squared_error(labels, predictions))
        self.summary_file.write("The coefficient of determination:"
            + " %s\n" % self.best_model.score(dataset.values, labels))
        self.summary_file.write("The Spearman correlation coefficient and p-value:" \
            " %s, %s \n" % stats.spearmanr(labels, predictions))
        r_value, pval_r = \
            stats.pearsonr(labels, predictions)
        self.summary_file.write("The Pearson correlation coefficient and p-value: " \
                " %s, %s \n" % (r_value, pval_r))
        self.summary_file.write("The plus/minus 1 dilution factor accuracy (for MICs):" \
            " %s \n\n" % self.within_1_tier_accuracy(
                labels, predictions
                )
            )

    def model_performance_classifier(self, dataset, labels, predictions):
            self.summary_file.write("Mean accuracy: %s\n" % self.best_model.score(dataset, labels))
            self.summary_file.write("Sensitivity: %s\n" % \
                    recall_score(labels, predictions))
            self.summary_file.write("Specificity: %s\n" % \
                    recall_score(
                        list(map(lambda x: 1 if x == 0 else 0, labels)), 
                        list(map(lambda x: 1 if x == 0 else 0, predictions
                        ))))
            self.summary_file.write("AUC-ROC: %s\n" % \
                roc_auc_score(labels, predictions, average="micro"))
            self.summary_file.write("Average precision: %s\n" % \
                average_precision_score(
                    labels, 
                    self.best_model.predict_proba(dataset)[:,1]
                    )                        )
            self.summary_file.write("MCC: %s\n" % \
                matthews_corrcoef(labels, predictions))
            self.summary_file.write("Cohen kappa: %s\n" %\
                cohen_kappa_score(labels, predictions))
            self.summary_file.write("Very major error rate: %s\n" %\
                self.VME(labels, predictions))
            self.summary_file.write("Major error rate: %s\n" %\
                self.ME(labels, predictions))
            self.summary_file.write('Classification report:\n\n %s\n' % classification_report(
                labels, predictions, 
                target_names=["sensitive", "resistant"]
                ))
            cm = confusion_matrix(labels, predictions)
            self.summary_file.write("Confusion matrix:\n")
            self.summary_file.write("Predicted\t0\t1:\n")
            self.summary_file.write("Actual\n")
            self.summary_file.write("0\t\t%s\t%s\n" % tuple(cm[0]))
            self.summary_file.write("1\t\t%s\t%s\n\n" % tuple(cm[1]))


    def write_model_coefficients_to_file(self):
        self.coeff_file.write("K-mer\tcoef._in_" + self.model_name_short + \
            "_model\tNo._of_samples_with_k-mer\tSamples_with_k-mer\n")
        df_for_coeffs = self.ML_df.iloc[:,0:-2]
        if self.model_name_short == "lin_reg":
            df_for_coeffs.loc['coefficient'] = \
                self.best_model.best_estimator_.coef_
        elif self.model_name_short in ("RF"):
            df_for_coeffs.loc['coefficient'] = \
                self.best_model.feature_importances_
        elif self.model_name_short in ("XGBR", "XGBC"):
            df_for_coeffs.loc['coefficient'] = \
                self.best_model.feature_importances_
        elif self.model_name_short in ("SVM", "log_reg"):
            if self.kernel != "rbf":
                df_for_coeffs.loc['coefficient'] = \
                    self.best_model.best_estimator_.coef_[0]
        for kmer in df_for_coeffs:
            if self.kernel == "rbf" or self.model_name_short == "NB":
                kmer_coef = "NA"
            else:
                kmer_coef = df_for_coeffs[kmer].loc['coefficient']
            samples_with_kmer = \
                df_for_coeffs.loc[df_for_coeffs[kmer] == 1].index.tolist()
            self.coeff_file.write("%s\t%s\t%s\t| %s\n" % (
                kmer, kmer_coef,
                len(samples_with_kmer), " ".join(samples_with_kmer)
                ))

    # ---------------------------------------------------------
    # Self-implemented performance measure functions
    @staticmethod
    def VME(targets, predictions):
        # Function to calculate the very major error (VME) rate
        VMEs = 0
        for item in izip(targets, predictions):
            if item[0] == 1 and item[1] == 0:
                VMEs += 1
        VME = str(float(VMEs)/len(targets)*100)+"%"
        return VME

    @staticmethod
    def ME(targets, predictions):
        # Function to calculate the major error (ME) rate
        MEs = 0
        for item in izip(targets, predictions):
            if item[0] == 0 and item[1] == 1:
                 MEs += 1
        ME = str(float(MEs)/len(targets)*100)+"%"
        return ME

    @staticmethod
    def within_1_tier_accuracy(targets, predictions):
        # Calculate the plus/minus one dilution factor accuracy
        # for predicted antibiotic resistance values.
        within_1_tier = 0
        for item in izip(targets, predictions):
            if abs(item[0]-item[1]) <= 1:
                within_1_tier +=1
        accuracy = float(within_1_tier)/len(targets)
        return accuracy

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
    # Takes kmer_list as an Input. Generates pairwise permutations of 
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

def assembling(kmers_passed_all_phenotypes, phenotypes_to_analyze):
    # Assembles the input k-mers and writes assembled sequences
    # into "assembled_kmers.txt" file in FastA format.

    if len(phenotypes_to_analyze) > 1:
        sys.stderr.write(
            "Assembling the k-mers used in regression modeling of:\n"
            )
    elif Samples.headerline:
        sys.stderr.write("Assembling the k-mers used in modeling of " 
            +  Samples.phenotypes[0] + " data...\n")
    else:
        sys.stderr.write(
            "Assembling the k-mers used in modeling...\n"
            )

    for j, k in enumerate(phenotypes_to_analyze):
        phenotype = Samples.phenotypes[k]
        if len(kmers_passed_all_phenotypes[j]) == 0:
            f1.write("No k-mers passed the step of k-mer selection for \
                assembling.\n")
            continue
        #Open files to write the results of k-mer assembling
        if Samples.headerline:
            f1 = open("assembled_kmers_" + phenotype + ".fasta", "w+")
            if len(phenotypes_to_analyze) > 1:
                sys.stderr.write("\t" + phenotype + "...\n")
        elif Samples.no_phenotypes > 1:
            f1 = open("assembled_kmers_" + phenotype + ".fasta", "w+")
            sys.stderr.write("\tphenotype " + phenotype + "...\n")
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
    # The main function of "phenotypeseeker modeling"

    # Processing the input data
    Input.get_input_data(args.inputfile, args.take_logs)
    Input.Input_args(
        args.alphas, args.alpha_min, args.alpha_max, args.n_alphas,
        args.gammas, args.gamma_min, args.gamma_max, args.n_gammas,
        args.min, args.max, args.mpheno, args.length, args.cutoff,
        args.num_threads, args.pvalue, args.n_kmers, args.FDR, 
        args.Bonferroni, args.binary_classifier, args.regressor, 
        args.penalty, args.max_iter, args.tol, args.l1_ratio,
        args.testset_size, args.kernel, args.n_iter, args.n_splits
        )
    Input.get_multithreading_parameters()

    # Operations with samples
    sys.stderr.write("Generating the k-mer lists for input samples:\n")
    Input.pool.map(
        lambda x: x.get_kmer_lists(), Input.samples.values()
        )
    sys.stderr.write("\nGenerating the k-mer feature vector.\n")
    Samples.get_feature_vector()
    sys.stderr.write("Mapping samples to the feature vector space:\n")
    stderr_print.currentSampleNum.value = 0
    Input.pool.map(
        lambda x: x.map_samples(), Input.samples.values()
        )
    if args.weights == "+":
        sys.stderr.write("\nEstimating the Mash distances between samples...\n")
    	Input.pool.map(
	        lambda x: x.get_mash_sketches(), Input.samples.values()
	        )
        Samples.get_weights()

    # Analyses of phenotypes
    phenotypes.preparations_for_kmer_testing()
    map(
        lambda x:  x.test_kmers_association_with_phenotype(), 
        Input.phenotypes_to_analyse.values()
        )
    sys.stderr.write("Filtering the k-mers by p-value:\n")
    map(
        lambda x:  x.get_kmers_filtered(), 
        Input.phenotypes_to_analyse.values()
        )
    phenotypes.preparations_for_modeling()
    map(
        lambda x: x.machine_learning_modelling(),
        Input.phenotypes_to_analyse.values()
        )

    if args.assembly == "+":
        assembling(kmers_passed_all_phenotypes, args.mpheno)