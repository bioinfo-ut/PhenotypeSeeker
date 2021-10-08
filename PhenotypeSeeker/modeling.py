#!/usr/bin/env python3

__author__ = "Erki Aun"
__version__ = "0.7.3"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

from itertools import chain, permutations, groupby
from subprocess import call, Popen, PIPE, check_output, run
import math
import os
import sys
import warnings
import pkg_resources
import joblib
import matplotlib
import matplotlib.pyplot as plt
import xgboost as xgb
import Bio
import numpy as np
import pandas as pd

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from collections import OrderedDict
from ete3 import Tree
from multiprocess import Manager, Pool, Value
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier, plot_tree
from sklearn.linear_model import (Lasso, LogisticRegression, Ridge, ElasticNet,
    SGDClassifier)
from sklearn.naive_bayes import BernoulliNB, GaussianNB
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import (
    classification_report, r2_score, mean_squared_error, recall_score,
    roc_auc_score, average_precision_score, matthews_corrcoef, cohen_kappa_score,
    confusion_matrix, accuracy_score, f1_score
    )
from sklearn.model_selection import (
    RandomizedSearchCV, GridSearchCV, train_test_split, StratifiedKFold,
    KFold
    )
from functools import partial

matplotlib.use('agg')
warnings.showwarning = lambda *args, **kwargs: None
pkg_resources.require(
    "numpy==1.18.1", "Biopython==1.76", "pandas==1.0.1", "xgboost==1.0.1", "scipy==1.4.1",
    "scikit-learn==0.22.1", "ete3==3.1.1", "multiprocess==0.70.9"
    )

import time
def timer(f):
    def wrapper(*args):
        start = time.time()
        f(*args)
        with open("log.txt", "a") as log:
            log.write(f"Func {f} took {time.time() - start} secs\n")
    return wrapper

class Input():

    samples = OrderedDict()
    phenotypes_to_analyse = OrderedDict()
    mpheno_to_index = []
    lock = Manager().Lock()

    jump_to = None
    num_threads = 8
    
    @classmethod
    def get_input_data(cls, inputfilename, take_logs, mpheno):
        # Read the data from inputfile into "samples" directory
        Samples.take_logs = take_logs
        with open(inputfilename) as inputfile:
            header = inputfile.readline().split()
            Samples.phenotypes = header[2:]
            Samples.no_phenotypes = len(header)-2
            for pheno in Samples.phenotypes:
                try:
                    float(pheno)
                    sys.stderr.write("\x1b[1;33mWarning! It seems that the input file " \
                        "is missing header row!\x1b[0m\n")
                    sys.stderr.flush()
                    break
                except ValueError:
                    pass
            for line in inputfile:
                if line.strip():
                    sample_name = line.split()[0]
                    cls.samples[sample_name] = (
                        Samples.from_inputfile(line)
                        )
        cls._get_phenotypes_to_analyse(mpheno)
        cls._set_phenotype_values(take_logs)

    @classmethod
    def _get_phenotypes_to_analyse(cls, mpheno):
        if not mpheno:
            cls.mpheno_to_index = range(Samples.no_phenotypes)
        else: 
            cls.mpheno_to_index = map(lambda x: x-1, mpheno)
        for item in cls.mpheno_to_index:
            cls.phenotypes_to_analyse[Samples.phenotypes[item]] = \
                phenotypes(Samples.phenotypes[item])

    # # ---------------------------------------------------------
    # # Set parameters for multithreading
    # @classmethod
    # def _get_multithreading_parameters(cls):
    #     cls.lock = Manager().Lock()

    @classmethod
    def _set_phenotype_values(cls, take_logs):
        for sample in cls.samples.values():
            for phenotype in cls.phenotypes_to_analyse.values():
                if phenotypes.pred_scale == "continuous":
                    try:
                        sample.phenotypes[phenotype.name] = float(sample.phenotypes[phenotype.name])
                        if take_logs:
                            sample.phenotypes[phenotype.name] = math.log(
                                sample.phenotypes[phenotype.name], 2
                                )
                    except:
                        phenotype.no_samples -= 1
                elif phenotypes.pred_scale == "binary":
                    try:
                        sample.phenotypes[phenotype.name] = int(sample.phenotypes[phenotype.name])
                    except:
                        phenotype.no_samples -= 1

    @classmethod
    def pop_phenos_out_of_kmers(cls):
        [Input.phenotypes_to_analyse.pop(pt) for pt in phenotypes.no_results]
        if len(Input.phenotypes_to_analyse) == 0:
            sys.stderr.write("\x1b[1;33mThere are no k-mers left for modelling for any phenotype.\x1b[0m\n")
            sys.stderr.write("\x1b[1;33mExiting PhenotypeSeeker\x1b[0m\n")
            sys.stderr.write("\n\x1b[1;1;101m######          PhenotypeSeeker modeling finished          ######\x1b[0m\n")
            raise SystemExit()

    # ---------------------------------------------------------
    # Functions for processing the command line input arguments

    @classmethod
    def Input_args(
            cls, alphas, alpha_min, alpha_max, n_alphas,
            gammas, gamma_min, gamma_max, n_gammas, 
            min_samples, max_samples, kmer_length,
            cutoff, num_threads, pvalue_cutoff, kmer_limit,
            FDR, B, binary_classifier, regressor, penalty, max_iter,
            tol, l1_ratio, n_splits_cv_outer, kernel, n_iter,
            n_splits_cv_inner, testset_size, train_on_whole,
            logreg_solver, jump_to, pca
            ):
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
        Input.num_threads = num_threads
        Input.jump_to = jump_to
        phenotypes.pvalue_cutoff = pvalue_cutoff
        phenotypes.kmer_limit = kmer_limit
        phenotypes.FDR = FDR
        phenotypes.B = B
        phenotypes.penalty = penalty.upper()
        phenotypes.max_iter = max_iter
        phenotypes.tol = tol
        phenotypes.l1_ratio = l1_ratio      
        phenotypes.kernel = kernel
        phenotypes.n_iter = n_iter
        phenotypes.testset_size = testset_size
        phenotypes.train_on_whole = train_on_whole
        cls.get_model_name(regressor, binary_classifier)
        phenotypes.n_splits_cv_outer = n_splits_cv_outer
        phenotypes.n_splits_cv_inner = n_splits_cv_inner
        phenotypes.logreg_solver = cls.get_logreg_solver(
            logreg_solver)
        phenotypes.pca = pca

    @staticmethod
    def get_model_name(regressor, binary_classifier):
        if phenotypes.pred_scale == "continuous":
            if regressor == "lin":
                phenotypes.model_name_long = "linear regression"
                phenotypes.model_name_short = "linreg"
            elif regressor == "XGBR":
                phenotypes.model_name_long = "XGBRegressor"
                phenotypes.model_name_short = "XGBR"
        elif phenotypes.pred_scale == "binary":
            if binary_classifier == "log":
                phenotypes.model_name_long = "logistic regression"
                phenotypes.model_name_short = "log_reg"
            elif binary_classifier == "SVM":
                phenotypes.model_name_long = "support vector machine"
                phenotypes.model_name_short = "SVM"
            elif binary_classifier == "RF":
                phenotypes.model_name_long = "random forest"
                phenotypes.model_name_short = "RF"
            elif binary_classifier == "DT":
                phenotypes.model_name_long = "decision tree"
                phenotypes.model_name_short = "DT"
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

    @staticmethod
    def get_logreg_solver(logreg_solver):
        if phenotypes.pred_scale == "binary":
            if phenotypes.model_name_short == "log_reg":
                if phenotypes.penalty == "L1":
                    if logreg_solver == None:
                        return "liblinear"
                    elif logreg_solver in ("liblinear", "saga"):
                        return logreg_solver
                    else:
                        raise SystemExit("Logistic Regression with L1 penalty supports only " +
                            "solvers in ['liblinear', 'saga'], got {}.".format(logreg_solver))
                elif phenotypes.penalty == "L2":
                    if logreg_solver == None:
                        return "lbfgs"
                    elif logreg_solver in ('liblinear', 'newton-cg', 'lbfgs', 'sag', 'saga'):
                        return logreg_solver
                    else:
                        raise SystemExit("Logistic Regression with L2 penalty supports only " +
                            "solvers in ['liblinear', 'newton-cg', 'lbfgs', 'sag', 'saga'], " +
                            "got {}.".format(logreg_solver))

class Samples():

    no_samples = 0
    no_phenoypes = 0
    phenotypes = []
    take_logs = None

    kmer_length = None
    cutoff = None
    min_samples = None
    max_samples = None

    tree = None

    mash_distances_args = []
    union_output = Manager().list()

    def __init__(self, name, address, phenotypes, weight=1):
        self.name = name
        self.address = address
        self.phenotypes = phenotypes
        self.weight = weight

        Samples.no_samples += 1

    @classmethod
    def from_inputfile(cls, line):
        sample_phenotypes = {}
        name, address, phenotype_list = \
            line.split()[0], line.split()[1], line.split()[2:]
        if not all(x == "0" or x == "1" or x == "NA" for x in phenotype_list):
            phenotypes.pred_scale = "continuous"
        for i,j in zip(cls.phenotypes, phenotype_list):
            sample_phenotypes[i] = j
        return cls(name, address, sample_phenotypes)

    def get_kmer_lists(self):
        # Makes "K-mer_lists" directory where all lists are stored.
        # Generates k-mer lists for every sample in names_of_samples variable 
        # (list or dict).
        run(["mkdir", "-p", "K-mer_lists"])
        process = run(
            ["glistmaker " + self.address + " -o K-mer_lists/" + 
            self.name + "_0" + " -w " + self.kmer_length + " -c " + self.cutoff], shell=True
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
                "glistquery K-mer_lists/" + self.name + "_0_" + self.kmer_length +
                ".list -l K-mer_lists/feature_vector.list"
                ]
                , shell=True, stdout=outputfile)
            call(
                [
                "rm " + "K-mer_lists/" + self.name + "_0_" + self.kmer_length + ".list" 
                ]
                , shell=True)
        call(
            [
            "split -a 5 -d -n r/" + str(Input.num_threads) + \
            " K-mer_lists/" + self.name + "_mapped.txt " + \
            "K-mer_lists/" + self.name + "_mapped_"
            ],
            shell=True
            )
        call(["rm K-mer_lists/{}_mapped.txt".format(self.name)], shell=True)

        Input.lock.acquire()
        stderr_print.currentSampleNum.value += 1
        Input.lock.release()
        stderr_print.print_progress("samples mapped.")

    @classmethod
    def get_feature_vector(cls):
        # Feature vector loop
        iterate_to_union = [[x] for x in list(Input.samples.values())]
        for i in range(math.log(cls.no_samples, 2).__trunc__()):
            iterate_to_union = [x[0] for x in iterate_to_union]
            iterate_to_union = [
                iterate_to_union[j: j + 4 if len(iterate_to_union) < j + 4 else j + 2] for j in range(0, len(iterate_to_union), 2) if j + 2 <= len(iterate_to_union)
                ]
            with Pool(Input.num_threads) as p:
                p.map(partial(cls.get_union, round=i), iterate_to_union)
        call(["mv %s K-mer_lists/feature_vector.list" % cls.union_output[-1]], shell=True)
        [(lambda x: call(["rm -f {}".format(x)], shell=True))(union) for union in cls.union_output[:-1]]

    @classmethod
    def get_union(cls, lists_to_unite, round):
        glistcompare_args = "glistcompare -u -o K-mer_lists/" + lists_to_unite[0].name + "_" + str(round + 1) + \
            "".join([ " K-mer_lists/" + sample.name + "_" + str(round) + "_" + Samples.kmer_length + ("_union" if round > 0 else "") + ".list" \
                for sample in lists_to_unite])
        call(glistcompare_args, shell=True)
        cls.union_output.append("K-mer_lists/%s_%s_%s_union.list" % (lists_to_unite[0].name, str(round + 1), Samples.kmer_length))

    # -------------------------------------------------------------------
    # Functions for calculating the mash distances and GSC weights for
    # input samples.

    def get_mash_sketches(self):
        mash_args = "mash sketch -r " + self.address + " -o K-mer_lists/" + self.name
        process = Popen(mash_args, shell=True, stderr=PIPE, universal_newlines=True)
        for line in iter(process.stderr.readline, ''):
             stderr_print(line.strip())

    @classmethod
    def get_weights(cls):
        cls.get_mash_distances()
        cls._mash_output_to_distance_matrix(list(Input.samples.keys()), "mash_distances.mat")
        dist_mat = cls._distance_matrix_modifier("distances.mat")
        cls._distance_matrix_to_phyloxml(list(Input.samples.keys()), dist_mat)   
        cls._phyloxml_to_newick("tree_xml.txt")
        stderr_print("\x1b[1;32mCalculating the GSC weights from mash distance matrix...\x1b[0m")
        weights = cls.GSC_weights_from_newick("tree_newick.txt", normalize="mean1")
        for key, value in weights.items():
            Input.samples[key].weight = value

    @classmethod
    def get_mash_distances(cls):
        mash_args = "mash paste reference.msh K-mer_lists/*.msh"
        process = Popen(mash_args, shell=True, stderr=PIPE, universal_newlines=True)
        for line in iter(process.stderr.readline, ''):
            stderr_print(line.strip())
        call(["rm K-mer_lists/*.msh"], shell=True)
        with open("mash_distances.mat", "w+") as f1:
            call(["mash dist reference.msh reference.msh"], shell=True, stdout=f1)

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
                                "\n" + names_of_samples[counter//cls.no_samples]
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
        for i in range(len(distancematrix)):
            for j in range(len(distancematrix[i])):
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

    @classmethod
    def GSC_weights_from_newick(cls, newick_tree, normalize="sum1"):
        # Calculating Gerstein Sonnhammer Coathia weights from Newick 
        # string. Returns dictionary where sample names are keys and GSC 
        # weights are values.
        cls.tree=Tree(newick_tree, format=1)
        cls.clip_branch_lengths(cls.tree)
        cls.set_branch_sum(cls.tree)
        cls.set_node_weight(cls.tree)

        weights = {}
        for leaf in cls.tree.iter_leaves():
            weights[leaf.name] = leaf.NodeWeight
        if normalize == "mean1":
            weights = {k: v*len(weights) for k, v in weights.items()}
        return(weights)

    @staticmethod
    def clip_branch_lengths(tree, min_val=1e-9, max_val=1e9): 
        for branch in tree.traverse("levelorder"):
            if branch.dist > max_val:
                branch.dist = max_val
            elif branch.dist < min_val:
                branch.dist = min_val

    @classmethod
    def set_branch_sum(cls, tree):
        total = 0
        for child in tree.get_children():
            cls.set_branch_sum(child)
            total += child.BranchSum
            total += child.dist
        tree.BranchSum = total

    @classmethod
    def set_node_weight(cls, tree):
        parent = tree.up
        if parent is None:
            tree.NodeWeight = 1.0
        else:
            tree.NodeWeight = parent.NodeWeight * \
                (tree.dist + tree.BranchSum)/parent.BranchSum
        for child in tree.get_children():
            cls.set_node_weight(child)



class stderr_print():
    # --------------------------------------------------------
    # Functions and variables necessarry to show the progress 
    # information in standard error.

    currentSampleNum = Value("i", 0)
    currentKmerNum = Value("i", 0)
    previousPercent = Value("i", 0)

    def __init__(self,data):
        sys.stderr.write("\r\x1b[K\x1b[1;32m"+data.__str__()+"\x1b[0m")
        sys.stderr.flush()

    @classmethod
    def check_progress(cls, totalKmers, text, phenotype=""):
        currentPercent = int((cls.currentKmerNum.value/float(totalKmers))*100)
        if currentPercent > cls.previousPercent.value:
            if currentPercent != 100:
                output = f"\t{phenotype} \x1b[1;91m{currentPercent}% \x1b[1;32m{text}"
            else:
                output = f"\t{phenotype} {currentPercent}% {text}"
            cls.previousPercent.value = currentPercent
            cls(output)

    @classmethod
    def print_progress(cls, txt):
        if cls.currentSampleNum.value != Samples.no_samples:
            output = f"""\t\x1b[1;91m{cls.currentSampleNum.value}\x1b[1;32m of {Samples.no_samples} {txt}"""
        else:
            output = f"""\t{cls.currentSampleNum.value} of {Samples.no_samples} {txt}"""            
        cls(output)

class phenotypes():

    model_package = {}

    pred_scale = "binary"

    model_name_long = None
    model_name_short = None

    # Multithreading parameters
    vectors_as_multiple_input = []
    progress_checkpoint = int()
    no_kmers_to_analyse = int()

    # Filtering parameters
    pvalue_cutoff = None
    kmer_limit = None
    FDR = None
    B = None

    # Machine learning parameters
    penalty = None
    max_iter = None
    tol = None
    l1_ratio = None
    hyper_parameters = None
    alphas = None
    gammes = None
    kernel = None
    n_iter = None
    n_splits_cv_outer = None
    n_splits_cv_inner = None
    testset_size = None

    # PCA
    pca = None

    no_results = []

    def __init__(self, name):
        self.name = name
        self.no_samples = Samples.no_samples
        self.n_splits_cv_outer = None
        self.n_splits_cv_inner = None
        self.pvalues = None
        self.kmers_for_ML = {}
        self.skl_dataset = None
        self.ML_df = None
        self.ML_dict = dict()
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
        self.model_fitted = None
        self.test_output = None
        self.ttest_statistics = []
        self.ttest_pvalues = []

        #PCA
        self.PCA_df = None
        self.scaled_df = None
        self.pca_components = None
        self.pca_explained_variance_ = None
        self.pca_explained_variance_ratio_ = None

        self.metrics_dict_train = {
            "MSE": [], "CoD": [], "SpCC": [], "Sp_pval": [], "PeCC": [], "Pe_pval": [],
            "DFA": [], "Acc": [], "Sn": [], "Sp": [], "AUCROC": [], "Pr": [], "MCC": [],
            "kappa": [],"VME": [], "ME": [], "F1_sc": []
            }
        self.metrics_dict_test = {
            "MSE": [], "CoD": [], "SpCC": [], "Sp_pval": [], "PeCC": [], "Pe_pval": [],
            "DFA": [], "Acc": [], "Sn": [], "Sp": [], "AUCROC": [], "Pr": [], "MCC": [],
            "kappa": [],"VME": [], "ME": [], "F1_sc": []
            }
        self.model = None
        self.best_model = None
        # ML output file holders
        self.summary_file = None
        self.coeff_file = None
        self.model_file = None

    # -------------------------------------------------------------------
    # Functions for calculating the association test results for kmers.
    @classmethod
    def kmer_testing_setup(cls):
        if phenotypes.pred_scale == "continuous":
            sys.stderr.write("\n\x1b[1;32mConducting the k-mer specific Welch t-tests:\x1b[0m\n")
            sys.stderr.flush()
        else:
            sys.stderr.write("\n\x1b[1;32mConducting the k-mer specific chi-square tests:\x1b[0m\n")
            sys.stderr.flush()

        # Read the numbers of k-mers from feature_vector.list and delete file thereafter
        ps = Popen("glistquery K-mer_lists/feature_vector.list", shell=True, stdout=PIPE)
        output = check_output(('wc', '-l'), stdin=ps.stdout)
        ps.wait()
        cls.no_kmers_to_analyse = int(output)
        cls.progress_checkpoint = int(
            math.ceil(cls.no_kmers_to_analyse/(100*Input.num_threads))
            )
        call(["rm K-mer_lists/feature_vector.list"], shell=True)

        # Set up split up vectors as multiple input list
        for sample in Input.samples:
            cls.vectors_as_multiple_input.append(
                [
                "K-mer_lists/" + sample + "_mapped_%05d" % i \
                for i in range(Input.num_threads)
                ]
                )

    def test_kmers_association_with_phenotype(self):
        stderr_print.currentKmerNum.value = 0
        stderr_print.previousPercent.value = 0
        with Pool(Input.num_threads) as p:
            pvalues_from_all_threads = p.map(
                self.get_kmers_tested, zip(*self.vectors_as_multiple_input)
            )
        # self.pvalues = \
        #     sorted(list(chain(*pvalues_from_all_threads)))
        sys.stderr.write("\n")
        sys.stderr.flush()
        # self.concatenate_test_files(self.name)
        self.ML_df = pd.concat(pvalues_from_all_threads, axis=1)
        del pvalues_from_all_threads
        self.ML_df.columns.name = "K-mer"
        self.ML_df.index = ['chi-square statistic', 'p-value'] + list(Input.samples.keys())
        print(self.ML_df.T)
        print(sys.getsizeof(self.ML_df))
        self.ML_df.T.to_csv('pd_matrix.tsv', sep='\t')
        if self.kmer_limit:
            self.ML_df = self.ML_df.sort_values('p-value', axis=1).iloc[:,self.kmer_limit]
            self.ML_df.T.to_csv('pd_matrix_limited.tsv', sep='\t')

    def get_kmers_tested(self, split_of_kmer_lists):

        kmer_matrix = pd.DataFrame()
        counter = 0

        mt_code = split_of_kmer_lists[0][-5:]
        if phenotypes.pred_scale == "continuous":
            test_results_file = open(
                "t-test_results_" + self.name + "_" + mt_code + ".txt", "w"
                )
        else:
            test_results_file = open(
                "chi-squared_test_results_" + self.name + "_" + mt_code + ".txt", "w"
                )
        for line in zip(*[open(item) for item in split_of_kmer_lists]):
            counter += 1
            if counter%self.progress_checkpoint == 0:
                Input.lock.acquire()
                stderr_print.currentKmerNum.value += self.progress_checkpoint
                Input.lock.release()
                stderr_print.check_progress(
                    self.no_kmers_to_analyse, "tests conducted.", self.name + ": "
                )
            kmer = line[0].split()[0]
            kmer_presence_vector = [j.split()[1].strip() for j in line]

            if phenotypes.pred_scale == "binary":
                kmer, chisquare, pvalue, kmer_presence = self.conduct_chi_squared_test(
                    kmer, kmer_presence_vector,
                    test_results_file, Input.samples.values()
                    )
            elif phenotypes.pred_scale == "continuous":
                pvalue = self.conduct_t_test(
                    kmer, kmer_presence_vector,
                    test_results_file, Input.samples.values()
                    )
            if pvalue:
                kmer_matrix[kmer] = [chisquare] + [pvalue] + kmer_presence
        Input.lock.acquire()
        stderr_print.currentKmerNum.value += counter%self.progress_checkpoint
        Input.lock.release()
        stderr_print.check_progress(
            self.no_kmers_to_analyse, "tests conducted.", self.name + ": "
        )
        test_results_file.close()
        return(kmer_matrix)

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
                    y.append(sample_phenotype)
                    y_weights.append(sample.weight)
                else:
                    x.append(sample_phenotype)
                    x_weights.append(sample.weight)
                    samples_w_kmer.append(sample.name)

    @staticmethod
    def t_test(x, y, x_weights, y_weights):
        #Parametes for group containig the k-mer
        wtd_mean_y = np.average(y, weights=y_weights)
        sumofweightsy = sum(y_weights)
        ybar = np.float64(sum([i*j for i,j in zip(y, y_weights)])/sumofweightsy)
        vary = sum([i*j for i,j in zip(y_weights, (y - ybar)**2)])/(sumofweightsy-1)
        
        #Parameters for group not containig the k-mer
        wtd_mean_x = np.average(x, weights=x_weights)
        sumofweightsx = sum(x_weights)
        xbar = np.float64(sum([i*j for i,j in zip(x, x_weights)])/sumofweightsx)
        varx = sum([i*j for i,j in zip(x_weights, (x - xbar)**2)])/(sumofweightsx-1)

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
            return None, None, None, None
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
        # test_results_file.write(
        #     kmer + "\t%.2f\t%.2E\t" % chisquare_results 
        #     + str(no_samples_w_kmer)  +"\t| " + " ".join(samples_w_kmer) + "\n"
        #     )
        chisquare, pvalue = chisquare_results
        if self.B and pvalue < (self.pvalue_cutoff/self.no_kmers_to_analyse):
            return kmer, round(chisquare,2), round(pvalue,2), kmer_presence
        elif pvalue < self.pvalue_cutoff:
            return kmer, round(chisquare,2), round(pvalue,2), kmer_presence
        else:
            return None, None, None, None

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
            if sample.phenotypes[self.name] == 1:
                if (kmers_presence_vector[index] != "0"):
                    with_pheno_with_kmer += sample.weight 
                    samples_w_kmer.append(sample.name)
                else:
                    with_pheno_without_kmer += sample.weight
                    no_samples_wo_kmer += 1
            elif sample.phenotypes[self.name] == 0:
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
        if phenotypes.pred_scale == "continuous":
            test_results = "t-test_results_"
        else:
            test_results = "chi-squared_test_results_"
        self.test_output = test_results + phenotype + ".txt"
        call(
            [
            "cat " + test_results + phenotype + "_* > " +
            self.test_output
            ],
            shell=True
            )
        if not int(check_output(["wc", "-l", self.test_output]).split()[0].decode()):
            with open(self.test_output, "w+") as test_out:
                test_out.write(
                    "No k-mer had a suitable distribution to conduct the test."
                    )
            self.no_results.append(phenotype)
        for l in range(Input.num_threads):
            call(
                [
                "rm " + test_results + phenotype +
                "_%05d.txt" % l
                ],
                shell=True
                )

    # -------------------------------------------------------------------
    # Functions for filtering the k-mers based on the p-values of
    # conducted tests.
    @timer
    def get_kmers_filtered(self):
        # Filters the k-mers by their p-value achieved in statistical 
        # testing.
        phenotype = self.name
        nr_of_kmers_tested = float(len(self.pvalues))
        self.get_pvalue_cutoff(self.pvalues, nr_of_kmers_tested)
        if self.kmer_limit:
            pval_limit = float('{:.2e}'.format(self.pvalues[self.kmer_limit]))
        # while self.pvalues[self.kmer_limit-counter] == reference:
        #     counter +=1
        # pvalues_for_ML_kmers = self.pvalues[:self.kmer_limit]
        # pvalues_for_ML_kmers = [float("%.2E" % x) for x in pvalues_for_ML_kmers]
        del self.pvalues

        stderr_print.currentKmerNum.value = 0
        stderr_print.previousPercent.value = 0
        checkpoint = int(math.ceil(nr_of_kmers_tested/100))
        inputfile = open(self.test_output)
        outputfile = open("k-mers_filtered_by_pvalue_" + self.name + ".txt", "w")
        self.write_headerline(outputfile)

        counter = 0
        unique_patterns = set()
        drop_collinearity = False
        for line in inputfile:
            counter += 1
            line_to_list = line.split()
            kmer, p_val = line_to_list[0], float(line_to_list[2])
            samples_with_kmer = set(line.split("|")[1].split())
            if p_val < self.pvalue_cutoff:
                outputfile.write(line)             
                # if p_val in pvalues_for_ML_kmers:
                if drop_collinearity:
                    if line.split("|")[1] not in unique_patterns:
                        unique_patterns.add(line.split("|")[1])
                        self.kmers_for_ML[kmer] = [
                            1 if sample in samples_with_kmer else 0 for sample in Input.samples.keys()
                            ] + [p_val]
                elif self.kmer_limit and (p_val <= pval_limit):
                    self.kmers_for_ML[kmer] = [
                            1 if sample in samples_with_kmer else 0 for sample in Input.samples.keys()
                            ] + [p_val]
                elif not self.kmer_limit:
                    self.kmers_for_ML[kmer] = [
                            1 if sample in samples_with_kmer else 0 for sample in Input.samples.keys()
                            ] + [p_val]
                # pvalues_for_ML_kmers.remove(p_val)
            if counter%checkpoint == 0:
                stderr_print.currentKmerNum.value += checkpoint
                stderr_print.check_progress(
                    nr_of_kmers_tested, "k-mers filtered.", self.name + ": "
                )
        stderr_print.currentKmerNum.value += counter%checkpoint
        stderr_print.check_progress(
            nr_of_kmers_tested, "k-mers filtered.", self.name + ": "
            )
        sys.stderr.write("\n")
        sys.stderr.flush()
        if len(self.kmers_for_ML) == 0:
            outputfile.write("\nNo k-mers passed the filtration by p-value.\n")
            self.no_results.append(phenotype)
        inputfile.close()
        outputfile.close()

    def get_pvalue_cutoff(self, pvalues, nr_of_kmers_tested):
        if self.B:
            self.pvalue_cutoff = (self.pvalue_cutoff/nr_of_kmers_tested)
        elif self.FDR:
            pvalue_cutoff_by_FDR = 0
            for index, pvalue in enumerate(pvalues):
                if  (pvalue  < (
                        (index+1) 
                        / nr_of_kmers_tested) * self.pvalue_cutoff
                        ):
                    pvalue_cutoff_by_FDR = pvalue
                elif pvalue > self.pvalue_cutoff:
                    break
            self.pvalue_cutoff = pvalue_cutoff_by_FDR

    @staticmethod
    def write_headerline(outputfile):
        if phenotypes.pred_scale == "continuous":
            outputfile.write(
                "K-mer\tWelch's_t-statistic\tp-value\t+_group_mean\
                \t-_group_mean\tNo._of_samples_with_k-mer\
                \tSamples_with_k-mer\n"
                )
        elif phenotypes.pred_scale == "binary":
            outputfile.write(
                "K-mer\tChi-square_statistic\tp-value\
                \tNo._of_samples_with_k-mer\tSamples_with_k-mer\n"
                )

    def machine_learning_modelling(self):
        sys.stderr.write("\x1b[1;32m\t" + self.name + ".\x1b[0m\n")
        sys.stderr.flush()
        self.set_model()
        self.set_hyperparameters()
        self.get_outputfile_names()
        self.get_ML_dataframe()
        if self.pca:
            self.PCA_analysis()
        if phenotypes.n_splits_cv_outer:
            self.assert_n_splits_cv_outer(phenotypes.n_splits_cv_outer, self.ML_df)
            self.assert_n_splits_cv_inner(phenotypes.n_splits_cv_inner, self.ML_df)
            if phenotypes.pred_scale == "continuous":
                kf = KFold(n_splits=self.n_splits_cv_outer)               
            elif phenotypes.pred_scale == "binary":
                kf = StratifiedKFold(n_splits=self.n_splits_cv_outer)
            fold = 0
            for train_index, test_index in kf.split(
                    self.ML_df, self.ML_df['phenotype'].values
                ):
                fold += 1
                self.ML_df_train, self.ML_df_test = (
                    self.ML_df.iloc[train_index], self.ML_df.iloc[test_index]
                    )
                self.X_train, self.y_train = self.split_df(self.ML_df_train)
                self.X_test, self.y_test = self.split_df(self.ML_df_test)

                self.get_best_model()
                self.fit_model()
                self.summary_file.write(
                    "\n##### Train/test split nr.%d: #####\n" % fold
                    )
                self.cross_validation_results()
                self.summary_file.write('\nTraining set:\n')
                self.predict(self.X_train, self.y_train, self.metrics_dict_train)
                self.summary_file.write('\nTest set:\n')
                self.predict(self.X_test, self.y_test, self.metrics_dict_test)
            
            if not self.train_on_whole:
                self.summary_file.write(
                '''\n### Outputting the last model to a model file! ###\n'''
                )

            if self.pred_scale == "continuous":
                self.summary_file.write(
                    "\nMean performance metrics over all train splits: \n\n"
                    )
                self.mean_model_performance_regressor(self.metrics_dict_train)
                self.summary_file.write(
                    "\nMean performance metrics over all test splits: \n\n"
                    )
                self.mean_model_performance_regressor(self.metrics_dict_test)
            elif self.pred_scale == "binary":
                self.summary_file.write(
                    "\nMean performance metrics over all train splits: \n\n"
                    )
                self.mean_model_performance_classifier(self.metrics_dict_train)
                self.summary_file.write(
                    "\nMean performance metrics over all test splits: \n\n"
                    )
                self.mean_model_performance_classifier(self.metrics_dict_test)

        elif self.testset_size:
            if phenotypes.pred_scale == "continuous":
                stratify = None
            elif phenotypes.pred_scale == "binary":
                stratify = self.ML_df['phenotype'].values
            (
            self.ML_df_train, self.ML_df_test
            ) = train_test_split(
            self.ML_df, test_size=self.testset_size, random_state=55,
            stratify=stratify
            )
            self.X_train, self.y_train= self.split_df(self.ML_df_train)
            self.X_test, self.y_test = self.split_df(self.ML_df_test)
            self.assert_n_splits_cv_inner(
                phenotypes.n_splits_cv_inner, self.ML_df, self.y_train
                )
            self.get_best_model()
            self.fit_model()
            self.cross_validation_results()
            self.summary_file.write('\nTraining set:\n')
            self.predict(self.X_train, self.y_train, self.metrics_dict_train)
            self.summary_file.write('\nTest set:\n')
            self.predict(self.X_test, self.y_test, self.metrics_dict_test)

            if not self.train_on_whole:
                self.summary_file.write(
                '\n### Outputting the model to a file! ###\n'
                )

        if (not phenotypes.n_splits_cv_outer and not self.testset_size) or self.train_on_whole:
            if phenotypes.n_splits_cv_outer or self.testset_size:
                self.summary_file.write(
                '\nThe final output model training on the whole dataset:\n'
                )
            self.X_train, self.y_train = self.split_df(self.ML_df)
            self.assert_n_splits_cv_inner(
                phenotypes.n_splits_cv_inner, self.ML_df, self.y_train
                )
            self.get_best_model()
            self.fit_model()
            self.cross_validation_results()
            self.predict(self.X_train, self.y_train)
            if phenotypes.n_splits_cv_outer or self.testset_size:
                self.summary_file.write(
                '\n### Outputting the last model trained on whole data to a model file! ###\n'
                )
            else:                
                self.summary_file.write(
                '\n### Outputting the model to a model file! ###\n'
                )

        # Set and dump model package
        self.model_package['model'] = self.model_fitted
        self.model_package['pca'] = self.pca
        self.model_package['pred_scale'] = self.pred_scale
        joblib.dump(self.model_package, self.model_file)

        self.write_model_coefficients_to_file()

        if phenotypes.model_name_long == "decision tree":
            self.visualize_model()

        self.summary_file.close()
        self.coeff_file.close()
        self.model_file.close()

    @classmethod
    def split_df(cls, df):
        return df.iloc[:,0:-2], df.iloc[:,-1]

    def set_model(self):
        if self.pred_scale == "continuous":
            if self.model_name_short == "linreg":
                # Defining linear regression parameters

                if self.penalty == 'L1':
                    self.model = Lasso(max_iter=self.max_iter, tol=self.tol)        
                if self.penalty == 'L2':
                    self.model = Ridge(max_iter=self.max_iter, tol=self.tol)
                if self.penalty == 'elasticnet' or self.penalty == "L1+L2":
                    self.model = ElasticNet(
                        l1_ratio=self.l1_ratio, max_iter=self.max_iter, tol=self.tol
                        )
            elif self.model_name_short == "XGBR":
                self.model = xgb.XGBRegressor()
        elif self.pred_scale == "binary":
            if self.model_name_long == "logistic regression":
                #Defining logistic regression parameters
                if self.penalty == "L1":
                    self.model = LogisticRegression(
                        penalty='l1', solver=self.logreg_solver,
                        max_iter=self.max_iter, tol=self.tol
                        )        
                elif self.penalty == "L2":
                    self.model = LogisticRegression(
                        penalty='l2', solver=self.logreg_solver,
                        max_iter=self.max_iter, tol=self.tol
                        )
                elif self.penalty == "elasticnet" or "L1+L2":
                    self.model = SGDClassifier(
                        penalty='elasticnet', l1_ratio=self.l1_ratio,
                        max_iter=self.max_iter, tol=self.tol, loss='log'
                        )
            elif self.model_name_long == "support vector machine":
                self.model = SVC(
                    kernel=self.kernel, probability=True,
                    max_iter=self.max_iter, tol=self.tol
                    ) 
            elif self.model_name_long == "random forest":
                self.model = RandomForestClassifier()
            elif self.model_name_long == "decision tree":
                self.model = DecisionTreeClassifier()
            elif self.model_name_long == "Naive Bayes":
                self.model = BernoulliNB()
            elif self.model_name_short == "XGBC":
                self.model = xgb.XGBClassifier()

    def set_hyperparameters(self):
        if self.pred_scale == "continuous":
            if self.model_name_short == "linreg":
                # Defining linear regression parameters    
                self.hyper_parameters = {'alpha': self.alphas}
        elif self.pred_scale == "binary":
            if self.model_name_long == "logistic regression":
                #Defining logistic regression parameters
                if self.penalty == "L1" or "L2":
                    Cs = list(map(lambda x: 1/x, self.alphas))
                    self.hyper_parameters = {'C':Cs}
                elif penalty == "elasticnet":
                    self.hyper_parameters = {'alpha': self.alphas}
            elif self.model_name_long == "support vector machine":
                Cs = list(map(lambda x: 1/x, self.alphas))
                Gammas = list(map(lambda x: 1/x, self.gammas))
                if self.kernel == "linear":
                    self.hyper_parameters = {'C':Cs}
                if self.kernel == "rbf":
                    self.hyper_parameters = {'C':Cs, 'gamma':Gammas}
            elif self.model_name_long == "random forest":
                self.hyper_parameters = {
                    'bootstrap': [True, False],
                    'max_depth': [4, 5, 6, 7, 8, 10, 20, 100, None],
                    'max_features': [None, 'sqrt', 'log2'],
                    'min_samples_leaf': [1, 2, 4],
                    'min_samples_split': [2, 5, 10],
                    'n_estimators': [
                        10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200
                        ],
                    'criterion' :['gini', 'entropy']
                    }
            elif self.model_name_long == "decision tree":
                self.hyper_parameters = {
                    'max_depth': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    'criterion' :['gini', 'entropy']
                    }

    def get_best_model(self):
        if self.pred_scale == "continuous":
            if self.model_name_short == "linreg":
                self.best_model = GridSearchCV(
                    self.model, self.hyper_parameters, cv=self.n_splits_cv_inner
                    )
            elif self.model_name_short == "XGBR":
                self.best_model = self.model
        elif self.pred_scale == "binary":
            if self.model_name_long == "logistic regression":
                self.best_model = GridSearchCV(
                    self.model, self.hyper_parameters, cv=self.n_splits_cv_inner
                    )
            elif self.model_name_long == "support vector machine":
                if self.kernel == "linear":
                    self.best_model = GridSearchCV(
                        self.model, self.hyper_parameters, cv=self.n_splits_cv_inner
                        )
                if self.kernel == "rbf":
                    self.best_model = RandomizedSearchCV(
                        self.model, self.hyper_parameters,
                        n_iter=self.n_iter, cv=self.n_splits_cv_inner
                        )
            elif self.model_name_long == "random forest":
                self.best_model = RandomizedSearchCV(
                    self.model, self.hyper_parameters, n_iter=self.n_iter, cv=self.n_splits_cv_inner
                    )
            elif self.model_name_long == "decision tree":
                self.best_model = GridSearchCV(
                    self.model, self.hyper_parameters, cv=self.n_splits_cv_inner
                    )
            elif self.model_name_short in ("NB", "XGBC"):
                self.best_model = self.model

    def get_outputfile_names(self):
        self.summary_file = open("summary_of_" + self.model_name_short + "_analysis_" \
            + self.name + ".txt", "w")
        self.coeff_file = open("k-mers_and_coefficients_in_" + self.model_name_short \
            + "_model_" + self.name + ".txt", "w")
        self.model_file = open(self.model_name_short + "_model_" + self.name + ".pkl", "wb")

    @timer
    def get_ML_dataframe(self):
        if Input.jump_to == "modelling":
            self.ML_df = pd.read_csv(
                self.name + "_" + self.model_name_short + "_df.csv", index_col=0
                )
            self.ML_df.index = self.ML_df.index.astype(str)
            self.model_package['kmers'] = self.ML_df.columns[:-2]
        elif len(self.kmers_for_ML) == 0:
            self.summary_file.write("No k-mers passed the step of k-mer filtering for " \
                "machine learning modelling.\n")
            return
        else:
            index = list(Input.samples.keys()) + ['p_val']
            self.ML_df = pd.DataFrame(self.kmers_for_ML, index=index)
            if self.kmer_limit:
                self.ML_df.sort_values('p_val', axis=1, ascending=True, inplace=True)
                self.ML_df = self.ML_df.iloc[:,:self.kmer_limit]
            self.ML_df.drop('p_val', inplace=True)
            self.ML_df['weights'] = [
                sample.weight for sample in Input.samples.values()
                ]
            self.ML_df['phenotype'] = [
                sample.phenotypes[self.name] for sample in Input.samples.values()
                ]
            self.ML_df = self.ML_df.loc[self.ML_df.phenotype != 'NA']
            self.ML_df.phenotype = self.ML_df.phenotype.apply(pd.to_numeric)
            self.ML_df.to_csv(self.name + "_" + self.model_name_short + "_df.csv")
            self.model_package['kmers'] = self.ML_df.columns[:-2]

    @timer
    def PCA_analysis(self):

        # Strandardization
        df_to_scale = self.ML_df.drop(['phenotype', 'weights'], axis=1)
        scaler = StandardScaler()
        scaler.fit(df_to_scale)
        scaled_data = scaler.transform(df_to_scale)

        # PCA transformation
        pca = PCA()
        pca.fit(scaled_data)
        self.PCA_df = pd.DataFrame(
            pca.transform(scaled_data),
            index=self.ML_df.index,
            )
        self.PCA_df.columns = [
            'PC_' + str(i) for i in  range(1, 1 + self.PCA_df.shape[1])
            ]

        del scaled_data
        np.set_printoptions(threshold=sys.maxsize)
        pd.set_option('display.max_rows', 1000)

        # Filter PCs by explained variance
        PCs_to_keep = pca.explained_variance_ > 0.9
        self.PCA_df = self.PCA_df.loc[:, PCs_to_keep]
        self.pca_components_ = pca.components_[PCs_to_keep]
        self.pca_explained_variance_ = pca.explained_variance_[PCs_to_keep]
        self.pca_explained_variance_ratio_ = pca.explained_variance_ratio_[PCs_to_keep]

        # Conduct the t-test analysis between PCs and phenotypes
        for idx, column in enumerate(self.PCA_df):
            x = self.PCA_df[column][self.ML_df.phenotype == 1].values
            y = self.PCA_df[column][self.ML_df.phenotype == 0].values
            x_weights = self.ML_df['weights'][self.ML_df.phenotype == 1].values
            y_weights = self.ML_df['weights'][self.ML_df.phenotype == 0].values
            t_statistic, pvalue, mean_x, mean_y = self.t_test(
                    x, y, x_weights, y_weights
                )
            # if pvalue < 0.05/self.PCA_df.shape[1]:
            if pvalue < 0.01:
                PCs_to_keep[idx] = True
                self.ttest_statistics.append(t_statistic)
                self.ttest_pvalues.append(pvalue)
            else:
                PCs_to_keep[idx] = False
        self.model_package['PCs_to_keep'] = PCs_to_keep

        # Filter PCs by association with phenotype
        PCs_to_keep = PCs_to_keep[:self.PCA_df.shape[1]]
        self.PCA_df = self.PCA_df.loc[:, PCs_to_keep]
        self.pca_components_ = self.pca_components_[PCs_to_keep]
        self.pca_explained_variance_ = self.pca_explained_variance_[PCs_to_keep]
        self.pca_explained_variance_ratio_ = self.pca_explained_variance_ratio_[PCs_to_keep]        

        # Set up the outputs
        self.ML_df = pd.concat([self.PCA_df, self.ML_df.iloc[:,-2:]], axis=1)
        self.model_package['scaler'] = scaler
        self.model_package['pca_model'] = pca

    def fit_model(self):
        if self.pred_scale == "continuous":
            if self.model_name_short == "linreg":
                if self.penalty in ("L1", "elasticnet"):
                    self.model_fitted = self.best_model.fit(self.X_train.values, self.y_train.values.flatten())
                elif self.penalty == "L2":
                    self.model_fitted = self.best_model.fit(self.X_train.values, self.y_train.values.flatten())
            elif self.model_name_short == "XGBR":
                self.model_fitted = self.best_model.fit(self.X_train.values, self.y_train.values.flatten())
        elif self.pred_scale == "binary":
            if self.model_name_short == "XGBC":
                self.model_fitted = self.best_model.fit(
                    self.X_train.values, self.y_train.values.flatten()
                    )
            else:
                self.model_fitted = self.best_model.fit(
                    self.X_train, self.y_train.values.flatten()
                    )


    def cross_validation_results(self):
        if self.model_name_short not in ("NB", "XGBC", "XGBR"):
            self.summary_file.write('Parameters:\n%s\n\n' % self.model)
            if self.pred_scale == "continuous":
                self.summary_file.write("Grid scores (R2 score) on development set: \n")
            elif self.pred_scale == "binary":
                self.summary_file.write("Grid scores (mean accuracy) on development set: \n")
            means = self.model_fitted.cv_results_['mean_test_score']
            stds = self.model_fitted.cv_results_['std_test_score']
            params = self.model_fitted.cv_results_['params']
            for mean, std, param in zip(
                    means, stds, params
                    ):
                self.summary_file.write(
                    "%0.3f (+/-%0.03f) for %r \n" % (mean, std * 2, param)
                    )
            self.summary_file.write("\nBest parameters found on development set: \n")
            for key, value in self.model_fitted.best_params_.items():
                self.summary_file.write(key + " : " + str(value) + "\n")

    def predict(self, dataset, labels, metrics_dict=None):
        predictions = self.model_fitted.predict(dataset.values)
        self.summary_file.write("\nModel predictions on samples:\nSample_ID " \
            "Acutal_phenotype Predicted_phenotype\n")
        for index, row in dataset.iterrows():
                self.summary_file.write('%s %s %s\n' % (
                    index, labels.loc[index],
                    self.model_fitted.predict(row.values.reshape(1, -1))[0]
                    ))
        self.summary_file.write('\n')

        if self.pred_scale == "continuous":
            self.model_performance_regressor(dataset, labels.values.flatten(), predictions, metrics_dict)
        elif self.pred_scale == "binary":
            self.model_performance_classifier(dataset, labels.values.flatten(), predictions, metrics_dict)

    def model_performance_regressor(self, dataset, labels, predictions, metrics_dict):

        MSE = mean_squared_error(labels, predictions).round(2)
        self.summary_file.write('\nMean squared error: %s\n' % MSE)
        if metrics_dict:
            metrics_dict["MSE"].append(MSE)

        CoD = round(self.model_fitted.score(dataset.values, labels), 2)
        self.summary_file.write("The coefficient of determination:"
            + " %s\n" % CoD)
        if metrics_dict:
            metrics_dict["CoD"].append(CoD)

        SpCC, Sp_pval = map(lambda x: round(x, 2), stats.spearmanr(labels, predictions))
        self.summary_file.write("The Spearman correlation coefficient and p-value:" \
            " %s, %s \n" % (SpCC, Sp_pval))
        if metrics_dict:
            metrics_dict["SpCC"].append(SpCC)
            metrics_dict["Sp_pval"].append(Sp_pval)

        PeCC, Pe_pval = map(lambda x: round(x, 2), stats.pearsonr(labels, predictions))
        self.summary_file.write("The Pearson correlation coefficient and p-value: " \
                " %s, %s \n" % (PeCC, Pe_pval))
        if metrics_dict:
            metrics_dict["PeCC"].append(PeCC)
            metrics_dict["Pe_pval"].append(Pe_pval)

        DFA = self.within_1_tier_accuracy(labels, predictions)
        self.summary_file.write(
            "The plus/minus 1 dilution factor accuracy (for MICs):" \
            " %s \n\n" % DFA
            )
        if metrics_dict:
            metrics_dict["DFA"].append(DFA)

    def mean_model_performance_regressor(self, metrics_dict):

        MSE = np.mean(metrics_dict["MSE"]).round(2)
        self.summary_file.write('\nMean squared error: %s\n' % MSE)

        CoD = np.mean(metrics_dict["CoD"]).round(2)
        self.summary_file.write("The coefficient of determination:"
            + " %s\n" % CoD)

        SpCC = np.mean(metrics_dict["SpCC"]).round(2)
        Sp_pval = np.mean(metrics_dict["Sp_pval"]).round(2)
        self.summary_file.write("The Spearman correlation coefficient and p-value:" \
            " %s, %s \n" % (SpCC, Sp_pval))

        PeCC = np.mean(metrics_dict["PeCC"]).round(2)
        Pe_pval = np.mean(metrics_dict["Pe_pval"]).round(2)
        self.summary_file.write("The Pearson correlation coefficient and p-value: " \
                " %s, %s \n" % (PeCC, Pe_pval))

        DFA = np.mean(metrics_dict["DFA"]).round(2)
        self.summary_file.write(
            "The plus/minus 1 dilution factor accuracy (for MICs):" " %s \n\n" % DFA
            )

    def model_performance_classifier(self, dataset, labels, predictions, metrics_dict):

        F1_sc = f1_score(labels, predictions).round(2)
        self.summary_file.write("F1-score of positive class: %s\n" % F1_sc)
        if metrics_dict:
            metrics_dict["F1_sc"].append(F1_sc)

        Acc = self.model_fitted.score(dataset, labels).round(2)
        self.summary_file.write("Mean accuracy: %s\n" % Acc)
        if metrics_dict:
            metrics_dict["Acc"].append(Acc)

        Sn = recall_score(labels, predictions).round(2)
        self.summary_file.write("Sensitivity: %s\n" % Sn)
        if metrics_dict:
            metrics_dict["Sn"].append(Sn)

        Sp = recall_score(
                    list(map(lambda x: 1 if x == 0 else 0, labels)), 
                    list(map(lambda x: 1 if x == 0 else 0, predictions))
                    ).round(2)
        self.summary_file.write("Specificity: %s\n" % Sp)
        if metrics_dict:
            metrics_dict["Sp"].append(Sp)

        AUCROC = roc_auc_score(labels, predictions, average="micro").round(2)
        self.summary_file.write("AUC-ROC: %s\n" % AUCROC)
        if metrics_dict:
            metrics_dict["AUCROC"].append(AUCROC)

        Pr = average_precision_score(
                labels, self.model_fitted.predict_proba(dataset)[:,1]
                ).round(2)
        self.summary_file.write("Average precision: %s\n" % Pr)
        if metrics_dict:
            metrics_dict["Pr"].append(Pr)

        MCC = round(matthews_corrcoef(labels, predictions), 2)
        self.summary_file.write("MCC: %s\n" % MCC)
        if metrics_dict:
            metrics_dict["MCC"].append(MCC)

        kappa = cohen_kappa_score(labels, predictions).round(2)
        self.summary_file.write("Cohen kappa: %s\n" % kappa)
        if metrics_dict:
            metrics_dict["kappa"].append(kappa)

        VME = self.VME(labels, predictions)
        self.summary_file.write("Very major error rate: %s\n" % VME)
        if metrics_dict:
            metrics_dict["VME"].append(VME)

        ME = self.ME(labels, predictions)
        self.summary_file.write("Major error rate: %s\n" % ME)
        if metrics_dict:
            metrics_dict["ME"].append(ME)

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

    def mean_model_performance_classifier(self, metrics_dict):

        F1_sc = np.mean(metrics_dict["F1_sc"]).round(2)
        self.summary_file.write("F1-score of positive class: %s\n" % F1_sc)

        Acc = np.mean(metrics_dict["Acc"]).round(2)
        self.summary_file.write("Mean accuracy: %s\n" % Acc)

        Sn = np.mean(metrics_dict["Sn"]).round(2)
        self.summary_file.write("Sensitivity: %s\n" % Sn)

        Sp = np.mean(metrics_dict["Sp"]).round(2)
        self.summary_file.write("Specificity: %s\n" % Sp)

        AUCROC = np.mean(metrics_dict["AUCROC"]).round(2)
        self.summary_file.write("AUC-ROC: %s\n" % AUCROC)

        Pr = np.mean(metrics_dict["Pr"]).round(2)
        self.summary_file.write("Average precision: %s\n" % Pr)

        MCC = np.mean(metrics_dict["MCC"]).round(2)
        self.summary_file.write("MCC: %s\n" % MCC)

        kappa = np.mean(metrics_dict["kappa"]).round(2)
        self.summary_file.write("Cohen kappa: %s\n" % kappa)

        VME = np.mean(metrics_dict["VME"]).round(2)
        self.summary_file.write("Very major error rate: %s\n" % VME)

        ME = np.mean(metrics_dict["ME"]).round(2)
        self.summary_file.write("Major error rate: %s\n" % ME)           

    def write_model_coefficients_to_file(self):
        if self.pca == True:
            self.coeff_file.write(
                "PC\tcoef._in_" + self.model_name_short + \
                "_model\texplained_variance\texplained_variance_ratio" + \
                "\tt-test_statistic\tt-test_pvalue\n"
                )
        else:
            self.coeff_file.write(
                "K-mer\tcoef._in_" + self.model_name_short + \
                "_model\tNo._of_samples_with_k-mer\tSamples_with_k-mer\n"
                )
        self.ML_df.drop(['phenotype', 'weights'], axis=1, inplace=True)
        if self.model_name_short == "linreg":
            dself.ML_df.loc['coefficient'] = \
                self.model_fitted.best_estimator_.coef_
        elif self.model_name_short in ("RF", "DT"):
            self.ML_df.loc['coefficient'] = \
                self.model_fitted.best_estimator_.feature_importances_
        elif self.model_name_short in ("XGBR", "XGBC"):
            self.ML_df.loc['coefficient'] = \
                self.model_fitted.feature_importances_
        elif self.model_name_short in ("SVM", "log_reg"):
            if self.kernel != "rbf":
                self.ML_df.loc['coefficient'] = \
                    self.model_fitted.best_estimator_.coef_[0]

        for idx, predictor in enumerate(self.ML_df):
            # Get coefficients
            if self.kernel == "rbf" or self.model_name_short == "NB":
                coef = "NA"
            else:
                coef = self.ML_df.loc['coefficient', predictor]

            if self.pca:
                self.coeff_file.write(
                    f"{predictor}\t{coef}\t{self.pca_explained_variance_[idx]}" + \
                    f"\t{self.pca_explained_variance_ratio_[idx]}" + \
                    f"\t{self.ttest_statistics[idx]}\t{self.ttest_pvalues[idx]}\n"
                    )
            else:
                samples_with_kmer = self.ML_df.index[self.ML_df[predictor] == 1].tolist()
                self.coeff_file.write(
                    f"{predictor}\t{coef}\t{len(samples_with_kmer)}\t| {' '.join(samples_with_kmer)}\n"
                    )

    def visualize_model(self):
        plt.figure(figsize=(15,10))
        plot = plot_tree(
            self.model_fitted.best_estimator_, feature_names=self.ML_df.columns,
            filled=True, 
            )
        plt.savefig(self.model_name_short + "_model_" + self.name + '_plot.png')

    # ---------------------------------------------------------
    # Self-implemented performance measure functions
    @staticmethod
    def VME(targets, predictions):
        # Function to calculate the very major error (VME) rate
        VMEs = 0
        for item in zip(targets, predictions):
            if item[0] == 1 and item[1] == 0:
                VMEs += 1
        VME = round(float(VMEs)/len(targets), 2)
        return VME

    @staticmethod
    def ME(targets, predictions):
        # Function to calculate the major error (ME) rate
        MEs = 0
        for item in zip(targets, predictions):
            if item[0] == 0 and item[1] == 1:
                 MEs += 1
        ME = round(float(MEs)/len(targets), 2)
        return ME

    @staticmethod
    def within_1_tier_accuracy(targets, predictions):
        # Calculate the plus/minus one dilution factor accuracy
        # for predicted antibiotic resistance values.
        within_1_tier = 0
        for item in zip(targets, predictions):
            if abs(item[0]-item[1]) <= 1:
                within_1_tier +=1
        accuracy = round(float(within_1_tier)/len(targets), 2)
        return accuracy

    def assert_n_splits_cv_outer(self, n_splits_cv_outer, ML_df):
        if self.pred_scale == "continuous" and (n_splits_cv_outer > self.no_samples // 2):
            self.n_splits_cv_outer = self.no_samples // 2
            sys.stderr.write("\x1b[1;33mWarning! The 'n_splits_cv_outer' parameter is too high to \n" \
                    "leave the required 2 samples into test set for each split!\x1b[0m\n")
            sys.stderr.write("\x1b[1;33mSetting number of train/test splits equal to " \
                + str(self.n_splits_cv_outer) + "!\x1b[0m\n\n")
        elif self.pred_scale == "binary" and np.min(np.bincount(ML_df['phenotype'].values)) < n_splits_cv_outer:
            self.n_splits_cv_outer = np.min(np.bincount(ML_df['phenotype'].values))
            sys.stderr.write("\x1b[1;33mSetting number of train/test splits \
                equal to minor phenotype count - " + str(self.n_splits_cv_outer) + "!\x1b[0m\n\n")
        else:
            self.n_splits_cv_outer = n_splits_cv_outer

    def assert_n_splits_cv_inner(self, n_splits_cv_inner, ML_df, y_train=None):
        if self.pred_scale == "continuous":
            if self.n_splits_cv_outer:
                min_cv_inner = self.no_samples - math.ceil(self.no_samples / self.n_splits_cv_outer)
            else:
                min_cv_inner = len(y_train)
        elif self.pred_scale == "binary":
            if self.n_splits_cv_outer:
                min_class = np.min(np.bincount(ML_df['phenotype'].values))
                min_cv_inner = (min_class - math.ceil(min_class / self.n_splits_cv_outer))
            else:
                min_cv_inner = np.min(np.bincount(y_train))
        self.n_splits_cv_inner = np.min([min_cv_inner, n_splits_cv_inner])

    # Assembly methods
    def ReverseComplement(self, kmer):
        # Returns the reverse complement of kmer
        seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
        return("".join([seq_dict[base] for base in reversed(kmer)]))

    def string_set(self, string_list):
        # Removes subsequences from kmer_list
        return set(i for i in string_list
                   if not any(i in s for s in string_list if i != s))

    def overlap(self, a, b, min_length=3):
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

    def pick_overlaps(self, reads, min_olap):
        # Takes kmer_list as an Input. Generates pairwise permutations of 
        # the kmers in kmer list. Finds the overlap of each pair. Returns 
        # the lists of kmers and overlap lengths of the pairs which overlap
        # by min_olap or more nucleotides.
        reada, readb, olap_lens = [], [], []
        for a, b in permutations(reads, 2):
            olap_len = self.overlap(a, b, min_length=min_olap)
            if olap_len > 0:
                reada.append(a)
                readb.append(b)
                olap_lens.append(olap_len)
        return reada, readb, olap_lens

    def kmer_assembler(self, min_olap=None):
        # Assembles the k-mers in kmer_list which overlap by at least 
        # min_olap nucleotides.
        if min_olap == None:
            min_olap = int(Samples.kmer_length)-1
        assembled_kmers = []

        # Adding the reverse-complement of each k-mer
        kmer_list = list(self.kmers_for_ML) + list(map(
            self.ReverseComplement, self.kmers_for_ML
            ))

        # Find the overlaping k-mers
        kmers_a, kmers_b, olap_lens = self.pick_overlaps(kmer_list, min_olap)

        while olap_lens != []:
            set_a = set(kmers_a)
            set_b = set(kmers_b)

            # Picking out the assembled k-mers which have no sufficient
            # overlaps anymore.
            for item in kmer_list:
                if (item not in set_a and item not in set_b
                        and self.ReverseComplement(item) not in assembled_kmers):
                    assembled_kmers.append(item)

            # Generating new kmer_list, where overlaping elements from previous
            # kmer_list are assembled.
            kmer_list = []
            for i, olap in enumerate(olap_lens):
                kmer_list.append(kmers_a[i] + kmers_b[i][olap:])

            # Removing substrings of other elements from kmer_list.
            kmer_list = list(self.string_set(kmer_list))

            # Find the overlaping elements in new generated kmer_list.
            kmers_a, kmers_b, olap_lens = self.pick_overlaps(kmer_list, min_olap)

        for item in kmer_list:
            # Picking out the assembled k-mers to assembled_kmers set.
            if (self.ReverseComplement(item) not in assembled_kmers):
                assembled_kmers.append(item)
        return(assembled_kmers)

    def assembling(self):
        # Assembles the input k-mers and writes assembled sequences
        # into "assembled_kmers.txt" file in FastA format.
        #Open files to write the results of k-mer assembling
        
        f1 = open("assembled_kmers_" + self.name + ".fasta", "w+")
        sys.stderr.write("\x1b[1;32m\t" + self.name + " data.\x1b[0m\n")
        sys.stderr.flush()
        
        assembled_kmers = sorted(
            self.kmer_assembler(), key = len
            )[::-1]
        for i, item in enumerate(assembled_kmers):
            f1.write(">seq_" + str(i+1) + "_length_" 
                + str(len(item)) + "\n" + item + "\n")
        f1.close()

def modeling(args):
    # The main function of "phenotypeseeker modeling"

    sys.stderr.write("\x1b[1;1;101m######                   PhenotypeSeeker                   ######\x1b[0m\n")
    sys.stderr.write("\x1b[1;1;101m######                      modeling                       ######\x1b[0m\n\n")

    # Processing the input data
    Input.get_input_data(args.inputfile, args.take_logs, args.mpheno)
    Input.Input_args(
        args.alphas, args.alpha_min, args.alpha_max, args.n_alphas,
        args.gammas, args.gamma_min, args.gamma_max, args.n_gammas,
        args.min, args.max, args.kmer_length, args.cutoff,
        args.num_threads, args.pvalue, args.n_kmers, args.FDR, 
        args.Bonferroni, args.binary_classifier, args.regressor, 
        args.penalty, args.max_iter, args.tolerance, args.l1_ratio,
        args.n_splits_cv_outer, args.kernel, args.n_iter, args.n_splits_cv_inner,
        args.testset_size, args.train_on_whole, args.logreg_solver, args.jump_to,
        args.pca
        )

    if not Input.jump_to:
        #  Operations with samples
        sys.stderr.write("\x1b[1;32mGenerating the k-mer lists for input samples:\x1b[0m\n")
        sys.stderr.flush()

        with Pool(Input.num_threads) as p:
            p.map(
                lambda x: x.get_kmer_lists(), Input.samples.values()
            )

        sys.stderr.write("\n\x1b[1;32mGenerating the k-mer feature vector.\x1b[0m\n")
        sys.stderr.flush()
        Samples.get_feature_vector()
        sys.stderr.write("\x1b[1;32mMapping samples to the feature vector space:\x1b[0m\n")
        sys.stderr.flush()
        stderr_print.currentSampleNum.value = 0
        with Pool(Input.num_threads) as p:
            p.map(
                lambda x: x.map_samples(), Input.samples.values()
            )
        if not args.no_weights:
            mash_files = ["distances.mat", "reference.msh", "mash_distances.mat"]
            for mash_file in mash_files:
                if os.path.exists(mash_file):
                    os.remove(mash_file)
                    sys.stderr.write("\n\x1b[1;32mDeleting the existing " + mash_file + " file...\x1b[0m")
            sys.stderr.write("\n\x1b[1;32mEstimating the Mash distances between samples...\x1b[0m\n")
            sys.stderr.flush()
            with Pool(Input.num_threads) as p:
                p.map(
                    lambda x: x.get_mash_sketches(), Input.samples.values()
                )
            Samples.get_weights()

        # Analyses of phenotypes
        phenotypes.kmer_testing_setup()
        list(map(
            lambda x:  x.test_kmers_association_with_phenotype(), 
            Input.phenotypes_to_analyse.values()
            ))

        # Remove phenotypes with no results
        Input.pop_phenos_out_of_kmers()

        sys.stderr.write("\x1b[1;32mFiltering the k-mers by p-value:\x1b[0m\n")
        sys.stderr.flush()
        list(map(
            lambda x:  x.get_kmers_filtered(), 
            Input.phenotypes_to_analyse.values()
            ))
        # Remove phenotypes with no results
        Input.pop_phenos_out_of_kmers()

    if not Input.jump_to or ''.join(i for i, _ in groupby(Input.jump_to)) in ["modeling"]:
        sys.stderr.write("\x1b[1;32mGenerating the " + phenotypes.model_name_long + " model for phenotype: \x1b[0m\n")
        sys.stderr.flush()
        with Pool(Input.num_threads) as p:
            p.map(
                lambda x: x.machine_learning_modelling(),
                Input.phenotypes_to_analyse.values()
            )

    call(['rm', '-rf', 'K-mer_lists'])

    if args.assembly:
        sys.stderr.write("\x1b[1;32mAssembling the k-mers used in modeling of: \x1b[0m\n")
        sys.stderr.flush()
        with Pool(Input.num_threads) as p:
            p.map(
                lambda x: x.assembling(),
                Input.phenotypes_to_analyse.values()
            )

    sys.stderr.write("\n\x1b[1;1;101m######          PhenotypeSeeker modeling finished          ######\x1b[0m\n")
