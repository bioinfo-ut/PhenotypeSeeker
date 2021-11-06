__author__ = "Erki Aun"
__version__ = "1.0.0"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

import os
import sys
import joblib
import glob

from subprocess import run, PIPE, DEVNULL
from collections import OrderedDict
from multiprocess import Manager, Pool, Value

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
    def print_progress(cls, txt):
        if cls.currentSampleNum.value != Samples.no_samples:
            output = f"""\t\x1b[1;91m{cls.currentSampleNum.value}\x1b[1;32m of {Samples.no_samples} {txt}"""
        else:
            output = f"""\t{cls.currentSampleNum.value} of {Samples.no_samples} {txt}"""            
        cls(output)

class Input():

    samples = OrderedDict()
    pool = None
    lock = None
    num_threads = 8

    lock = Manager().Lock()

    kmers = None
    kmer_lenght = None
    
    @classmethod
    def get_input_data(cls, inputfilename):
        # Read the data from inputfile into "samples" directory
        with open(inputfilename) as inputfile:
            header = inputfile.readline().split()
            if os.path.exists(header[1]):
                inputfile.seek(0)
            for line in inputfile:
                if line.strip():
                    sample_name = line.split()[0]
                    cls.samples[sample_name] = (
                        Samples.from_inputfile(line)
                        )

    @classmethod
    def get_kmers(cls, model_pkg):
        model_pkg = joblib.load(model_pkg)
        if model_pkg['LR']:
            cls.kmers = model_pkg['kmers'].loc[model_pkg['kmers_to_keep']]
        else:
            cls.kmers = model_pkg['kmers']
        cls.kmer_length = str(len(cls.kmers.index.values[0]))

class Samples():

    no_samples = 0

    kmer_length = None
    cutoff = None
    min_samples = None
    max_samples = None

    tree = None

    vectors_as_multiple_input = Manager().list()

    def __init__(self, name, address):
        self.name = name
        self.address = address
    
        Samples.no_samples += 1

    @classmethod
    def from_inputfile(cls, line):
        name, address = line.split()[0], line.split()[1]
        return cls(name, address)

    def call_prokka(self):
        run([f"/storage8/erkia/prokka/bin/prokka --kingdom Bacteria --outdir prokka/prokka_{self.name} \
            --genus Enterococcus --locustag {self.name} {self.address}"], shell=True)
        Input.lock.acquire()
        stderr_print.currentSampleNum.value += 1
        Input.lock.release()
        stderr_print.print_progress("genomes annotated.")

    def get_kmer_indexes(self):
        # Makes "K-mer_lists" directory where all lists are stored.
        # Generates k-mer lists for every sample in names_of_samples variable 
        # (list or dict).
        run(["mkdir", "-p", "K-mer_lists"])
        process = run(
            ["glistmaker " + self.address + " -o K-mer_lists/" + 
            self.name + " -w " + Input.kmer_length + " --index"], shell=True,
            stderr=DEVNULL, stdout=DEVNULL
            )
        Input.lock.acquire()
        stderr_print.currentSampleNum.value += 1
        Input.lock.release()
        stderr_print.print_progress("lists generated.")

class annotate():

    genome_annotations = {}
    kmer_annotations = OrderedDict()

    @classmethod
    def read_in_prokka_results(cls):
        for sample in Input.samples.values():
            with open(glob.glob(f"prokka/prokka_{sample.name}/PROKKA*.gff")[0], 'r') as prokka_res:
                for line in prokka_res:
                    if f"ID={sample.name}" in line:
                        line2list = line.split('\t')
                        contig = line2list[0]
                        strand = line2list[6]
                        if strand == "+":
                            gene_start = int(line2list[3])
                            gene_end = int(line2list[4])
                        elif strand == "-":
                            gene_start = int(line2list[4])
                            gene_end = int(line2list[3])                          
                        product_name = line2list[-1].strip().split("product=")[-1]
                        if 'gene' in line2list[-1]:
                            gene_name = line2list[-1].split("gene=")[-1].split(";")[0]
                        else:
                            gene_name = "-"
                        data = {'gene_start': gene_start, 'gene_name': gene_name, 'gene_end': gene_end,
                                'strand': strand, 'product_name': product_name
                            }
                        if sample.name not in cls.genome_annotations:
                            cls.genome_annotations[sample.name] = {contig : {
                                gene_start : data,
                                gene_end : data
                                }}
                        else:
                            if contig not in cls.genome_annotations[sample.name]:
                                cls.genome_annotations[sample.name][contig] = {
                                    gene_start : data,
                                    gene_end : data
                                }
                            else:
                                cls.genome_annotations[sample.name][contig][gene_start] = data
                                cls.genome_annotations[sample.name][contig][gene_end] = data

    @classmethod
    def get_kmer_annotations(cls):
        for kmer, strains in Input.kmers.items():
            for strain in strains:
                contig_mapper = {}
                query_seqs = run(
                        ["glistquery", "--sequences",
                        f"K-mer_lists/{strain}_{Input.kmer_length}.index"
                        ]
                        , capture_output=True, text=True)
                for line in query_seqs.stdout.strip().split("\n"):
                    contig_mapper[line.split()[1]] = line.split()[2]
                returncode = -1
                while returncode != 0:
                    indexes = run(
                        ["glistquery", "--locations", "-q", kmer,
                        f"K-mer_lists/{strain}_{Input.kmer_length}.index"
                        ]
                        , capture_output=True, text=True)
                    returncode = indexes.returncode
                for line in indexes.stdout.strip().split("\n")[1:]:
                    _, contig, pos, _ = line.split()
                    cls.annotate_kmers(
                        kmer, strain, contig_mapper[contig], int(pos)+1)


    @classmethod
    def annotate_kmers(cls, kmer, strain, contig, pos):
        # Find the nearest position
        print(kmer, strain, contig, pos)
        relative_pos = "-"
        gene = "-"
        product = "-"
        if contig in cls.genome_annotations[strain]:
            nearest = min(cls.genome_annotations[strain][contig], key=lambda x:abs(x-pos))
            print(cls.genome_annotations[strain][contig][nearest])
            gene = cls.genome_annotations[strain][contig][nearest]['gene_name']
            product = cls.genome_annotations[strain][contig][nearest]['product_name']
            if cls.genome_annotations[strain][contig][nearest]['strand'] == '+':
                if (pos >= cls.genome_annotations[strain][contig][nearest]['gene_start'] and
                   pos <= cls.genome_annotations[strain][contig][nearest]['gene_end']):
                    relative_pos = 'inside'
                elif pos < cls.genome_annotations[strain][contig][nearest]['gene_start']:
                    relative_pos = 'preceding'
                elif pos > cls.genome_annotations[strain][contig][nearest]['gene_end']:
                    relative_pos = 'succeeding'
            elif cls.genome_annotations[strain][contig][nearest]['strand'] == '-':
                if (pos <= cls.genome_annotations[strain][contig][nearest]['gene_start'] and
                   pos >= cls.genome_annotations[strain][contig][nearest]['gene_end']):
                    relative_pos = 'inside'
                elif pos > cls.genome_annotations[strain][contig][nearest]['gene_start']:
                    relative_pos = 'preceding'
                elif pos < cls.genome_annotations[strain][contig][nearest]['gene_end']:
                    relative_pos = 'succeeding'
        if f"{kmer}\t{relative_pos}\t{gene}\t{product}" not in cls.kmer_annotations:
            cls.kmer_annotations[f"{kmer}\t{relative_pos}\t{gene}\t{product}"] = [strain]
        else:
            cls.kmer_annotations[f"{kmer}\t{relative_pos}\t{gene}\t{product}"].append(strain)

    @classmethod
    def write_results(cls):
        with open('kmer_annotations.txt', 'w') as out:
            out.write(f'kmer\trelative_position\tgene\tproduct\tsupporting_strains')
            for key, value in cls.kmer_annotations.items():
                out.write(f"key\t{' '.join(value)}")


def annotation(args):
    Input.get_input_data(args.inputfile)
    Input.get_kmers(args.model_file)
    # sys.stderr.write("\x1b[1;32mGenerating the k-mer indexes in input samples:\x1b[0m\n")
    # with Pool(Input.num_threads) as p:
    #     p.map(
    #         lambda x: x.get_kmer_indexes(),
    #         Input.samples.values()
    #     )
    # sys.stderr.write("\x1b[1;32mAnnotating the k-mer genomes with prokka:\x1b[0m\n")
    # stderr_print.currentSampleNum.value += 0
    # with Pool(Input.num_threads) as p:
    #     p.map(
    #         lambda x: x.call_prokka(),
    #         Input.samples.values()
    #     )
    sys.stderr.write("\x1b[1;32m\nReading in prokka annotations.\x1b[0m")
    annotate.read_in_prokka_results()
    sys.stderr.write("\x1b[1;32m\nAnnotating the k-mers:\x1b[0m\n")
    annotate.get_kmer_annotations()
    sys.stderr.write("\x1b[1;32m\nWriting results to 'kmer_annotations.txt':\x1b[0m\n")
    annotate.write_results()




