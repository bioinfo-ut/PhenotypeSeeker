__author__ = "Erki Aun"
__version__ = "1.0.0"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

import os

from subprocess import run

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

    def get_kmer_lists(self):
        # Makes "K-mer_lists" directory where all lists are stored.
        # Generates k-mer lists for every sample in names_of_samples variable 
        # (list or dict).
        run(["mkdir", "-p", "K-mer_lists"])
        process = run(
            ["glistmaker " + self.address + " -o K-mer_lists/" + 
            self.name + "_0" + " -w " + self.kmer_length + "--index"], shell=True
            )
        Input.lock.acquire()
        stderr_print.currentSampleNum.value += 1
        Input.lock.release()
        stderr_print.print_progress("lists generated.")

def annotation(args):
    Input.get_input_data(args.inputfile)
    sys.stderr.write("\x1b[1;32mGenerating the k-mer lists in input samples:\x1b[0m\n")
    with Pool(Input.num_threads) as p:
        p.map(
            lambda x: x.get_kmer_indexes(),
            Input.samples.values()
        )
    # with Pool(Input.num_threads) as p:
    #     p.map(
    #         lambda x: x.call_prokka(),
    #         Input.samples.values()
    #     )
    sys.stderr.write("\x1b[1;32mGenerating the k-mer indexes in input samples:\x1b[0m\n")