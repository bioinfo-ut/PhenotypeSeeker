__author__ = "Erki Aun"
__version__ = "0.7.0"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

import os

from collections import OrderedDict

from multiprocess import Manager, Pool, Value

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

def annotation(args):
	Input.get_input_data(args.inputfile)
	for key, value in Input.samples:
		print(value.address)