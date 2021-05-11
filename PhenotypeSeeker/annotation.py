__author__ = "Erki Aun"
__version__ = "0.7.0"
__maintainer__ = "Erki Aun"
__email__ = "erki.aun@ut.ee"

import os

from collections import OrderedDict

class Input():

    samples = OrderedDict()
    phenotypes_to_analyse = OrderedDict()
    pool = None
    lock = None

    jump_to = None
    num_threads = 8
    
    @classmethod
    def get_input_data(cls, inputfilename):
        # Read the data from inputfile into "samples" directory
        Samples.take_logs = take_logs
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

def annotation(args):
	Input.get_input_data(args.inputfile)
	print(Input.samples)