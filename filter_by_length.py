#!/usr/bin/env python3

###################################
# script:       FilterByLength.py
# date:         WED 8 dec 2021 18:26:13 CET
# author:       Abhijeet Singh

###################################
#
version = ': version (1.0.1)'
program = 'filter_by_length.py'

###################################
# imports and check
import sys
import argparse
import subprocess
###################################
try:
    from Bio import SeqIO
except  ImportError:
    print("Biopython missing, attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.78"])
    from Bio import SeqIO
####
try:
    import plotext as plt
except  ImportError:
    print("Plotext missing, attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "plotext"])
    import plotext as plt

###################################
# Help
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description="Filter multifasta sequences by sequence length")
parser.add_argument("-l", "--length_min", dest='length_min', required=True, type=int, help='minimum length of sequence to retain')
parser.add_argument("-m", "--length_max", dest='length_max', required=False, type=int, help='max length of sequence to retain')
parser.add_argument("-i", "--input", dest='input', required=True, type=str, help='input fasta file')
parser.add_argument("-o", "--output", dest='output', required=False, type=str, help='output fasta file')
parser.add_argument("-v", "--verbose", dest='verbose', metavar='Y/y or N/n', default='Y', type=str, help='print progress to the terminal (default:verbose)')
parser.add_argument('-V', '--version', action='version', version='%(prog)s'+ str(version))
args = parser.parse_args()

###################################
# min length
input_length_min = args.length_min
min_length = int(float(input_length_min))

# input file
input_file = args.input
input = open(input_file, 'r')

# output file
if args.output is not None:
    output_file = args.output
    out = open(output_file, 'w')
else:
    pass

# verbosity
if args.verbose is not None:
    verbosity=args.verbose
else:
    verbosity=args.verbose("Y")#default
verbosity=verbosity.upper()    

###################################
# max length
if args.length_max is not None:
    input_length_max = args.length_max
    max_length = int(float(input_length_max))
else:
    max_length = int() 
    
###################################
# main 
length_list = []
filter_list = []
#
for sequence in SeqIO.parse(input, 'fasta'):
    #
    length = len(sequence.seq)
    length_list.append(length)
    min_length_def = min(length_list)
    max_length_def = max(length_list)
    data_size = len(length_list)
    if args.length_max is None:
        #
        max_length = max(length_list)
        input_length_max =  max_length               
    #
    if args.output is not None:
        #
        if  min_length <= len(sequence.seq) and len(sequence.seq) <= max_length:
            SeqIO.write(sequence, out, 'fasta')
            filter_list.append(sequence.id)
            if verbosity == "Y":
                print(">" + sequence.id + "\t" + str(len(sequence.seq)) + "bp" + "\n" + sequence.seq)
    #
    else:
        #
        if  min_length <= len(sequence.seq) and len(sequence.seq) <= max_length:
            filter_list.append(sequence.id)
            #
            if verbosity == "Y":
                #
                
                print(">" + sequence.id + "\t" + str(len(sequence.seq)) + "bp" + "\n" + sequence.seq)
                #print(len(filter_list))

###################################
# count number of filtered sequences
filter_data_size = len(filter_list)

###################################
print('-' * 45)
print('[Number of total sequences]:\t  ' + str(data_size))
print('[Minimum sequence length]:\t  ' + str(min_length_def) + ' bp')
print('[Maximum sequence length]:\t  ' + str(max_length_def) + ' bp')
print('-' * 45)
print('[Number of filtered sequences]:\t  ' + str(filter_data_size))
print('[Minimum sequence length (user)]: ' + str(min_length) + ' bp')
print('[Maximum sequence length (user)]: ' + str(input_length_max) + ' bp')
print('-' * 45)

###################################
# plot in terminal
import plotext as plt
plt.scatter(length_list)
plt.plotsize(40, 20)
plt.show()

###################################
# input file close
input.close()
# output file close
if args.output is not None:
    out.close()

###################################
# End of script
