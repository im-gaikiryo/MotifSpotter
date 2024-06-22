#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, csv, os, sys, datetime, subprocess
from argparse import (ArgumentParser, FileType)

try:
    import regex
except ImportError:
    print("The package 'regex' is necesary to run.\nPlease try running 'pip install regex' or see help with '-h' or '--help' option.")
    exit(1)

#=======================================================VERSION=======================================================#
_version = '1.0, Jun-19-2024'
#====================================================END OF VERSION===================================================#

docinfo = '''
=======================================================================================================================
MotifSpotter: A regular expressions supported motif searching script
Author: Gaiki Ryo/Kaiqi Liang/Carter Leung
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MIT License

Copyright (c) 2024 Gaiki Ryo/Kaiqi Liang/Carter Leung

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TESTING ENVIRONMENT:
- - - - - - - - - - 
* Python 3.12.2
  Theoretically, this script should works in Python 3.5+ environment.
* macOS 12.4 Monterey
  Testing for Windows and Linux environment is needed.

USAGE:
- - - 
python3 motifspotter -f $INPUT_FASTA_FILE -m $MORIF_SEQUENCE -t $MOTIF_TYPE [-o] [$OUTPUT_FILE]

For further information and detailed decription of each parameter, please run `python3 motifspotter -h`
=======================================================================================================================
'''

#======================================================VARIABLES======================================================#
result = []
now = datetime.datetime.now()
IUPAC_alphabet = re.compile(r'[^a-zA-Z-. \n]+')
IUPAC_dna_table = str.maketrans({
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "[AC]",
    "R": "[AG]",
    "W": "[AT]",
    "S": "[CG]",
    "Y": "[CT]",
    "K": "[GT]",
    "V": "[ACG]",
    "H": "[ACT]",
    "D": "[AGT]",
    "B": "[CGT]",
    "X": "[GATC]",
    "N": "[GATC]",
}) 
IUPAC_rna_table = str.maketrans({
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "[AC]",
    "R": "[AG]",
    "W": "[AU]",
    "S": "[CG]",
    "Y": "[CU]",
    "K": "[GU]",
    "V": "[ACG]",
    "H": "[ACU]",
    "D": "[AGU]",
    "B": "[CGU]",
    "X": "[GAUC]",
    "N": "[GAUC]",
})
IUPAC_amino_table = str.maketrans({
    "A": "A", #Ala
    "C": "C", #Cys
    "D": "D", #Asp
    "E": "E", #Glu
    "F": "F", #Phe
    "G": "G", #Gly
    "H": "H", #His
    "I": "I", #Ile
    "K": "K", #Lys
    "L": "L", #Leu
    "M": "M", #Met
    "N": "N", #Asn
    "P": "P", #Pro
    "Q": "Q", #Gln
    "R": "R", #Arg
    "S": "S", #Ser
    "T": "R", #Thr
    "V": "V", #Val
    "W": "W", #Trp
    "Y": "Y", #Tyr
    "B": "[DN]", #Asx
    "Z": "[EQ]", #Glx
    "J": "[LI]", #Xle
    "X": "[A-Z]", #Xxx
})
#Create translate table for translate() using `str.maketrans()`
#==================================================END OF VARIABLES==================================================#

#===================================================CORE FUNCTIONS===================================================#
def parse_args():
    """Parse the input arguments, use '-h' for help"""
    commands = ArgumentParser(prog='MotifSpotter', description='Regular expressions supported motif searching scripts.')
    commands.add_argument('-f', metavar='FILE', type=str, required=True, help='Source FASTA file for processing.')
    commands.add_argument('-m', metavar='MOTIF', type=str, required=True, help='Target motif sequence. Regular expression is supported.')
    commands.add_argument('-t', type=str, required=True, choices=['dna','rna','amino'], help='Type of motif sequence.') 
    #commands.add_argument('-e', metavar='INTEGER', type=int, required=True,  help='Maximum allowed ambiguity. Number of maximum allowed mismatch should be specified')
    commands.add_argument('-o', metavar='OUTPUT', type=str, required=False, default='Matched_Seq_' + now.strftime('%Y%m%d_%H%M%S') + '.csv', help='Location and name for output file.')
    commands.add_argument('-v', action='version', version='%(prog)s {}'.format(_version))
    return commands.parse_args()
args = parse_args()

def pre_parser(_input_file):
    """Checking if input file is in valid FASTA format and remove blank line."""
    with open(_input_file, 'r') as f, open('_tmp', 'w') as o:
        ff = list(f)
        if not(ff[0].startswith('>')):
               raise IOError("NOT A VALID FASTA FILE")
        # To make sure input file is a valid FASTA file
        for line in ff:
            if line.strip():
                o.write(line)
        # Empty objects are considered false, and string that is not empty will be considered to be true in a boolean.
        # The intention of using `if str.strip()` here is to filter out empty lines and execute further commands only for normal strings.
        # Lines containing only spaces and empty strings are filtered out.
        # The result is written to ./_tmp for further processing.

def search(_motif, _type):
    """Searching for target motif in each sequence and export the result as group."""
    if _type == 'dna':
        regexed_motif = _motif.translate(IUPAC_dna_table)
    elif _type == 'rna':
        regexed_motif = _motif.translate(IUPAC_rna_table)
    elif _type == 'amino':
        regexed_motif = _motif.translate(IUPAC_amino_table)
    else:
        raise IOError('NOT VALID MOTIF TYPE')
    # Translate motif table into corresponded regex form.

    
    with open('./_tmp', 'r') as f:
        for line in f:
            if line.startswith(">"):
                _head =line.strip("\n>")
                _seq = next(f, '').strip('').replace(" ","")
                # Extract and assign identifier and sequences to variable `_head` and `_seq`, respectively.
                
                if IUPAC_alphabet.search(_seq):
                    raise IOError("NOT A VALID FASTA FILE")
                # Make sure _seq do not contains any invalid string
                
                for _matched in regex.finditer(regexed_motif, _seq, regex.BESTMATCH, overlapped=True):
                    _result = [_head, _matched.span(), _matched.group(), _matched.fuzzy_counts]
                    result.append(_result)
                # Search for motif matching sequences in `_seq` and extract as iterator.
                # Assign the identifier, the position of the matching sequence, and the sequence itself to the `_result` variable as a list.
                
        
def write(result, _output_csv):
    """Receive the result and export to specified output file"""
    with open(_output_csv, mode='w') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        writer.writerow(['Results of spotting motif {}'.format(args.m)])
        writer.writerow(['Identifier', 'Position', 'Sequence', 'Fuzzy Count (substitution, insertion, deletion)'])
        writer.writerows(result)
#================================================END OF CORE FUNCTIONS===============================================#

#========================================================MAIN========================================================#
def main():
    try:
        pre_parser(args.f)
        search(args.m, args.t)
        write(result, args.o)
        os.remove('./_tmp')
        subprocess.run(['open', args.o], check=True)
        return 0
    except OSError:
        print >>sys.stderr, "Execution failed:"

if __name__ == '__main__':
    sys.exit(main())
#====================================================END OF MAIN=====================================================#