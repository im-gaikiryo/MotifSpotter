#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, csv, os, sys, datetime, subprocess, platform
from argparse import (ArgumentParser)

try:
    import regex
except ImportError:
    print("The package 'regex' is necesary to run.\nPlease try running 'pip install regex' or see help with '-h' or '--help' option.")
    exit(1)

try:
    from pyexcelerate import Workbook, Style, Font, Alignment, Color
except ImportError:
    print("The package 'pyexcelerate' is necesary to run.\nPlease try running 'pip install pyexcelerate' or see help with '-h' or '--help' option.")
    exit(1)

#=======================================================VERSION=======================================================#
_version = '1.1, Jul-09-2024'
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
* PyExcelerate 0.12.0
* regex 2024.5.15  
* macOS 12.4 Monterey
  Testing for Windows and Linux environment is needed.

USAGE:
- - - 
python3 motifspotter -f $INPUT_FASTA_FILE -m $MORIF_SEQUENCE -e $ALLOWED_MISMATCH -t $MOTIF_TYPE -o $OUTPUT_FORMAT 

For further information and detailed decription of each parameter, please run `python3 motifspotter -h`
=======================================================================================================================
'''

#=======================================================CONSTANT======================================================#
now = datetime.datetime.now()
pf = platform.system()
output_prefix='Matched_Seq_' + now.strftime('%Y%m%d_%H%M%S')
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
# Create translate table for translate() using `str.maketrans()`
#===================================================END OF CONSTANT==================================================#

#===================================================CORE FUNCTIONS===================================================#
def parse_args():
    """Parse the input arguments, use '-h' for help"""
    commands = ArgumentParser(prog='MotifSpotter', description='Regular expressions supported motif searching scripts.')
    commands.add_argument('-f', metavar='FILE', type=str, required=True, help='Source FASTA file for processing.')
    commands.add_argument('-m', metavar='MOTIF', type=str, required=True, help='Target motif sequence. IUPAC ambiguity codes and regular expression is supported.')
    commands.add_argument('-e', metavar='INTEGER', type=int, required=True, help='Maximum allowed ambiguity. Number of maximum allowed mismatch (include subsittution, insertion and deletion) should be specified.')
    commands.add_argument('-t', type=str, required=True, choices=['dna','rna','amino'], help='Type of motif sequence.') 
    commands.add_argument('-o', type=str, choices=['csv','excel'], required=True, help='Desired output format.')
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

motif = "(%s)" % args.m + "{e<=%d}" % args.e
# Compile motif

result = [['Results of spotting motif {}'.format(motif)], ['Identifier', 'Position', 'Matched Sequence', 'Fuzzy Count (substitution, insertion, deletion)', 'Relative change position (substitution, insertion, deletion)']]
# Define a `result` variable with header to store the result.

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
                    
                    rel_span = list(_matched.fuzzy_changes)
                    '''Adding asterisks above and below the matched sequences to indicate location of substitution and insertion'''

                    def get_rel_pos(x):
                        return x - _matched.span()[0]
                    # Define a function to calculate the relative position of the mismatched base in each found sequence
                    # Substracting the mismatched base's absolute position from the first base's absolute position of the matched seuqnce.

                    for i in range(len(rel_span)):
                        if not(len(rel_span[i])==0):
                            rel_span[i] = list(map(get_rel_pos, rel_span[i]))
                    # Using map() function to get relative position of mismatched base
                    

                    if rel_span[0]:
                        max_substitution_position = max(rel_span[0])
                        _substitution_mark = [' '] * (max_substitution_position + 1)
                        # If the first section of the rel_span(representing the relative position of base substitution) isn't null
                        # Assign the largest value in rel_span[0] to `max_substitution_position`
                        # Create a blank list `_substitution_mark` and specify the list size to the `max_substitution_position + 1` 
                        # In this way, the whole `rel_span[0]` can be fitted whithin the list

                        for pos in rel_span[0]:
                            _substitution_mark[pos] = '*'
                        substitution_mark = ''.join(_substitution_mark)
                        # As for the integer `pos` described in the `rel_span[0]`, assign `*` to the nth of the `_substitution_mark` list.
                        # Join the looped result together as string and define as `substitution_mark`

                    else:
                        substitution_mark = ''
                        # If the first section of the rel_span is null, return null to `subsitution_mark`

                    if rel_span[1]:
                        max_insertion_position = max(rel_span[1])
                        _insertion_mark = [' '] * (max_insertion_position + 1)

                        for pos in rel_span[1]:
                            _insertion_mark[pos] = '^'
                        insertion_mark = ''.join(_insertion_mark)
                    else:
                        insertion_mark = ''
                    # Define `insertion_mark` as the location of `insertion`

                    _matched_seq_with_asterisk = substitution_mark + '\n' + _matched.group() + '\n' + insertion_mark


                    _result = [_head, _matched.span(), _matched_seq_with_asterisk, _matched.fuzzy_counts, rel_span]
                    result.append(_result)
                # Search for motif matching sequences in `_seq` and extract as iterator.
                # Assign the identifier, the position of the matching sequence, and the sequence itself to the `_result` variable as a list.
                
        
def csv_write(result):
    """Receive the result and export to specified output file"""

    global output_filename
    output_filename = output_prefix + ".csv"
    with open(output_filename, mode='w') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        writer.writerows(result)

def xlsx_write(_result):
    """Write the result and apply style on it, export as Excel file"""

    global output_filename  
    max = len(result)
    output_filename = output_prefix + ".xlsx"
    # Declare necessary variables and constants

    wb=Workbook()
    ws = wb.new_sheet("Results", data=_result)
    # Write the result

    ws.range("A1", "E2").style.fill.background = Color(189,192,191)
    ws.range("A1", "E2").style.font = Font(bold = True, size = 14)
    ws.range("A1", "E2").style.alignment = Alignment(horizontal= "center", vertical= "center", wrap_text = True)
    ws.range("A1", "E2").style.borders.bottom.color = Color(0, 0, 0)
    ws.range("A2", "E2").style.borders.right.color = Color(0, 0, 0)
    # Configure the style for header

    ws.range("A3", "A%d" % max).style.fill.background = Color(189,192,191)
    ws.range("A3", "A%d" % max).style.font = Font(bold = True, size = 12)
    ws.range("A3", "A%d" % max).style.alignment = Alignment(horizontal= "center", vertical= "center", wrap_text = True)
    ws.range("A3", "A%d" % max).style.borders.right.style = '_'
    ws.range("A3", "A%d" % max).style.borders.bottom.style = '_'
    ws.set_col_style(1, Style(size = -1))
    ws.set_col_style(2, Style(alignment=Alignment(wrap_text=True, horizontal='center', vertical='center'), font=Font(size = 12), size = -1))
    if pf == 'Windows':
        ws.set_col_style(3, Style(alignment=Alignment(wrap_text=True, horizontal='left', vertical='center'), font=Font(family="Consolas", size = 12), size = -1))
    elif pf == 'Darwin':
        ws.set_col_style(3, Style(alignment=Alignment(wrap_text=True, horizontal='left', vertical='center'), font=Font(family="Monaco", size = 12), size = -1))
    elif pf == 'Linux':
        ws.set_col_style(3, Style(alignment=Alignment(wrap_text=True, horizontal='left', vertical='center'), font=Font(family="DejaVu Sans Mono", size = 12), size = -1))
    # Change font type base on 
    ws.set_col_style(4, Style(alignment=Alignment(wrap_text=True, horizontal='center', vertical='center'), font=Font(size = 12), size = -1))
    ws.set_col_style(5, Style(alignment=Alignment(wrap_text=True, horizontal='center', vertical='center'), font=Font(size = 12), size = -1))
    # Configure the style of the remaning part of the sheet

    wb.save(output_filename)


#================================================END OF CORE FUNCTIONS===============================================#

#========================================================MAIN========================================================#
def main():
    try:
        pre_parser(args.f)
        search(motif, args.t)
        if args.o == "csv":
            csv_write(result)
        elif args.o == "excel":
            xlsx_write(result)
        else:
            assert False, f"Invalid argumant {args.o}"

        if os.name == 'nt':
            subprocess.Popen(['start', output_filename], shell=True)
        else:
            subprocess.run(['open', output_filename], check=True)
        # Excute command based on opering system
        
        os.remove('./_tmp')
        return 0

    except OSError:
        print >>sys.stderr, "Execution failed:"

if __name__ == '__main__':
    sys.exit(main())
#====================================================END OF MAIN=====================================================#