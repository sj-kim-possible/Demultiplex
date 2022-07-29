#!/usr/bin/env python
# Author: SJ Kim skim6@uoregon.edu
# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module
'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''
__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning
DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def validate_base_seq(seq: str, RNAflag=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def convert_phred(letter: str) -> int:
    '''Converts single-character quality scores into an int phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """Takes a string representing phred scores and calculates the average quality score of the whole string."""
    sum = 0
    for score in phred_score:
        sum += convert_phred(score)
    return sum/len(phred_score)


def gc_content(dna: str) -> float:
    '''This fx takes a string and calculates the GC content of a DNA sequence.'''
    assert validate_base_seq(dna), "String contains invalid characters."
    dna.isalpha()
    dna=dna.upper()
    gcCount = dna.count("G") + dna.count("C")
    return gcCount/len(dna)

def oneline_fasta(originalFaFilename, intermediaryFaFilename : str):
    '''This function takes a .fa file and removes newline breaks in the protein sequence and
    creates an intermediary .fa file with 2 lines per record: the header and the sequence'''
    with open(originalFaFile, "r") as ogFile:
        with open(intermediaryFaFilename, "w") as noNewLineFile:
            for i, line in enumerate(ogFile):
                header = ""
                record = ""
                if line.startswith(">"):
                    header = line
                    if i == 0:
                        noNewLineFile.write(f"{header}")
                    else:
                        noNewLineFile.write(f"\n{header}")
                else:
                    line = line.strip("\n")
                    record += line
                    noNewLineFile.write(f"{record}")

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
print("Passed DNA and RNA tests")
assert gc_content("GCGCGC") == 1
assert gc_content("AATTATA") == 0
assert gc_content("GCATGCAT") == 0.5
print("correctly calculated GC content")
assert convert_phred("A") == 32, "Phred score incorrect"
assert convert_phred("@") == 31, "Phred score incorrect"
assert convert_phred("#") == 2, "Phred score incorrect"
print("Phred scores converted correctly")