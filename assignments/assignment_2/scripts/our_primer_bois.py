##########################################
# Part 0: Catch sequence between primers #
# Sat, 12:06, SPB, Oct-24-2019           #
# @author: M.V., ITMO, M4235             #
##########################################

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
from collections import namedtuple


### This is reimplementation of list var type, because there was continuous problem with one sequence
### which returned an empty list resulting in Index error. This implementation skips this error, and
### continues the loop.

class RecordList(list):
    def __getitem__(self, index):
        return super(RecordList, self).__getitem__(index)
        
        
def load_primers(file):
  primers = list()
  with open(file, 'r') as f:
    for item in f:
      if not item.startswith('>'):
        primers.append(item.strip('\n'))
        
  return primers
  
  
#modified dictionary with IUPAC codes for fuzzy search
IUPAC = {}
IUPAC['A'] = "A"
IUPAC['C'] = "C"
IUPAC['G'] = "G"
IUPAC['T'] = "T"
IUPAC['M'] = "[A|C]"
IUPAC['R'] = "[A|G]"
IUPAC['W'] = "[A|T]"
IUPAC['S'] = "[C|G]"
IUPAC['Y'] = "[C|T]"
IUPAC['K'] = "[G|T]"
IUPAC['V'] = "[A|C|G]"
IUPAC['H'] = "[A|C|T]"
IUPAC['D'] = "[A|G|T]"
IUPAC['B'] = "[C|G|T]"
IUPAC['X'] = "[G|A|T|C]"
IUPAC['N'] = "[G|A|T|C]"


def make_primers(primers):
    fwd_list = list(primers[0])
    fwd = *map(IUPAC.get, fwd_list),
    fwd_primer_1 = Seq(''.join(map(str, fwd)), generic_dna)
    fwd_primer = fwd_primer_1.transcribe()
    
    revs = (Seq(primers[1])).reverse_complement()
    rev_list = list(revs)
    rev = *map(IUPAC.get, rev_list),
    rev_primer_1 = Seq(''.join(map(str, rev)), generic_dna)
    rev_primer = rev_primer_1.transcribe()

    return fwd_primer, rev_primer 
        

# the real G, catchin' them primer bois in the hood

# Additionally it will also clean records with incomplete annotation

def primer_bois_matchmaking(fwd_primer, rev_primer, fasta_database_file, output_file):

    records = RecordList(SeqIO.parse(fasta_database_file, "fasta"))
    filter = 'Eukaryota|unidentified|Unknown|uncultured|Uncultured|Unclassified|sp.'
    
    for i in range(len(records)):
        condition_0 = bool(re.search(r"{}".format(filter), str(records[i].description)))
        if not condition_0 == True:
            condition_1 = bool(re.search(r"{}".format(fwd_primer), str(records[i].seq)))
            condition_2 = bool(re.search(r"{}".format(rev_primer), str(records[i].seq)))
            if condition_1 and condition_2 == True:
                x_seq = re.findall(r"(?<={})([^']*)(?={})".format(fwd_primer, rev_primer), str(records[i].seq))
                y_id = records[i].description
                if (len(x_seq) > 0) == True:
                    with open(output_file, 'a') as fa:
                        fa.write(">{}\n{}\n".format(y_id, x_seq[0].strip('\n')))

                    
            
def main():
  primers = load_primers(file = '/home/is5/assignment_2/primers.fna')
  fwd_primer, rev_primer = make_primers(primers)
  fwd = str(fwd_primer)
  rev = str(rev_primer)
  
  fasta_database_file='/home/is5/assignment_2/silva.fna.fasta'
  output_file = '/home/is5/assignment_2/clean_silva.fa'
  
  primer_bois_matchmaking(fwd, rev, fasta_database_file, output_file)
  
    
main()
