##########################################
# Part 1: Salvation for silva sequences  #
# Sat, 00:18, SPB, Oct-19-2019           #
# @author: M.V., ITMO, M4235             #
##########################################

# Problem statement:
# 
# 1.) I must parse a database. It has FASTA sequences.
# 2.) I should get a fasta file with headers, and a taxonomy file in the given format

from collections import namedtuple
import string

"""Properly parse fasta file"""

def parse_fasta(stream):
    Fasta = namedtuple('Fasta', 'header sequence')
    header = ''
    sequence = ''

    for line in stream:
        line = line.strip()

        if line.startswith('>'):
            if header:
                yield Fasta(header, sequence)

            header = line[1:].strip()
            sequence = ''
        else:
            sequence += line

    if header:
        yield Fasta(header, sequence)
  

"""Parse each label in the fasta file"""

def parse_label(label):
    new_header,taxonomy = label.split(' ',1)
    return new_header,taxonomy

"""Transcription i.e. convert T to U and remove any whitespaces"""

def rna_strand(strand):
    rna = ""
        # Generate the RNA string
    for items in strand.split():
        for i in items:
            # Replace all occurrences of T with U
            if i == "U":
                rna += "T"
            else:
                rna += i

        return rna
        
"""The real processing function"""

def process_silva(seqs, tax_out, seq_out):
    for label,seq in parse_fasta(seqs):
        new_header,taxonomy = parse_label(label)
        fixed_seq = rna_strand(seq)
        tax_out.write(new_header + '\t' + taxonomy + '\n')
        seq_out.write('>' + new_header + '\n' + fixed_seq + '\n')


def main(seqs,tax_out,seq_out):
    process_silva(seqs, tax_out, seq_out)
    tax_out.close()
    seq_out.close()


if __name__ == '__main__':
    from sys import argv
    USE = "\nUSE: python silvation.py silva_fasta_infile taxonomy_outfile sequence_outfile\n"
    
    if len(argv) != 4:
        print(USE)
    else:
        seqs = open(argv[1])
        tax_out = open(argv[2], 'w')
        seq_out = open(argv[3], 'w')
        main(seqs, tax_out, seq_out)
        
        
