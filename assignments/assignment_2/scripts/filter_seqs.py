import re
from collections import namedtuple


infile = 'output/taxon_fasta_after_silva_processing.fasta'

patterns = ["Bacteria"]

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
  

for seq_record in parse_fasta(infile):
    for p in patterns:
        if p.lower() in seq_record.description.lower():
            print(seq_record.format("fasta"))
            break
