#!/usr/bin/env python3
#
# Run:
#  python preprocess_moni_fasta.py \
#    < /home/stephen/Documents/projects/langmead_lab/omem/data/bacteria_pangeome_fasta/ecoli.fa \
#    >  ecoli_sterilize_w_rc.txt \
#    2> ecoli_sterilize_w_rc.bed

from Bio import SeqIO
import sys


def reverse_complement(seq):
    ''' Return the reverse complement of a nucleotide seq. '''
    complement = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G'}
    return ''.join([complement[base] for base in seq.upper()[::-1]])

def read_records(path):
    ''' Read fasta record. '''
    records = list(SeqIO.parse(path, "fasta"))
    headers, seqs, seq_lens = [], [], []
    for rec in records:
        headers.append(rec.id)
        seq = str(rec.seq).upper()
        seqs.append(seq)
        seq_lens.append(len(seq))
    return headers, seqs, seq_lens

def sterilize_seq(seqs):
    ''' Sterilze and concatenate seq and reverse complement of seq. '''
    num_seqs = len(seqs)
    for seq in seqs[:num_seqs]:
        seqs.append(reverse_complement(seq))
    return '.'.join(seqs)

def print_headers_as_bed(headers, seq_lens):
    ''' Print bed positions to stderr. '''
    cumlen = 0
    for header, slen in zip(headers + [header + '_rc' for header in headers], seq_lens * 2):
        print(header + '\t' + str(cumlen), end='\t', file=sys.stderr)
        cumlen += slen
        print(cumlen, file=sys.stderr)
        cumlen += 1
    return 1


def main():
    ''' Sterize and concatenate fasta file with '.' spacer for MONI input.

    Input:
        - stdin, fasta
    Output:
        - stdout: sterilized, concatenated seq with '.' spacer
        - stdout: bed file of indexes of each record
    '''
    path = sys.stdin
    headers, seqs, seq_lens = read_records(path)
    print(sterilize_seq(seqs))
    print_headers_as_bed(headers, seq_lens)

if __name__ == "__main__":
    main()
