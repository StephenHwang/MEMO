#!/usr/bin/env python3
#
# Run:
#  python preprocess_moni_fasta.py \
#    < /home/stephen/Documents/projects/langmead_lab/omem/data/bacteria_pangeome_fasta/ecoli.fa \
#    >  ecoli_sterilize_w_rc.fa \

from Bio import SeqIO
import sys


def reverse_complement(seq):
    ''' Return the reverse complement of a nucleotide seq. '''
    complement = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G',
                  'N': 'N'
                  }
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

def print_sterilize_seq(seqs, headers, rc=False):
    ''' Print sterilzized fasta input. '''
    num_seqs = len(seqs)
    for header, seq in zip(headers, seqs):
        if rc:
            print('>' + header + '_rc')
            print(reverse_complement(seq))
        else:
            print('>' + header)
            print(seq)


def main():
    ''' Sterize and concatenate fasta file with '.' spacer for MONI input.

    Input:
        - stdin, fasta
    Output:
        - stdout: sterilized fasta file of records and their reverse complement
    '''
    path = sys.stdin
    headers, seqs, seq_lens = read_records(path)
    print_sterilize_seq(seqs, headers)
    print_sterilize_seq(seqs, headers, rc=True)


if __name__ == "__main__":
    main()
