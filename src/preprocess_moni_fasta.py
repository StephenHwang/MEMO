#!/usr/bin/env python3
#
# Run:
#  python preprocess_moni_fasta.py \
#    < /home/stephen/Documents/projects/langmead_lab/omem/data/bacteria_pangeome_fasta/ecoli.fa \
#    >  ecoli_sterilize_w_rc.fa \

from Bio import SeqIO
import argparse
import sys
import textwrap
import warnings

def complement_seq(seq):
    ''' Return the reverse complement of a nucleotide seq. '''
    complement_dict = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G',
                  'N': 'N'
                  }
    return ''.join([complement_dict[base] for base in seq])

def read_records(path):
    ''' Read fasta record. '''
    records = list(SeqIO.parse(path, "fasta"))
    header_list, seq_list = [], []
    for rec in records:
        header_list.append(rec.id)
        seq_list.append(str(rec.seq).upper())
    return header_list, seq_list

def print_seq(header_list, seq_list, reverse=False, complement=False):
    ''' Print converted fasta input. '''
    for header, seq in zip(header_list, seq_list):
        if reverse and complement:
            # print('reverse and complement')
            header = header + '_reverse_complement'
            seq = complement_seq(seq[::-1])
        elif reverse and not complement:
            # print('reverse and no complement')
            header = header + '_reverse'
            seq = seq[::-1]
        elif not reverse and complement:
            # print('no reverse and complement')
            header = header + '_complement'
            seq = complement_seq(seq)
        elif not reverse and not complement:
            # print('no reverse and no complement')
            pass
        else:
            raise Exception('Impossible combination of flags')
        print('>' + header)
        print(textwrap.fill(seq, width=80))


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Reads fasta file from stdin. Output sterilzized sequence with optional rc.")
    parser.add_argument('-c', '--complement', action="store_true", default=False, dest='complement', help='Print reverse of sequence.')
    parser.add_argument('-r', '--reverse', action="store_true", default=False, dest='reverse', help='Print complement of sequence.')
    args = parser.parse_args()
    return args

def main(args):
    ''' Sterize fasta file.

    Input:
        - stdin: fasta
    Output:
        - stdout: sterilized fasta file of records and, optionally, their
                  reverse complement
    '''
    path = sys.stdin
    header_list, seq_list = read_records(path)
    print_seq(header_list, seq_list, reverse=args.reverse, complement=args.complement)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
