from Bio import SeqIO
from random import randrange
import numpy as np

################################################################################

def mutate_base(base):
    bases = set('atcg')
    bases.remove(base)
    mut_base = list(bases)[randrange(len(bases))]
    return mut_base

def gen_mutation_window(start_idx, end_idx, r):
    ''' Note: Use index 0 '''
    st = set()
    window_size = end_idx - start_idx
    while len(st) < window_size * r:
        st.add(randrange(start_idx, end_idx))
    return st

def mutate_seq(seq, start_idx, end_idx, r):
    seq = seq.lower()
    new_seq = list(seq)
    end_idx = min(len(seq), end_idx)
    st = gen_mutation_window(start_idx, end_idx, r)
    for pos in st:
        new_seq[pos] = mutate_base(seq[pos])
    return ''.join(new_seq)

def mutate_full_chr(seq, r, window_size):
    for pos in list(np.arange(0, len(seq), window_size)):
        start_idx = pos
        end_idx = pos + window_size
        seq = mutate_seq(seq, start_idx, end_idx, r)
    return seq

def check_mutation_rate(seq, seq_mutated):
    cnt = 0
    for x, y in zip(seq, seq_mutated):
        if x != y:
            cnt += 1
    return cnt / len(seq)

def print_seq(seq):
    for pos in list(np.arange(0, len(seq), 80)):
        start_idx = pos
        end_idx = min(len(seq), pos + 80)
        print(seq[start_idx:end_idx].upper())

################################################################################


path = '/home/stephen/Documents/projects/langmead_lab/omem/data/bacteria_pangeome_fasta/ecoli.fa'
divergence_rates = (0.00, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40)
r = divergence_rates[1]
window_size = 10_000

with open(path) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        header = record.id
        seq = str(record.seq).lower()
        seq_mutated = mutate_full_chr(seq, r, window_size)
        print('>' + header + '_' + str(r) + 'r')
        print_seq(seq_mutated)

