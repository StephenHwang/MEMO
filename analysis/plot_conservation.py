#!/usr/bin/env python
#
# Run: visualize sequence conservation
#  ./plot_conservation.py \
#    -i memo_conservation.txt \
#    -o out.png \
#    -n num_genomes_in_pangenome \
#    -b num_bins

import argparse
import pandas as pd
import numpy as np
from collections import Counter
from plotnine import \
    ggplot, aes, theme, themes, element_blank, element_line, element_text, \
    geom_bar, ggtitle, xlab, ylab, scale_x_continuous, scale_y_continuous, \
    scale_fill_gradient
from plotnine.options import figure_size


# plotnine theme
def theme_tufte_func(base_size=19, base_family="sans", figure_size=None):
    ''' Plotnine plotting theme. '''
    thm = themes.theme_bw(base_size=base_size, base_family = base_family) + \
        theme(
            legend_background = element_blank(),
            legend_key = element_blank(),
            panel_background = element_blank(),
            panel_border = element_blank(),
            strip_background = element_blank(),
            plot_background = element_blank(),
            panel_grid = element_blank(),
            axis_line = element_line(colour = "black", size = 1),
            axis_text_y = element_text(colour = "black")
        )
    if figure_size is not None:
        thm += theme(figure_size=figure_size)
    return thm


def fileReader(path):
    ''' Read file from path line-by-line. '''
    with open(path, 'r') as inFile:
        for line in inFile:
            yield line.strip()

def preprocess_data(path, n_docs, n_bins):
    ''' Preprocess conservation data for plotting. '''
    num_docs_per_pos = list(map(int, list(fileReader(path))))
    positions = len(num_docs_per_pos)

    per_bin_doc_composition_list = []
    bin_space = list(map(int,np.linspace(0, positions, n_bins)))

    for bin_idx, start_end in enumerate(list(zip(bin_space[:-1], bin_space[1:]))):
        bin_start, bin_end = start_end
        doc_count_per_order_in_bin = Counter(num_docs_per_pos[bin_start : bin_end])
        normalized_doc_count_per_order_in_bin = [(order, doc_count_per_order_in_bin[order]/sum(doc_count_per_order_in_bin.values())) for order in range(n_docs + 1) ]
        normalized_doc_count_per_order_in_bin_in_sorted_order_mem_order = sorted(normalized_doc_count_per_order_in_bin, key=lambda x: x[0])   # sorting by order
        per_bin_doc_composition_list.append([bin_idx]+[norm_cnt[1] for norm_cnt in normalized_doc_count_per_order_in_bin_in_sorted_order_mem_order])

    cnames = ['pos'] + list(range(n_docs + 1))
    per_bin_doc_composition_df = pd.DataFrame(per_bin_doc_composition_list, columns=cnames)

    per_bin_doc_composition_df = pd.melt(per_bin_doc_composition_df, id_vars=['pos'], value_vars=cnames[1:])
    per_bin_doc_composition_df.columns = ['bin','No. Genomes','value']
    per_bin_doc_composition_df['No. Genomes'] = pd.Categorical(per_bin_doc_composition_df['No. Genomes'], categories=cnames[:0:-1])
    per_bin_doc_composition_df['No. Genomes'] = per_bin_doc_composition_df['No. Genomes'].astype('float64')
    return per_bin_doc_composition_df[per_bin_doc_composition_df['No. Genomes'] != n_docs]


def plot_sequence_conservation(data, n_docs, n_bins):
    ''' Plot MEM-version of Panagram. '''
    p = (
        ggplot(data, aes(x='bin', y='value', fill='No. Genomes')) +
        geom_bar(stat="identity", width=1) +
        ggtitle("K-mer Conservation") +
        xlab("Genomic bin (n =" + str(n_bins) + ")" ) +
        ylab("Proportion of\nconserved k-mers") +
        scale_y_continuous(
            breaks=np.linspace(0,1,5),
            labels=['0','0.25','0.50','0.75', '1'],
            expand=(0,0),
            limits=(0,1)
        ) +
        scale_x_continuous(expand=(0,0)) +
        scale_fill_gradient(low='#000000', high='#c6dbef',
                            limits=(1,n_docs-1)
                           ) +
        theme_tufte_func(base_size=18, figure_size=[20,4])
    )
    return p

def save_plot(p, out_file, dpi=600):
    ''' Save plot to png. '''
    p.save(filename=out_file, dpi=dpi)


###############################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Extract and query overlap MEMs for k-mer presence/absence.")
    parser.add_argument('-i', '--in_file', dest='in_file', help='in file', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='output file', required=True)
    parser.add_argument('-n', '--ndocs', dest='n_docs', help='total number of genomes in the pangenome', required=True)
    parser.add_argument('-b', '--num_bins', dest='num_bins', help='number of histogram bins', required=True)
    parser.add_argument('-d', '--dpi', dest='dpi', default=600, help='Plot dpi', required=False)
    args = parser.parse_args()
    return args

def main(args):
    ''' Plot sequence conservation. '''
    in_file = args.in_file
    out_file = args.out_file
    n_docs = int(args.n_docs)
    n_bins = int(args.num_bins)
    dpi = int(args.dpi)
    data = preprocess_data(in_file, n_docs, n_bins)
    plot = plot_sequence_conservation(data, n_docs, n_bins)
    save_plot(plot, out_file, dpi)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
