#!/usr/bin/env python3
# Plotly + Dash for order MEM panagram

from dash import Dash, dcc, html, Input, Output
import subprocess
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from collections import Counter



def call_query(k):
    cmd = [
        "/home/stephen/Documents/projects/langmead_lab/omem/src/query.sh",
        "-k", str(k),
        "-n", "4",
        "-o", "/home/stephen/Documents/projects/langmead_lab/omem/data/example_dap",
        "-b", "omem.bed.gz",
        "-r", "NZ_CP015023.1:0-5506800",
        #"-i", # save intermediate
        # "-r", "omem_overlaps.bed",
        "-p"
        ]
    subprocess.check_call(cmd)
    return 1


def update_data(k, n_bins, num_docs):
    call_query(k)

    bed_dir = '/home/stephen/Documents/projects/langmead_lab/omem/data/example_dap/'
    MEM_bed_path = bed_dir + 'omem_' + str(k) + 'mer.bed'
    mem_bed = pd.read_csv(MEM_bed_path, sep='\t', header=None, names=['chrm', 'start', 'end', 'order'])

    positions = list(mem_bed.index)
    num_orders = len(set(mem_bed['order']))
    mem_bed_order_matrix = np.zeros((num_orders+1, len(positions)))

    for order in range(num_orders):
        mem_bed_order_subset = mem_bed[mem_bed['order'] == order+1]
        for start, end in zip(mem_bed_order_subset['start'], mem_bed_order_subset['end']):
            mem_bed_order_matrix[order, start:end] = 1
    mem_bed_order_matrix[order+1,:] = 1
    num_docs_per_pos = np.argmax(mem_bed_order_matrix, axis=0)

    per_bin_doc_composition_list = []
    bin_space = list(map(int,np.linspace(0, len(positions), n_bins)))

    for bin_idx, start_end in enumerate(list(zip(bin_space[:-1], bin_space[1:]))):
        bin_start, bin_end = start_end
        doc_count_per_order_in_bin = Counter(num_docs_per_pos[bin_start : bin_end])
        for doc_no_count in set(range(0,num_docs+1)) - set(doc_count_per_order_in_bin.keys()):
            doc_count_per_order_in_bin[doc_no_count] = 0
        normalized_doc_count_per_order_in_bin = [(order, cnt/sum(doc_count_per_order_in_bin.values())) for order, cnt in doc_count_per_order_in_bin.items()]
        normalized_doc_count_per_order_in_bin_in_sorted_order_mem_order = sorted(normalized_doc_count_per_order_in_bin, key=lambda x: x[0])   # sorting by order
        per_bin_doc_composition_list.append([bin_idx]+[norm_cnt[1] for norm_cnt in normalized_doc_count_per_order_in_bin_in_sorted_order_mem_order])

    cnames = ['pos', '1', '2', '3', '4', '5']
    per_bin_doc_composition_df = pd.DataFrame(per_bin_doc_composition_list, columns=cnames)
    per_bin_doc_composition_df = pd.melt(per_bin_doc_composition_df, id_vars=['pos'], value_vars=cnames[1:])
    per_bin_doc_composition_df.columns = ['bin','Num docs','value']
    per_bin_doc_composition_df['Num docs'] = pd.Categorical(per_bin_doc_composition_df['Num docs'], categories=cnames[:0:-1])

    return per_bin_doc_composition_df

def view(n_bins):
    ''' '''
    colors = ['#fde725', '#5cc863', '#21908d', '#3b518b', '#440154']
    num_docs = 4

    app = Dash(__name__)

    app.layout = html.Div([
        html.H4('Interactive plot with custom data source'),
        dcc.Graph(id="graph"),
        html.P("k:"),
        dcc.Slider(id="slider", min=10, max=20, value=12, step=1),
    ])

    @app.callback(
        Output("graph", "figure"),
        Input("slider", "value"))
    def update_bar_chart(k):
        # run query
        data = update_data(k, n_bins, num_docs)

        fig = go.Figure()
        for doc_idx in range(1, num_docs+2):
            doc_data = data[data['Num docs'] == str(doc_idx)]
            x, y = doc_data['bin'], doc_data['value']
            color = dict(color=colors[doc_idx-1])
            name = str(doc_idx)

            fig.add_trace(go.Bar(x=x, y=y, name = name,
                                 legendgroup="Num docs",
                                 legendgrouptitle_text="Num docs",
                                 marker=color,
                                 marker_line=color))

        fig.update_layout(barmode='stack', bargap=0.0, template="simple_white",
                          title={'text': "Conserved " + str(k) + '-mers',
                                'x':0.5, 'xanchor': 'center',
                                'yanchor': 'bottom'})
        fig.update_xaxes(title_text="Bin (n=" + str(n_bins) +')')
        fig.update_yaxes(title_text="Proportion conserved")

        return fig

    app.run_server(debug=True)


def main():
    ''' Run view. '''
    n_bins = 350*2
    view(n_bins)


if __name__ == "__main__":
    main()
