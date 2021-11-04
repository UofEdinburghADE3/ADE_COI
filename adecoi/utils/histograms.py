import csv
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns


def make_read_length_histogram_svg(read_lengths,min_len,max_len,outfile):
    new_rc_params = {'text.usetex': False,
    "svg.fonttype": 'none',
    'figure.figsize':[15,5],
    'font.family': 'Helvetica',
    'font.weight' : 'light',
    }

    mpl.rcParams.update(new_rc_params)

    fig,ax= plt.subplots(figsize=(15,5),facecolor='w',frameon=False)
    [ax.spines[loc].set_visible(False) for loc in ['top','right']]
    sns.histplot(data=read_lengths, x="Read length",bins=50)
    plt.axvline(x=min_len, color='dimgrey', linestyle='--')
    plt.axvline(x=max_len, color='dimgrey', linestyle='--')

    plt.xlabel("Read length")
    plt.ylabel("Count")
    plt.savefig(outfile)
