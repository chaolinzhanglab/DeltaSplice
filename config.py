import numpy as np
import os

# the custom path of data files.

Fapath = "hgfile/"  # the directory contains .fa file of reference genomes
AnnoPath = "data/annotation/"  # the path to save formated expression data
NpyPath = "data/Npy"  # the path to save processed train/test/valid data
# the bed file of gene expression. the table file should in 1-base.
TableFile = "data/gene_dataset.tsu.txt"


# remain the following unchanged
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
SeqTable = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4, "-": 5}
BASESET = set(["A", "C", "G", "T"])
IN_MAP = np.eye(6)[:, 1:]
CL = 5000
EL = 30000
npydtype = float

Train_Chromes = [
    "chr15",
    "chr17",
    "chr19",
    "chr21",
    "chr4",
    "chr6",
    "chr8",
    "chr10",
    "chr12",
    "chr14",
    "chr16",
    "chr18",
    "chr20",
    "chr22",
    "chrX",
    "chrY",
    "chr2",
]
Valid_Chromes = ["chr11", "chr13"]
Test_Chromes = ["chr3", "chr5", "chr9", "chr7", "chr1"]
