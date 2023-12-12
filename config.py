import numpy as np
import os

# the custom path of data files.

Fapath = "/home/xuchencheng/Data/hgfile/"  # the directory contains .fa file of reference genomes
AnnoPath = "/temp/xuchencheng/eSplicedata/annotation/"  # the path to save formated expression data
NpyPath = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy"  # the path to save processed train/test/valid data
TableFile = "/home/xuchencheng/code/WorkDir/gene_dataset.tsu.txt"  # the bed file of gene expression. the table file should in 1-base.


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
