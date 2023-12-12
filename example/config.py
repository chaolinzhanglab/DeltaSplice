import numpy as np
import pyBigWig
import os

# the custom path of data files.

Fapath = "data/Fa"  # the directory contains .fa file of reference genomes
AnnoPath = "data/annotation/"  # the path to save formated expression data
NpyPath = "data/Npy"  # the path to save processed train/test/valid data
WGPath = "data/WGFile"  # the prefix of conservation file
WGFiles = {"hg19": "placentalMammals.phyloP46way.bw", "mm10": "mm10.60way.phyloP60way.bw"}  # {genome: its conservation file in .bw format}
TableFile = "example/gene_dataset.tsu.txt"  # the bed file of gene expression. the table file should in 1-base.


# remain the following unchanged
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
SeqTable = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4, "-": 5}
BASESET = set(["A", "C", "G", "T"])
IN_MAP = np.eye(6)[:, 1:]
CL = 5001
EL = 10000
npydtype = float


def getWGFile(WGPath=WGPath):
    return {key: pyBigWig.open(os.path.join(WGPath, WGFiles[key])) for key in WGFiles}


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
