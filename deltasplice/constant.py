import numpy as np
import torch
from deltasplice.models.delta_pretrain import MainModel
from functools import partial
default_model_paths=[f"pretrained_models/DeltaSplice_models/model.ckpt-{i}" for i in range(5)]
default_model_human_paths=[f"pretrained_models/DeltaSplice_human/model.ckpt-{i}" for i in range(5)]

repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
SeqTable = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4, "-": 5}
BASESET = set(["A", "C", "G", "T"])
IN_MAP = np.eye(6)[:, 1:]
CL = 5000
EL = 30000
Fapath="fafiles/"
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


W = [
    11,
    11,
    11,
    11,
    19,
    19,
    19,
    19,
    25,
    25,
    25,
    25,
    33,
    33,
    33,
    33,
    43,
    43,
    85,
    85,
    85,
    85,
    85,
    85,
]
AR = [
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    8,
    8,
    8,
    8,
    16,
    16,
    16,
    16,
    16,
    16,
    32,
    32,
]

optimizer = partial(torch.optim.Adam, lr=2e-3)
model = MainModel(64, W, AR, 0.3, EL=EL, optimizer=optimizer, scheduler=partial(torch.optim.lr_scheduler.StepLR, step_size=1, gamma=0.5))
