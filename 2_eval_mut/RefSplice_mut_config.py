import torch
from models.model_utils import (
    get_available_gpus,
)
from models.delta_pretrain import MainModel
from functools import partial
import os
L = 64
# convolution window size in residual units
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
EL = 30000
withcons = False
optimizer = partial(torch.optim.Adam, lr=1e-3)
model = MainModel(L, W, AR, 0.3, EL=EL, optimizer=optimizer)

path = "experiments/2_eval_mut/"


class config:
    log_path = path + "train.log"
    save_path = path
    seed = 321
    EL = EL
    batch_size = 14 * get_available_gpus()
    withcons = withcons

    # for training, can be ignored in test process, except the path to test data
    testjsonfile = [
    ]

    mut_data = [
        "data/vexseq/data.json",
        "data/mfass/data.json",
    ]
    model = model
    epoch_num = 3
    num_workers = 5  # number of dataloader workers

    # for evaluation , can be ignored in training process
    model_path = [os.path.join("RefSplice_models/", x)
                  for x in os.listdir("RefSplice_models/")]
    print_summary = True
    is_train = False
    save_files = True
