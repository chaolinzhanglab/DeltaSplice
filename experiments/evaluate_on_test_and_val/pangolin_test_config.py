import torch
from models.model_utils import (
    get_available_gpus,
)
from models.pangolin import MainModel
from functools import partial
import os
L = 32
# convolution window size in residual units
W = [11, 11, 11, 11, 11, 11, 11, 11,
     21, 21, 21, 21, 41, 41, 41, 41]
# atrous rate in residual units
AR = [1, 1, 1, 1, 4, 4, 4, 4,
      10, 10, 10, 10, 25, 25, 25, 25]
EL = 30000
withcons = False
optimizer = partial(torch.optim.Adam, lr=1e-3)
model = MainModel(L, W, AR, 0.3, EL=EL, optimizer=optimizer)

path = "experiments/1_evaluate_on_test_and_val/"


class config:
    log_path = path + "train.log"
    save_path = path
    seed = 321
    EL = EL
    batch_size = 14 * get_available_gpus()
    withcons = withcons

    # for training, can be ignored in test process, except the path to test data
    testjsonfile = [
        "data/Npy/valid/data.json",
        "data/Npy/test/data.json",
    ]

    mut_data = None
    model = model
    epoch_num = 3
    num_workers = 5  # number of dataloader workers

    # for evaluation , can be ignored in training process
    model_path = [os.path.join("pangolin_human_models/", x)
                  for x in os.listdir("pangolin_human_models/")]
    print_summary = True
    is_train = False
    save_files = True
