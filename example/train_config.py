import torch
from models.model_utils import (
    get_available_gpus,
)
from models.delta_pretrain import MainModel
from functools import partial

EL = 10000
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

withcons = True
optimizer = partial(torch.optim.Adam, lr=1e-3)
model = MainModel(64, W, AR, 0.3, EL=EL, optimizer=optimizer, scheduler=partial(torch.optim.lr_scheduler.StepLR, step_size=1, gamma=0.5), withcons=withcons)

path = "tasks/Pretrain_withcons/"

# must provide a config


class config:
    log_path = path + "train.log"  # the path to save the training log
    save_path = path  # the path to save models & summary during training
    seed = 321
    EL = EL  # the length of extra length. For loading data
    batch_size = 16 * get_available_gpus()
    withcons = withcons  # whether use the conservation score. True is recommended

    # for training, can be ignored in test process, except the path to test data
    trainjsonfile = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/train/data.json"
    validjsonfile = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/valid/human.json"
    testjsonfile = None  # data_path+"dataset_test.hdf5"
    mut_data = None
    model = model
    epoch_num = 3
    num_workers = 5  # number of dataloader workers

    # for evaluation , can be ignored in training process
    model_path = "tasks/test/model/model.ckpt-19"
    print_summary = True
    is_train = True
