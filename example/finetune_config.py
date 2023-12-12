import torch
from models.model_utils import (
    get_available_gpus,
)
from models.delta_ft import MainModel
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
model = MainModel(64, W, AR, 0.3, EL=EL, optimizer=optimizer, scheduler=partial(torch.optim.lr_scheduler.StepLR, step_size=1, gamma=0.5), withcons=withcons, numensemble=1)
state_dict = torch.load("tasks/Pretrain_withcons/models/model.ckpt-best")
for name in model.state_dict():
    if name in state_dict:
        model.state_dict()[name].copy_(state_dict[name])
    else:
        print(name)
path = "tasks/Finetune_withcons/"


class config:
    log_path = path + "train.log"
    save_path = path
    seed = 321
    EL = EL
    batch_size = 64 * get_available_gpus()
    withcons = withcons

    # for training, can be ignored in test process, except the path to test data
    trainjsonfile = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/train/human.json"
    validjsonfile = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/valid/human.json"
    testjsonfile = None  # data_path+"dataset_test.hdf5"
    mut_data = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/VexSeq/data.json"
    model = model
    epoch_num = 1
    num_workers = 5  # number of dataloader workers

    # for evaluation , can be ignored in training process
    model_path = "tasks/test/model/model.ckpt-19"
    print_summary = True
    is_train = True
