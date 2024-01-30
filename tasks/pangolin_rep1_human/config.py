import torch
from models.model_utils import (
    get_available_gpus,
)
from models.pangolin import MainModel
from functools import partial
L = 32
# convolution window size in residual units
W = [11, 11, 11, 11, 11, 11, 11, 11,
     21, 21, 21, 21, 41, 41, 41, 41]
# atrous rate in residual units
AR = [1, 1, 1, 1, 4, 4, 4, 4,
      10, 10, 10, 10, 25, 25, 25, 25]
EL = 10000
withcons = False
optimizer = partial(torch.optim.Adam, lr=1e-3)
model = MainModel(L, W, AR, 0.3, EL=EL, optimizer=optimizer)

path = "tasks/pangolin_rep1_human/"


class config:
    log_path = path + "train.log"
    save_path = path
    seed = 49
    EL = EL
    batch_size = 14 * get_available_gpus()
    withcons = withcons

    # for training, can be ignored in test process, except the path to test data
    trainjsonfile = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/train/human.json"
    validjsonfile = "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/valid/human.json"

    testjsonfile = None  # data_path+"dataset_test.hdf5"
    mut_data = None
    model = model
    epoch_num = 10
    num_workers = 5  # number of dataloader workers

    # for evaluation , can be ignored in training process
    model_path = "tasks/test/model/model.ckpt-19"
    print_summary = True
    is_train = True
