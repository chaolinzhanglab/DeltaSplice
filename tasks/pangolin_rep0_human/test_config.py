from .config import *
from .config import config as oldconfig


class config(oldconfig):
    # for training, can be ignored in test process, except the path to test data
    testjsonfile = [
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/valid/data.json",
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Npy/test/data.json",
    ]

    # for evaluation , can be ignored in training process
    model_path = path+"models/model.ckpt-2"
    print_summary = True
    save_files = True
    is_train = False
