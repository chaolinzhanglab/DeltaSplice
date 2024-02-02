from .config import *
from .config import config as oldconfig


class config(oldconfig):
    # for training, can be ignored in test process, except the path to test data
    mut_data = [
        "data/vexseq/data.json",
        "data/mfass/data.json",
        "data/autism_genome/data.json",
        "data/autism_exome/data.json",
    ]

    # for evaluation , can be ignored in training process
    model_path = path+"models/model.ckpt-temp"
    print_summary = True
    save_files = True
    is_train = False
