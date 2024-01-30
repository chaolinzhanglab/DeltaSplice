from .config import *
from .config import config as oldconfig


class config(oldconfig):
    # for training, can be ignored in test process, except the path to test data
    mut_data = [
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/VexSeq/data.json",
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/MFASS/data.json",
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Autism_genome/data.json",
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Autism_exome/data.json",
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Mapsy_vivo/data.json",
        "/temp/xuchencheng/eSplicedata/Species_37_ProcessedData/Mapsy_vitro/data.json",
    ]

    # for evaluation , can be ignored in training process
    model_path = path+"models/model.ckpt-2"
    print_summary = True
    save_files = True
    is_train = False
