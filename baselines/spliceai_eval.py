from tensorflow.keras.models import load_model
import numpy as np
import json
from pyfasta import Fasta
import os
import sys
import tensorflow as tf
from loguru import logger

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    # Restrict TensorFlow to only allocate 1GB of memory on the first GPU
    try:
        tf.config.experimental.set_virtual_device_configuration(
            gpus[0],
            [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=5096)])  # Notice here
        logical_gpus = tf.config.experimental.list_logical_devices('GPU')
        print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    except RuntimeError as e:
        # Virtual devices must be set before GPUs have been initialized
        print(e)

input_sequence = 'CGATCTGACGTGGGTGTCATCGCATTATCGATATTGCAT'
# Replace this with your custom sequence

# Fapath = "fafiles/"
Fapath = "py36hgfile/"
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
IN_MAP = np.eye(6)[:, 1:]
SeqTable = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4, "-": 5}
# Species = ["bosTau9", "susScr11", "mm10", "rheMac10", "rn6", "panTro5"]
context = 15000


class DataGenerator():
    def __init__(self,  jsonfile=None):
        with open(jsonfile, "r") as f:
            self.data = json.load(f)
        self.fafiles = {}
        self.EL = 30000

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d = self.data[index]
        species, chrom, start, end, strand, name = d["species"], d["chrom"], int(
            d["start"]), int(d["end"]), d["strand"], d["name"]
        if species not in self.fafiles:
            self.fafiles[species] = Fasta(os.path.join(Fapath, species+".fa"))

        seq = self.fafiles[species][chrom][start:end].upper()
        label = np.zeros([end-start-self.EL, 3])
        for v in d["label"]:
            idx, value = v
            idx = int(idx)
            if idx-self.EL//2 >= 0 and idx < end-start-self.EL//2:
                value = np.array([float(_) for _ in value])
                label[idx-self.EL//2] = value
        if strand == "-":
            seq = [repdict[_] for _ in seq][::-1]
            label = np.copy(label[::-1])
        seq = IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
        label[:, 0] = 1.-label[:, 1:].sum(-1)

        return {"X": seq, "single_pred_psi": label, "species": species, "chrom": chrom, "name": name, "strand": strand, "txstart": start, "txend": end}


paths = ('spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(os.path.join('baselines/spliceai_models', x))
          for x in paths]
logger.info("finish loading model")


dataset = DataGenerator(sys.argv[1])
logger.info("handling size is {}".format(len(dataset)))
save_path = []


def write_splice_site_file_header(fout):
    yt = ["Yt"]
    yp = ["Yp"]
    header = "\t".join(
        ["#species", "chrom", "chromStart", "chromEnd",
            "name", "score", "strand", "type"]
        + yt
        + yp
    )
    fout.write(header + os.linesep)
    fout.flush()


def write_splice_sites(fouts, CHROM, NAME, STRAND, TX_START, TX_END, SPECIES, Yp, Yt):
    output_class_labels = ["Null", "acceptor", "donor"]
    # The three neurons per output correspond to no splicing, splice acceptor (AG)
    # and splice donor (GT) respectively.

    num_row = 0

    if num_row == 0:
        num_row = Yt.shape[0]
    else:
        assert num_row == Yt.shape[0]
        assert num_row == Yp.shape[0]

    for i in range(num_row):
        chrom = CHROM[
            i
        ]  # .decode()  #TODO: there is compatibility issue with more recent version of h5py
        name = NAME[i]  # .decode()
        strand = STRAND[i]  # .decode()
        tx_start = TX_START[i]
        tx_end = TX_END[i]
        species = SPECIES[i]

        for c in [1, 2]:
            fo = fouts[c - 1]

            y_t = np.copy(Yt[i, :, c])
            y_t[np.isnan(y_t)] = 1
            idx = np.nonzero(y_t > 1e-10)[0]  # positions of true splice sites

            for j in idx:
                pos = (tx_start + tx_end)//2

                yt = [str(Yt[i, j, c])]
                yp = [str(Yp[i, j, c])]

                line = "\t".join(
                    map(
                        str,
                        [
                            species,
                            chrom,
                            pos,
                            pos+1,
                            name,
                            0,
                            strand,
                            output_class_labels[c],
                        ]
                        + yt
                        + yp,
                    )
                )
                fo.write(line + os.linesep)
                fo.flush()


Y_true_1 = []
Y_true_2 = []
Y_pred_1 = []
Y_pred_2 = []
save_path.append(open(sys.argv[2]+"_acceptor", "w"))
save_path.append(open(sys.argv[2]+"_donor", "w"))
[write_splice_site_file_header(_) for _ in save_path]
for idx, d in enumerate(dataset):
    if idx % 500 == 0:
        logger.info(idx)
    inp = np.array(d["X"])[None][:, :, :4]

    pred = np.mean([models[m].predict(inp) for m in range(5)], axis=0)

    assert pred.shape[1] == 25000, pred.shape
    pred = pred[:, 10000:15000]
    assert len(pred.shape) == 3

    label = d["single_pred_psi"][None]
    assert pred.shape == label.shape

    pred_sites = np.nonzero(label[:, :, 1:].sum(-1))

    CHROM = [d["chrom"]]
    NAME = [d["name"]]
    STRAND = [d["strand"]]
    TX_START = np.array([d["txstart"]])
    TX_END = np.array([d["txend"]])
    SPECIES = [d["species"]]
    write_splice_sites(save_path, CHROM, NAME, STRAND,
                       TX_START, TX_END, SPECIES, pred, label)

[_.close() for _ in save_path]
