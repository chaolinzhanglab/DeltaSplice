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
Fapath = "tempfafile/"
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
IN_MAP = np.eye(6)[:, 1:]
SeqTable = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4, "-": 5}
# Species = ["bosTau9", "susScr11", "mm10", "rheMac10", "rn6", "panTro5"]
context = 10000


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


class MutGenerator():
    def __init__(self,  EL, SavePath, species):
        with open(os.path.join(SavePath, "{}.json".format(species)), "r") as f:
            self.data = json.load(f)
        self.fafiles = {}
        self.EL = EL
        self.s = species

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d = self.data[index]
        hg19_seq = d["hg19_seq"]
        hg19_strand = d["hg19_strand"]
        end_seq = d["{}_seq".format(self.s)]
        end_strand = d["{}_strand".format(self.s)]
        hg19_label = d["hg19_label"]
        end_label = d["{}_label".format(self.s)]
        cidx = len(end_seq)//2
        if hg19_strand == "-":
            hg19_seq = "".join([repdict[_] for _ in hg19_seq.upper()][::-1])

        if end_strand == "-":
            end_seq = "".join([repdict[_] for _ in end_seq.upper()][::-1])
        distance = []

        hg19_seq = IN_MAP[np.array([SeqTable[_] for _ in hg19_seq.upper()])]
        end_seq = IN_MAP[np.array([SeqTable[_] for _ in end_seq.upper()])]
        transition_seq = np.array([hg19_seq, end_seq])

        return {"hg19_seq": hg19_seq, "transition_seq": [hg19_seq, end_seq], "end_seq": end_seq, "hg19_label": [float(_) for _ in hg19_label], "end_label": [float(_) for _ in end_label], "distance": distance}


paths = ('spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(os.path.join('baselines/spliceai_models', x))
          for x in paths]
logger.info("finish loading model")


for species in ["susScr11", "mm10", "rheMac10", "rn6", "panTro5", "bosTau9"]:
    hg19_pred = []
    hg19_gt = []
    end_pred = []
    end_gt = []
    if os.path.exists(os.path.join(sys.argv[2], "spliceai_{}_pred".format(species))):
        continue
    with open(os.path.join(sys.argv[2], "spliceai_{}_pred".format(species)), "w") as f:
        f.write("start")
    dataset = MutGenerator(context, sys.argv[1], species)
    logger.info("handling {}, size is {}".format(species, len(dataset)))
    for idx, d in enumerate(dataset):
        if idx % 500 == 0:
            logger.info(idx)
        inp = np.array(d["transition_seq"])[:, :, :4]

        pred = np.mean([models[m].predict(inp) for m in range(5)], axis=0)

        assert pred.shape[1] == 25001
        assert len(pred.shape) == 3
        didx = pred.shape[1]//2
        pred = pred[:, didx]

        hg19_label = d["hg19_label"]
        end_label = d["end_label"]
        if hg19_label[1] > 0:
            pred = pred[:, 1]
        else:
            pred = pred[:, 2]
        hg19_pred.append(pred[0])
        end_pred.append(pred[-1])
        hg19_gt.append(sum(hg19_label))
        end_gt.append(sum(end_label))
    with open(os.path.join(sys.argv[2], "spliceai_{}_pred".format(species)), "w") as f:
        for a, b, c, d in zip(hg19_pred, hg19_gt, end_pred, end_gt):
            f.writelines("{}\t{}\t{}\t{}\n".format(a, b, c, d))
            f.flush()
