from tensorflow.keras.models import load_model
from pkg_resources import resource_filename
from spliceai_utils import one_hot_encode
import numpy as np
import json
from pyfasta import Fasta
import os
import sys
import tensorflow as tf
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

# Fapath="fafiles/"
Fapath = "py36hgfile/"
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
context = 10000


class MutGenerator():
    def __init__(self,  EL, jsonfile):
        with open(jsonfile, "r") as f:
            self.data = json.load(f)
        self.fafiles = {}
        self.EL = EL

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d = self.data[index]
        if "mutpos" in d:
            species, chrom, start, end, strand, mutpos, oriv, mutv = d["species"], d["chrom"], int(
                d["start"]), int(d["end"]), d["strand"], int(d["mutpos"]), d["oriv"].upper(), d["mutv"].upper()
            if species not in self.fafiles:
                self.fafiles[species] = Fasta(
                    os.path.join(Fapath, species+".fa"))

            seq = self.fafiles[species][chrom][start:end].upper()
            mutseq = seq[:mutpos]+mutv+seq[mutpos+1:]
            assert seq[mutpos] == oriv
        else:
            chrom, start, end, strand, seq, mutseq = d["chrom"], int(d["start"]), int(
                d["end"]), d["strand"], d["seq"].upper(), d["mutseq"].upper()
            # assert strand=="+"

        label = np.zeros([end-start-self.EL, 3])
        mutlabel = np.zeros([end-start-self.EL, 3])

        for v in d["label"]:
            idx, value = v
            idx = int(idx)
            if idx-self.EL//2 >= 0 and idx < end-start-self.EL//2:
                value = np.array([float(_) for _ in value])
                label[idx-self.EL//2] = value

        for v in d["mutlabel"]:
            idx, value = v
            idx = int(idx)
            assert idx >= self.EL//2 and idx < end-start-self.EL//2
            value = np.array([float(_) for _ in value])
            mutlabel[idx-self.EL//2] = value

        if strand == "-":
            seq = "".join([repdict[_] for _ in seq][::-1])
            mutseq = "".join([repdict[_] for _ in mutseq][::-1])
            label = np.copy(label[::-1])
            mutlabel = np.copy(mutlabel[::-1])

        seq = one_hot_encode(seq)
        mutseq = one_hot_encode(mutseq)

        label[:, 0] = 1.-label[:, 1:].sum(-1)
        return {"X": seq, "mutX": mutseq, "single_pred_psi": label, "mutY": mutlabel}


paths = ('spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(os.path.join('baselines/spliceai_models', x))
          for x in paths]
dataset = MutGenerator(context, sys.argv[1])
print("the size of dataset is {} model size is {}".format(
    len(dataset), len(models)))
with open(sys.argv[2], "w") as f:
    for idx, d in enumerate(dataset):
        if idx % 500 == 0:
            print(idx)
        x, mutx, gt, rpsi = d["X"], d["mutX"], d["mutY"], d["single_pred_psi"]
        inp = np.concatenate([x[None, :], mutx[None, :]], 0)
        assert len(inp.shape) == 3

        y = np.mean([models[m].predict(inp) for m in range(5)], axis=0)
        y, muty = y[0], y[1]
        assert len(gt.shape) == len(y.shape)
        refy = y[:, 1:].reshape(-1)
        muty = muty[:, 1:].reshape(-1)
        gt = gt[:, 1:].reshape(-1)
        rpsi = rpsi[:, 1:].reshape(-1)

        gt_idx = np.nonzero(gt)
        refy = refy[gt_idx].mean()
        muty = muty[gt_idx].mean()
        gt = gt[gt_idx].mean()
        rpsi = rpsi[gt_idx].mean()
        f.writelines("{}\t{}\t{}\t{}\n".format(rpsi, gt, refy, muty-refy))
        f.flush()
