from pangolin_model import *
import numpy as np
import json
from pyfasta import Fasta
import os
import sys
IN_MAP = np.asarray([[0, 0, 0, 0],
                     [1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])

Fapath = "hgfile"
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
context = 10000
model_nums = [0, 2, 4, 6]
modelpath = "baselines/pangolin_models"
INDEX_MAP = {0: 1, 1: 2, 2: 4, 3: 5, 4: 7, 5: 8, 6: 10, 7: 11}
Table = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 0}


def one_hot_encode(seq):
    seq = seq.upper()
    seq = [Table[_] for _ in seq]
    return IN_MAP[seq]


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
        if "mutpos" in d and "species" in d:
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


models = []
for i in model_nums:
    for j in range(1, 6):
        model = Pangolin(L, W, AR)
        if torch.cuda.is_available():
            model.cuda()
        model.load_state_dict(torch.load(
            os.path.join(modelpath, "final.%s.%s.3" % (j, i))))
        model.eval()
        models.append(model)
dataset = MutGenerator(context, sys.argv[1])
print("the size of dataset is {} model size is {}".format(
    len(dataset), len(models)))
with open(sys.argv[2], "w") as f:
    for idx, d in enumerate(dataset):
        if idx % 500 == 0:
            print(idx)
        x, mutx, gt, rpsi = d["X"], d["mutX"], d["mutY"], d["single_pred_psi"]
        inp = np.concatenate([x[None, :], mutx[None, :]], 0)
        seq = torch.from_numpy(inp).permute(0, 2, 1).float()
        if torch.cuda.is_available():
            seq = seq.to(torch.device("cuda"))

        ref, mut = [], []
        gt = gt[:, 1:].sum(-1).reshape(-1)
        gt_idx = np.nonzero(gt)
        for j, model_num in enumerate(model_nums):
            with torch.no_grad():
                y = [m(seq)[:, INDEX_MAP[model_num], :].cpu().numpy()
                     for m in models[5*j:5*j+5]]
                y = np.mean(y, axis=0)
                ref_score = y[0][gt_idx].mean()
                mut_score = y[1][gt_idx].mean()
                for a, b in zip(y[0].shape, gt.shape):
                    assert a == b
                ref.append(ref_score)
                mut.append(mut_score)

        ref = np.array(ref).reshape(-1)
        mut = np.array(mut).reshape(-1)
        assert len(ref) == len(mut)
        assert len(ref) == 4
        ref_score = ref[np.argmax(abs(mut-ref), axis=0)]
        mut_score = mut[np.argmax(abs(mut-ref), axis=0)]

        gt = gt[gt_idx].mean()
        refgt = rpsi[:, 1:].sum(-1).reshape(-1)[gt_idx].mean()
        f.writelines("{}\t{}\t{}\t{}\n".format(
            refgt, gt, ref_score, mut_score-ref_score))
        f.flush()
