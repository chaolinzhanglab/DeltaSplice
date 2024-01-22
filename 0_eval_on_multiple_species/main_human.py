from loguru import logger
import os
import numpy as np
import torch
from torch.utils.data import Dataset
import json
import os
from config import repdict, SeqTable, IN_MAP
from models.delta_pretrain import MainModel
from functools import partial
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["axes.labelsize"] = 18
plt.rc("legend", fontsize=18)
plt.rc("figure", titlesize=20)
plt.rc("font", size=18)


DataPath = "data/Hg19VsOthers/"
SavePath = "experiments/0_eval_on_multiple_species/test_results"
ModelPrefix = "RefSplice_models_human"
ModelPath = [os.path.join(ModelPrefix, x) for x in os.listdir(ModelPrefix)]
Species = ["susScr11", "mm10", "rheMac10", "rn6", "panTro5", "bosTau9"]
revseqtable = ["A", "C", "G", "T"]


class DataGenerator(Dataset):
    def __init__(self, s):
        with open(os.path.join(DataPath, "{}.json".format(s)), "r") as f:
            self.data = json.load(f)
        self.s = s

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d = self.data[index]
        hg19_seq = d["hg19_seq"]
        transition_seq = d["transition"]
        hg19_strand = d["hg19_strand"]
        end_seq = d["{}_seq".format(self.s)]
        end_strand = d["{}_strand".format(self.s)]
        hg19_label = d["hg19_label"]
        end_label = d["{}_label".format(self.s)]
        cidx = len(end_seq)//2
        if hg19_strand == "-":
            hg19_seq = "".join([repdict[_] for _ in hg19_seq.upper()][::-1])
            transition_seq = [
                "".join([repdict[_] for _ in x.upper()][::-1]) for x in transition_seq]
        if end_strand == "-":
            end_seq = "".join([repdict[_] for _ in end_seq.upper()][::-1])
        distance = []
        '''for s in transition_seq:
            v = 0
            for i in range(cidx-200, cidx+200+1):
                if s[i].upper() != end_seq[i].upper():
                    v += 1
            distance.append(v)'''
        hg19_seq = IN_MAP[np.array([SeqTable[_] for _ in hg19_seq.upper()])]
        end_seq = IN_MAP[np.array([SeqTable[_] for _ in end_seq.upper()])]
        transition_seq = [
            IN_MAP[np.array([SeqTable[_] for _ in x.upper()])] for x in transition_seq]

        return {"hg19_seq": hg19_seq, "transition_seq": [hg19_seq, end_seq], "end_seq": end_seq, "hg19_label": [float(_) for _ in hg19_label], "end_label": [float(_) for _ in end_label], "distance": distance, "gene": d["hg19_name"], "chr": d["hg19_chrom"], "center": d["hg19_idx"]}


def main():
    # prepare models
    EL = 30000
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

    Models = []
    for m in ModelPath:
        model = MainModel(64, W, AR, 0.3, EL=EL,
                          optimizer=partial(torch.optim.Adam, lr=1e-3))
        model.load_state_dict(torch.load(m))
        Models.append(model)
    if not os.path.exists("PICS"):
        os.mkdir("PICS")

    for species in Species:
        if os.path.exists(os.path.join(SavePath, "human_{}_pred".format(species))):
            continue
        with open(os.path.join(SavePath, "human_{}_pred".format(species)), "w") as f:
            f.write("start")
        data = DataGenerator(species)
        logger.info("The number of data is {}".format(len(data)))

        hg19_pred = []
        hg19_gt = []
        end_pred = []
        end_gt = []
        BS = 16

        def eval(model, target, ref, reflabel):
            n = (len(target)+BS-1)//BS
            ret = []
            for i in range(n):
                ret.append(model.encode(
                    target[BS*i:BS*i+BS], ref[BS*i:BS*i+BS], reflabel[BS*i:BS*i+BS]))
            return torch.cat(ret, 0)

        with torch.no_grad():
            print(len(data))
            for i in range(len(data)):
                d = data[i]
                if i % 100 == 0:
                    print(i, flush=True)
                inp = torch.tensor(d["transition_seq"]).cuda().float()
                reflabel = torch.tensor(np.array(d["hg19_label"]).reshape(
                    1, 1, 3)).cuda().float()+torch.zeros_like(inp[:, :inp.shape[1]-EL, :3])
                reflabel[:, :, 0] = 1-reflabel[:, :, 1:].sum(-1)
                pred = sum([eval(m, inp[:, :, :4], inp[:1, :, :4]+torch.zeros_like(
                    inp[:, :, :4]), reflabel) for m in Models])/len(Models)
                assert inp.shape[1] % 2 == 1
                idx = pred.shape[1]//2
                pred = pred[:, idx].detach().cpu().numpy()
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
                N = 50
                assert len(hg19_gt) == i+1

                '''if len(pred) > 2:
                    fig=plt.figure(figsize=[16,6])
                    ax1=fig.add_subplot("111")
                    ax2=ax1.twiny()
                    hgseq=d["transition_seq"][0].argmax(-1)
                    targetseq=d["transition_seq"][-1].argmax(-1)
                    idx=len(hgseq)//2
                    hgseq=[revseqtable[_] for _ in hgseq[idx-N//2:idx+N//2+1]]
                    targetseq=[revseqtable[_] for _ in targetseq[idx-N//2:idx+N//2+1]]
                    for j in range(len(d["transition_seq"])-2):
                        seq=d["transition_seq"][j+1].argmax(-1)[idx-N//2:idx+N//2+1]
                        seq=[revseqtable[_] for _ in seq]
                        for k in range(len(seq)):
                            if seq[k]!=hgseq[k]:
                                assert seq[k]==targetseq[k]
                                ax1.scatter([k], [pred[j+1]], c="blue", s=100)
                    ax1.plot(range(N), [pred[0] for _ in range(N)], label="Predicted usage of the hg19 sequence", c="r", linestyle="dotted")
                    ax1.plot(range(N), [pred[-1] for _ in range(N)], label="Predicted usage of the {} sequence".format(species), c="black", linestyle="dotted")
                    ax1.plot(range(N), [hg19_gt[-1] for _ in range(N)], label="Real usage of the hg19 sequence", c="r")
                    ax1.plot(range(N), [end_gt[-1] for _ in range(N)], label="Real usage of the {} sequence".format(species), c="black")
                    ax1.set_xlabel("Sequence of hg19 {}:{}-{}".format(d["chr"], d["center"]-N//2, d["center"]+N//2))
                    ax2.set_xlabel("Homologous sequence on {}".format(species))
                    ax2.set_xlim(ax1.get_xlim())
                    ax2.set_xticks(range(N+1))
                    ax2.set_xticklabels(targetseq)
                    ax1.set_xticks(range(N+1))
                    ax1.set_xticklabels(hgseq)
                    ax1.legend()
                    plt.savefig("PICS/{}_{}_{}.png".format(species, i, d["gene"]))
                    plt.close()'''
        assert len(hg19_gt) == len(data)

        with open(os.path.join(SavePath, "human_{}_pred".format(species)), "w") as f:
            for a, b, c, d in zip(hg19_pred, hg19_gt, end_pred, end_gt):
                f.writelines("{}\t{}\t{}\t{}\n".format(a, b, c, d))


if __name__ == "__main__":
    main()
