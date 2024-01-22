from loguru import logger
import os
import numpy as np
import torch
from torch.utils.data import Dataset
import json
import os
from config import repdict, SeqTable, IN_MAP
from models.pangolin import MainModel
from functools import partial
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
DataPath = "data/Hg19VsOthers/"
SavePath = "experiments/0_eval_on_multiple_species/test_results"
ModelPath = [
    "pangolin_human_models/model.ckpt-rep{}".format(i) for i in range(5)]
Species = ["susScr11", "mm10", "rheMac10", "rn6", "panTro5", "bosTau9"]


class DataGenerator(Dataset):
    def __init__(self, s):
        with open(os.path.join(DataPath, "{}.json".format(s)), "r") as f:
            self.data = json.load(f)
        self.s = s
        print(len(self.data), s)

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
        for s in transition_seq:
            v = 0
            for i in range(cidx-200, cidx+200+1):
                if s[i].upper() != end_seq[i].upper():
                    v += 1
            distance.append(v)
        hg19_seq = IN_MAP[np.array([SeqTable[_] for _ in hg19_seq.upper()])]
        end_seq = IN_MAP[np.array([SeqTable[_] for _ in end_seq.upper()])]
        transition_seq = [
            IN_MAP[np.array([SeqTable[_] for _ in x.upper()])] for x in transition_seq]

        return {"hg19_seq": hg19_seq, "transition_seq": [hg19_seq, end_seq], "end_seq": end_seq, "hg19_label": [float(_) for _ in hg19_label], "end_label": [float(_) for _ in end_label], "distance": distance}


def main():
    # prepare models
    L = 32
    # convolution window size in residual units
    W = [11, 11, 11, 11, 11, 11, 11, 11,
         21, 21, 21, 21, 41, 41, 41, 41]
    # atrous rate in residual units
    AR = [1, 1, 1, 1, 4, 4, 4, 4,
          10, 10, 10, 10, 25, 25, 25, 25]
    EL = 30000

    Models = []
    for m in ModelPath:
        model = MainModel(32, W, AR, 0, EL=EL,
                          optimizer=partial(torch.optim.Adam, lr=1e-3))
        model.load_state_dict(torch.load(m))
        Models.append(model)
    if not os.path.exists("PICS"):
        os.mkdir("PICS")
    for species in Species:
        data = DataGenerator(species)
        logger.info("The number of data is {}".format(len(data)))

        hg19_pred = []
        hg19_gt = []
        end_pred = []
        end_gt = []

        with torch.no_grad():
            for i, d in enumerate(data):
                if i % 100 == 0:
                    print(i)
                inp = torch.tensor(d["transition_seq"]).cuda().float()
                pred = sum([m.encode(inp[:, :, :4])[1]
                           for m in Models])/len(Models)
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
                '''if len(pred)>2:
                    plt.figure()
                    plt.scatter(d["distance"], pred, label="{}/{}".format(round(sum(hg19_label), 4), round(sum(end_label), 4)))
                    plt.scatter(d["distance"][:1], [sum(hg19_label)], c="r")
                    plt.scatter([0], [sum(end_label)], c="r")
                    plt.legend()
                    plt.savefig("PICS/{}_{}.png".format(species,i))
                    plt.close()'''

        """plt.figure(figsize=[8*3, 6*1])
        f = plt.subplot(1, 3, 1)
        f.scatter(hg19_gt, hg19_pred, label="$R^2$={}".format(pearsonr(hg19_pred, hg19_gt)[0]))
        f.legend()
        f.set_xlabel("Gt of hg19")
        f.set_ylabel("Pred of hg19")
        f = plt.subplot(1, 3, 2)
        f.scatter(end_gt, end_pred, label="$R^2$={}".format(pearsonr(end_pred, end_gt)[0]))
        f.legend()
        f.set_xlabel("Gt of {}".format(species))
        f.set_ylabel("Pred of {}".format(species))
        f = plt.subplot(1, 3, 3)
        f.scatter([a-b for a, b in zip(hg19_gt, end_gt)], [a-b for a, b in zip(hg19_pred, end_pred)],
                  label="$R^2$={}".format(pearsonr([a-b for a, b in zip(hg19_gt, end_gt)], [a-b for a, b in zip(hg19_pred, end_pred)])[0]))
        f.legend()
        f.set_xlabel("Real $\Delta$usage")
        f.set_ylabel("Predicted $\Delta$usage")

        plt.savefig("PICS/pangolin_human_merge_{}.png".format(species), bbox_inches="tight")
        plt.close()"""

        with open(os.path.join(SavePath, "pangolin_human_{}_pred".format(species)), "w") as f:
            for a, b, c, d in zip(hg19_pred, hg19_gt, end_pred, end_gt):
                f.writelines("{}\t{}\t{}\t{}\n".format(a, b, c, d))


if __name__ == "__main__":
    main()
