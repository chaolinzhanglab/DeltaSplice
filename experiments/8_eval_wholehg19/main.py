from loguru import logger
import os
import torch
import os
from models.delta_pretrain import MainModel
from functools import partial
from pyfasta import Fasta
from config import Fapath, IN_MAP, repdict, SeqTable
import sys

SAVEPATH = "experiments/8_eval_wholehg19/test_results"
ModelPrefix = "RefSplice_models"
ModelPath = [os.path.join(ModelPrefix, x) for x in os.listdir(ModelPrefix)]


def main(species):
    # prepare models
    fasta = Fasta(Fapath+"/{}".format(species))
    SavePath=os.path.join(SAVEPATH, species)
    if not os.path.exists(SAVEPATH):
        os.mkdir(SAVEPATH)
    if not os.path.exists(SavePath):
        os.mkdir(SavePath)
    EL = 10000
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
    [_.eval() for _ in Models()]
    if not os.path.exists("PICS"):
        os.mkdir("PICS")

    Pred = []
    Mutv = []
    Gpred = []

    for chrom in fasta:
        if os.path.exists(SavePath+"/{}_res.csv".format(chrom)):
            continue
        save_file = open(SavePath+"/{}_res.csv".format(chrom), "w")
        save_file.writelines(
            "chrom, seq, strand,idx, pred_prob_class_0, pred_prob_class_1, pred_prob_class_2\n")

        print(chrom)
        chrom_length = len(fasta[chrom])
        starts = [_*5000 for _ in range(chrom_length//5000)]
        for start in starts:
            if start % 100000 == 0:
                print(start, chrom, chrom_length, flush=True, )
            ori_start = start
            seq_length = 35000//2
            start = start-15000
            end = start+35000+1
            seq = fasta[chrom][max(start, 0):min(
                end, chrom_length)].upper()
            if start < 0:
                seq = "".join(["N" for _ in range(-start)])+seq
            if end > chrom_length:
                seq = seq+"".join(["N" for _ in range(end-chrom_length)])
            # seq to tensor
            neg_seq = "".join([repdict[_] for _ in seq][::-1])
            str_seq = seq
            str_neg_seq = neg_seq

            seq = IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
            neg_seq = IN_MAP[[SeqTable[_] for _ in neg_seq]][:, :4]
            seq = torch.tensor(seq).cuda()[None]
            neg_seq = torch.tensor(neg_seq).cuda()[None]
            d = {"X": torch.cat([seq, neg_seq], 0)}
            TPRED = [m.predict(d) for m in Models]
            gpred = sum([v["single_pred_psi"] for v in TPRED])/len(TPRED)
            neg_gpred = gpred[1][::-1]
            gpred = gpred[0]
            for index in range(ori_start, ori_start+5000):
                if fasta[chrom][index:index+2].upper() in ["AG", "GT", "CT", "AC"]:
                    if fasta[chrom][index:index+2].upper() in ["AG", "GT"]:
                        idx = index-ori_start
                        assert str_seq[idx+15000:idx+2+15000] in ["AG",
                                                                  "GT"], str_seq[idx+15000:idx+2+15000]
                        if sum(gpred[idx][1:].tolist()) > 1e-3:
                            save_file.writelines("{},{},{},{},{},{},{}\n".format(
                                chrom, str_seq[idx+15000],"+", index, *gpred[idx].tolist()))
                        if sum(gpred[idx+1][1:].tolist()) > 1e-3:
                            save_file.writelines("{},{},{},{},{},{},{}\n".format(
                                chrom, str_seq[idx+15000+1],"+", index+1, *gpred[idx+1].tolist()))
                    else:
                        idx = index-ori_start
                        assert str_seq[idx+15000:idx+2+15000] in ["CT",
                                                                  "AC"], str_seq[idx+15000:idx+2+15000]
                        if sum(neg_gpred[idx][1:].tolist()) > 1e-3:
                            save_file.writelines("{},{},{},{},{},{},{}\n".format(
                                chrom, str_seq[idx+15000],"-", index, *neg_gpred[idx].tolist()))
                        if sum(neg_gpred[idx+1][1:].tolist()) > 1e-3:
                            save_file.writelines("{},{},{},{},{},{},{}\n".format(
                                chrom, str_seq[idx+15000+1],"-", index+1, *neg_gpred[idx+1].tolist()))
        save_file.close()


if __name__ == "__main__":
    main(sys.argv[1])
