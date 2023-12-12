import sys
from scipy.stats import spearmanr, pearsonr
from loguru import logger
import numpy as np
from sklearn import metrics


def Cors(x, y):
    return spearmanr(x, y), pearsonr(x, y)


def loadpred(path):
    with open(path, "r") as f:
        content = f.readlines()
    content = [x.replace("\n", "").split("\t") for x in content if x[0] != "#"]
    Mgt = []
    Mpred = []
    for line in content:
        mpred, mgt = line[-2:]

        Mgt.append(np.float(mgt))
        Mpred.append(np.float(mpred))
    return Mpred, Mgt


def parse_line(line, i):
    line = line.replace("\n", "").split("\t")
    try:
        Name, _, Chr, Strand, Ts, Te, Js, Je, Psi, deltaPsi = line
        deltaPsi = float(deltaPsi)
        if deltaPsi == 0:
            deltaPsi = 1e-5
    except:
        Name, _, Chr, Strand, Ts, Te, Js, Je, deltaPsi = line
        deltaPsi = 1.0
        Psi = 1.0
    Pos, Ori, Mut = Name.split("_")[1:]
    Name = "{}_{}".format(i, Name)

    return Name, Chr, Strand, int(Ts)-1, int(Te)-1, int(Js)-1, int(Je)-1, float(Psi), float(deltaPsi), Ori, Mut, int(Pos)-1


def loadgt(path):
    with open(path, "r") as f:
        content = f.readlines()
    content = [x for x in content if x[0] != "#"]
    content = [parse_line(_, i) for i, _ in enumerate(content)]
    name = [x[0].split("_")[1] for x in content]
    content = [x[7:9] for x in content]
    Gt = []
    MGt = []
    for line in content:
        Gt.append(line[0])
        MGt.append(line[1])
    # print(name)
    return Gt, MGt, name


def topk(GT, Pred, k=50):
    thresholds = 0.5
    '''tgt, tpred = GT, Pred
    GT, Pred = [], []
    for g, p in zip(tgt, tpred):
        if g > thresholds:
            GT.append(g)
            Pred.append(p)
        elif g < -thresholds:
            GT.append(g)
            Pred.append(p)'''
    GT = np.array(GT)
    Pred = np.array(Pred)
    Pred, GT = np.abs(Pred), np.abs(GT)

    fpr, tpr, _ = metrics.roc_curve(GT >= thresholds, Pred)

    logger.info("ACC & AUC & AUPRC  {} & {} & {} with a threshold of {}".format(((GT >= thresholds) == (Pred >= thresholds)).astype(
        float).sum()/len(GT), metrics.auc(fpr, tpr), metrics.average_precision_score(GT >= thresholds, Pred), thresholds))
    Pred, GT = np.abs(Pred), np.abs(GT)
    idxs = np.argsort(Pred)
    GT = GT[idxs]
    logger.info("Top {} Average Value {}".format(k, GT[-k:].mean()))


def main():
    Mpred, PMgt = loadpred(sys.argv[1])
    if len(sys.argv) > 2:
        threshold = float(sys.argv[2])
    else:
        threshold = 0.05

    logger.info("The cor for mut is {}".format(Cors(PMgt, Mpred)))
    logger.info("The cor for mut with gt larger than {} is {} with number of {}".format(threshold, Cors([a for a in PMgt if abs(a) > threshold], [
                a for a, b in zip(Mpred, PMgt) if abs(b) > threshold]), len([a for a in PMgt if abs(a) > threshold])))
    topk(PMgt, Mpred)


if __name__ == "__main__":
    main()
