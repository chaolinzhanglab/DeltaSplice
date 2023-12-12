import json
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import pandas
GlobalDict = {}


def load_json(Path):
    with open(Path, "r") as f:
        return json.load(f)


def load_txt(Path):
    with open(Path, "r") as f:
        content = f.readlines()
    content = [list(map(float, x.strip().split("\t"))) for x in content]
    return content


def plot(jsonFile, txtFile, name):
    Json = load_json(jsonFile)
    Pred = load_txt(txtFile)
    assert len(Json) == len(Pred)
    GTL = []
    GTD = []
    PredD = []
    for G, P in zip(Json, Pred):
        predl, predd, gtd = P
        gtlabel = G["label"]
        mutseq = G["mutseq"]

        gtv = 0
        gtd = 0
        for v in gtlabel:
            gtv += sum(v[1])
        for v in G["mutlabel"]:
            gtd += sum(v[1])
        # if abs(gtd+gtv - len(gtlabel)) < 1e-3 or abs(gtv - len(gtlabel)) < 1e-3:
        #    continue

        gtv = gtv/len(gtlabel)
        GTL.append(gtv)
        PredD.append(predd)
        GTD.append(gtd/len(gtlabel))
        if mutseq not in GlobalDict:
            GlobalDict[mutseq] = [gtd/len(gtlabel)+gtv, []]
        GlobalDict[mutseq][1].append([predd+predl, gtv])

    plt.figure(figsize=[8*2, 6])
    f1, f2 = plt.subplot(1, 2, 1), plt.subplot(1, 2, 2)
    f1.scatter(GTL, GTD, s=2., label="GT", c="r")
    f1.scatter(GTL, PredD, s=2., label="Pred", c="black")
    f1.legend(fontsize=15)
    f1.set_xlabel("Start PSI", fontsize=12)
    f1.set_ylabel("$\Delta$PSI", fontsize=12)

    #f2.set_xlim([-1, 1])
    #f2.set_ylim([-1, 1])
    f2.scatter(GTD, PredD, s=2., c="black")
    sp = round(spearmanr(GTD, PredD)[0], 3)
    pear = round(pearsonr(GTD, PredD)[0], 3)
    #f2.text(0, 0, "SR={},PR={}".format(sp, pear), fontsize=20)
    f2.set_xlabel("GT $\Delta$PSI", fontsize=12)
    f2.set_ylabel("Pred $\Delta$PSI", fontsize=12)

    print(spearmanr(GTD, PredD), pearsonr(GTD, PredD))
    plt.savefig("Pic/{}.png".format(name))


def plot_global():
    acceptor = pandas.read_csv("example_merged.json+evaltasks_Pretrain_withoutcons_rep0_models_model.ckpt-2_acceptor", sep="\t")
    donor = pandas.read_csv("example_merged.json+evaltasks_Pretrain_withoutcons_rep0_models_model.ckpt-2_donor", sep="\t")
    X1 = acceptor["Yt"]
    Y1 = (acceptor["Yp"]+donor["Yp"])*0.5
    l = len(GlobalDict)
    m, n = (l+3)//4, 4
    F = plt.figure(figsize=(3*8, 6))
    idx = 1
    X, Y = [], []
    tpx, tpy = [], []
    for key in GlobalDict:
        gv, pv = GlobalDict[key]
        idx += 1
        pvx, pvy = [x[0] for x in pv], [x[1] for x in pv]
        if len(pvx) > len(tpx):
            tpx, tpy = pvx, pvy
        meanv = sum(pvy)/len(pvy)
        X.append(gv)
        Y.append(meanv)

    sp = round(spearmanr(X, Y)[0], 3)
    pear = round(pearsonr(X, Y)[0], 3)
    sp1 = round(spearmanr(X1, Y1)[0], 3)
    pear1 = round(pearsonr(X1, Y1)[0], 3)

    f1, f2, f3 = plt.subplot(1, 3, 1), plt.subplot(1, 3, 2), plt.subplot(1, 3, 3)
    f3.scatter(tpx, tpy, s=2., c="black")
    f1.scatter(X, Y, s=2., c="black")
    f2.scatter(X1, Y1, s=2., c="black")
    print(sp, pear)
    f1.text(0.2, 0.05, "SR={}, PR={}".format(sp, pear), fontsize=20)
    f2.text(0.2, 0.9, "SR={}, PR={}".format(sp1, pear1), fontsize=20)
    f1.set_xlim([0, 1])
    f1.set_ylim([0, 1])
    f2.set_xlim([0, 1])
    f2.set_ylim([0, 1])
    f1.set_xlabel("Real PSI")
    f1.set_ylabel("Pred PSI")
    f2.set_xlabel("Real PSI")
    f2.set_ylabel("Pred PSI")

    plt.savefig("Pic/Total.png")


if __name__ == "__main__":
    import os
    for v in os.listdir("FAS_exon6"):
        if not v.endswith("json"):
            continue
        jsonfile = "FAS_exon6/{}".format(v)
        predfile = "Mutation__temp_xuchencheng_eSplicedata_Species_37_ProcessedData_FAS_exon6_{}+tasks_Pretrain_withoutcons_rep0_models_model.ckpt-2.txt".format(v)
        plot(jsonfile, predfile, v.replace(".json", ""))
    plot_global()
