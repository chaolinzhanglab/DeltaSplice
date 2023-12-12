from loguru import logger
import os
import sys
import json
from config import CL, AnnoPath,  EL, Fapath
from bisect import bisect_left
import numpy as np
from pyfasta import Fasta


def parse_line(line, i):
    line = line.replace("\n", "").split("\t")
    try:
        Name, _, Chr, Strand, Ts, Te, Js, Je, Psi, deltaPsi = line  # the data should be in 1-base
        deltaPsi = float(deltaPsi)
        if deltaPsi == 0:
            deltaPsi = 1e-5
        deltaPsi_vitro = None
    except:
        Name, _, Chr, Strand, Ts, Te, Js, Je, Psi, deltaPsi, deltaPsi_vitro = line[:11]  # the data should be in 1-base
        deltaPsi = float(deltaPsi)
        Psi = float("nan")
        if deltaPsi == 0:
            deltaPsi = 1e-5
        deltaPsi_vitro = float(deltaPsi_vitro)
        if deltaPsi_vitro == 0:
            deltaPsi_vitro = 1e-5

    Pos, Ori, Mut = Name.split(";")[0].split("_")[1:]
    Name = "{}_{}".format(i, Name)

    return Name, Chr, Strand, int(Ts)-1, int(Te)-1, int(Js)-1, int(Je)-1, float(Psi), float(deltaPsi), deltaPsi_vitro, Ori, Mut, int(Pos)-1


def get_gene_table(File):
    with open(File, "r") as f:
        content = f.readlines()
    content = [x for x in content if x[0] != "#" and x[0] != "n"]
    content = [parse_line(_, i) for i, _ in enumerate(content)]
    return content


def main(File, species, SavePath):
    mutinfo = get_gene_table(File)
    fafile = Fasta(os.path.join(Fapath, species+".fa"))

    padsize = (EL+CL)//2
    numNonexp, NumTotal = 0, 0
    with open(os.path.join(AnnoPath, "data.json"), "r") as f:
        Annotation = json.load(f)
    SortedKeys = {}
    for sp in Annotation:
        SortedKeys[sp] = {}
        for chr in Annotation[sp]:
            SortedKeys[sp][chr] = {}
            SortedKeys[sp][chr]["+"] = sorted([int(_) for _ in Annotation[sp][chr]["+"].keys()])
            SortedKeys[sp][chr]["-"] = sorted([int(_) for _ in Annotation[sp][chr]["-"].keys()])
    SaveData = []
    SaveData_vitro = []
    for idx, line in enumerate(mutinfo):
        NumTotal += 1
        Name, Chr, Strand, Ts, Te, Js, Je, Psi, deltaPsi, deltaPsi_vitro, Ori, Mut, Pos = line

        center = Js
        start, end = max(center-padsize, 0), min(center+padsize+1, len(fafile[Chr]))
        assert end-start>10000

        d = {"species": species, "chrom": Chr, "start": start, "end": end, "strand": Strand, "mutpos": Pos-start, "oriv": Ori, "mutv": Mut, "label": [], "mutlabel": [], "name": Name}
        d_vitro = {"species": species, "chrom": Chr, "start": start, "end": end, "strand": Strand, "mutpos": Pos-start, "oriv": Ori, "mutv": Mut, "label": [], "mutlabel": [], "name": Name}

        startidx, endidx = bisect_left(SortedKeys[species][Chr][Strand], Js), bisect_left(SortedKeys[species][Chr][Strand], Je)
        try:
            assert startidx < len(SortedKeys[species][Chr][Strand]) and Js == SortedKeys[species][Chr][Strand][startidx]
            assert endidx < len(SortedKeys[species][Chr][Strand]) and Je == SortedKeys[species][Chr][Strand][endidx]
        except:
            logger.warning("The junction site not exists in annotation file")

        if Strand == "+":
            try:
                assert float(Annotation[species][Chr]["+"][str(Js)][2]) > 0 or np.isnan(float(Annotation[species][Chr]["+"][str(Js)][2]))
                assert float(Annotation[species][Chr]["+"][str(Je)][1]) > 0 or np.isnan(float(Annotation[species][Chr]["+"][str(Je)][1]))
            except:
                numNonexp += 1
                logger.info("Encounter nonexpression sites of {}/{}, will not pass ...".format(numNonexp, NumTotal))
            d["label"].append([Js-start, [0, 0, Psi]])
            d["label"].append([Je-start, [0, Psi, 0]])
            d["mutlabel"].append([Js-start, [0, 0, deltaPsi]])
            d["mutlabel"].append([Je-start, [0, deltaPsi, 0]])

            d_vitro["label"].append([Js-start, [0, 0, Psi]])
            d_vitro["label"].append([Je-start, [0, Psi, 0]])
            d_vitro["mutlabel"].append([Js-start, [0, 0, deltaPsi_vitro]])
            d_vitro["mutlabel"].append([Je-start, [0, deltaPsi_vitro, 0]])
        else:
            assert Strand == "-"
            # print(Annotation[species][Chr]["-"][str(Je)])
            try:
                assert float(Annotation[species][Chr]["-"][str(Js)][1]) > 0 or np.isnan(float(Annotation[species][Chr]["-"][str(Js)][1]))
                assert float(Annotation[species][Chr]["-"][str(Je)][2]) > 0 or np.isnan(float(Annotation[species][Chr]["-"][str(Je)][2]))
            except:
                numNonexp += 1
                logger.info("Encounter nonexpression sites of {}/{}, will not pass ...".format(numNonexp, NumTotal))
            d["label"].append([Js-start, [0, Psi, 0]])
            d["label"].append([Je-start, [0, 0, Psi]])
            d["mutlabel"].append([Js-start, [0, deltaPsi, 0]])
            d["mutlabel"].append([Je-start, [0, 0, deltaPsi]])

            d_vitro["label"].append([Js-start, [0, Psi, 0]])
            d_vitro["label"].append([Je-start, [0, 0, Psi]])
            d_vitro["mutlabel"].append([Js-start, [0, deltaPsi_vitro, 0]])
            d_vitro["mutlabel"].append([Je-start, [0, 0, deltaPsi_vitro]])
        SaveData.append(d)
        SaveData_vitro.append(d_vitro)
    logger.info("finish handling {}".format(len(SaveData)))
    with open(os.path.join(SavePath, "data.json"), "w") as f:
        f.write(json.dumps(SaveData))
    if deltaPsi_vitro is not None:
        with open(os.path.join(SavePath, "data_vitro.json"), "w") as f:
            f.write(json.dumps(SaveData_vitro))


if __name__ == "__main__":
    import sys
    File, SavePath, species = sys.argv[1:4]
    if not os.path.exists(SavePath):
        os.mkdir(SavePath)
    main(File, species, SavePath)
