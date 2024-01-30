from loguru import logger
import os
import sys
import json
from config import CL, AnnoPath,  EL, Fapath
from pyfasta import Fasta
SavePath = "data/clinvar"


def parse_line(line, i):
    line = line.replace("\n", "").split("\t")

    Name, _, Chr, Strand, Ts, Te, Js, Je, _ = line  # the data should be in 1-base
    deltaPsi = float("Nan")
    Psi = float("Nan")

    deltaPsi_vitro = None

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

    SaveData = []
    SaveData_vitro = []
    for idx, line in enumerate(mutinfo):
        NumTotal += 1
        Name, Chr, Strand, Ts, Te, Js, Je, Psi, deltaPsi, deltaPsi_vitro, Ori, Mut, Pos = line

        try:
            assert fafile[Chr][Pos].upper() == Ori.upper()
        except:
            print(1)
            continue
        center = Js
        start, end = max(center-padsize, 0), min(center +
                                                 padsize+1, len(fafile[Chr]))
        if end-start < padsize*2:
            continue
        assert end-start > 10000

        d = {"species": species, "chrom": Chr, "start": start, "end": end, "strand": Strand,
             "mutpos": Pos-start, "oriv": Ori, "mutv": Mut, "label": [], "mutlabel": [], "name": Name}
        d_vitro = {"species": species, "chrom": Chr, "start": start, "end": end, "strand": Strand,
                   "mutpos": Pos-start, "oriv": Ori, "mutv": Mut, "label": [], "mutlabel": [], "name": Name}

        if Strand == "+":
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
    File = "data/clinvar_GRCh37_20200316_splice.snps2inter_dataset.txt"
    if not os.path.exists(SavePath):
        os.mkdir(SavePath)

    main(File, "hg19", SavePath)
