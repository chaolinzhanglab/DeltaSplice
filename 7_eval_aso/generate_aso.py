from loguru import logger
import os
import sys
import json
from config import CL, AnnoPath,  EL, repdict, Fapath
from bisect import bisect_left
import numpy as np
from pyfasta import Fasta
import pandas
import copy


def refine_with_info(sequence, infofile):
    data = pandas.read_csv(infofile, sep="\t")
    ret = []
    # preprocess: mean all the ctrl.psi and dPSI of the same sequence

    for ctrl_psi, dpsi, aso1, aso2, aso3, aso4 in zip(data["ctrl.psi"], data["dPSI"], data["ext0_mask.seq.name"], data["ext1_mask.seq.name"], data["ext2_mask.seq.name"], data["ext3_mask.seq.name"]):
        for aso in [aso1, aso2, aso3, aso4]:
            if not isinstance(aso, str):
                continue
            d = copy.deepcopy(sequence[aso])
            if dpsi == 0:
                dpsi = dpsi+1e-3
            if ctrl_psi == 0:
                ctrl_psi = ctrl_psi+1e-3

            for x in d["label"]:
                for i in range(len(x[1])):
                    if np.isnan(x[1][i]):
                        x[1][i] = ctrl_psi/100.
            for x in d["mutlabel"]:
                for i in range(len(x[1])):
                    if np.isnan(x[1][i]):
                        x[1][i] = dpsi/100.
            ret.append(d)
    return ret


def parse_line(line, i):
    line = line.replace("\n", "").split("\t")
    # the data should be in 1-base
    Name, _, Chr, Strand, Ts, Te, Js, Je = line[:8]
    return Name, Chr, Strand, int(Ts)-1, int(Te)-1, int(Js)-1, int(Je)-1


def get_seq(File):
    with open(File, "r") as f:
        content = f.readlines()
    content = [x.replace("\n", "").split("\t")[1].upper() for x in content]
    return content


def get_gene_table(File):
    with open(File, "r") as f:
        content = f.readlines()
    content = [x for x in content if x[0] != "#" and x[0] != "n"]
    content = [parse_line(_, i) for i, _ in enumerate(content)]
    return content


def main(File, seqFile, infoFile, species, SavePath, hgpath, ref_values):
    mutinfo = get_gene_table(File)
    hgFile = Fasta(hgpath)

    if os.path.exists(seqFile):
        preseqs = get_seq(seqFile)
    else:
        logger.info("Not using seqfile")
        preseqs = None

    padsize = (EL+CL)//2
    numNonexp, NumTotal = 0, 0
    with open(os.path.join(AnnoPath, "data.json"), "r") as f:
        Annotation = json.load(f)
    SortedKeys = {}
    for sp in Annotation:
        SortedKeys[sp] = {}
        for chr in Annotation[sp]:
            SortedKeys[sp][chr] = {}
            SortedKeys[sp][chr]["+"] = sorted([int(_)
                                              for _ in Annotation[sp][chr]["+"].keys()])
            SortedKeys[sp][chr]["-"] = sorted([int(_)
                                              for _ in Annotation[sp][chr]["-"].keys()])
    SaveData = []
    SaveData_vitro = []
    NumNoInfo = 0
    for idx, line in enumerate(mutinfo):
        NumTotal += 1
        Name, Chr, Strand, Ts, Te, Js, Je = line
        Psi, deltaPsi = float("nan"), float("nan")

        center = (Js+Je)//2+1
        start, end = center-padsize, center+padsize+1
        seq = hgFile[Chr][start:end].upper()

        if preseqs is not None:
            assert len(preseqs[idx]) == 10500
            bpadsize = (len(seq)-10500)//2
            if idx == 0:
                baseseq = seq[bpadsize:bpadsize+10500]
                assert preseqs[idx][:100] == baseseq[:100]
                baseseq = preseqs[idx]
                seq = seq[:bpadsize]+baseseq+seq[bpadsize+10500:]
                assert seq[bpadsize:bpadsize +
                           10500].upper() == preseqs[idx].upper()
                assert len(seq) == bpadsize*2+10500+1
            else:
                bs, be = map(int, Name.split("||p")[-1].split("-"))
                seq = seq[:bpadsize]+baseseq[:bs] + \
                    "".join(["N" for _ in range(be-bs+1)]) + \
                    baseseq[be+1:]+seq[bpadsize+10500:]
                assert len(seq) == bpadsize*2+10500+1
                assert seq[bpadsize:bpadsize +
                           10500].upper() == preseqs[idx].upper()
        elif idx > 0:
            bpadsize = (35000-10500)//2
            bs, be = map(int, Name.split("||p")[-1].split("-"))
            seq = seq[:bpadsize+bs] + \
                "".join(["N" for _ in range(be-bs+1)])+seq[be+1+bpadsize:]
            assert len(seq) == bpadsize*2+10500+1

        if idx == 0:
            nonmutseq = seq
        d = {"species": species, "chrom": Chr, "start": start, "end": end,
             "strand": Strand,  "label": [], "mutlabel": [], "name": Name}
        d_vitro = {"species": species, "chrom": Chr, "start": start, "end": end,
                   "strand": Strand, "label": [], "mutlabel": [], "name": Name}
        d["mutseq"] = seq
        d["seq"] = nonmutseq
        startidx, endidx = bisect_left(SortedKeys[species][Chr][Strand], Js), bisect_left(
            SortedKeys[species][Chr][Strand], Je)
        try:
            assert startidx < len(
                SortedKeys[species][Chr][Strand]) and Js == SortedKeys[species][Chr][Strand][startidx]
            assert endidx < len(
                SortedKeys[species][Chr][Strand]) and Je == SortedKeys[species][Chr][Strand][endidx]
        except:
            logger.warning("The junction site not exists in annotation file")

        if Strand == "+":
            try:
                assert float(Annotation[species][Chr]["+"][str(Js)][2]) > 0 or np.isnan(
                    float(Annotation[species][Chr]["+"][str(Js)][2]))
                assert float(Annotation[species][Chr]["+"][str(Je)][1]) > 0 or np.isnan(
                    float(Annotation[species][Chr]["+"][str(Je)][1]))
            except:
                numNonexp += 1
                logger.info(
                    "Encounter nonexpression sites of {}/{}, will not pass ...".format(numNonexp, NumTotal))
            d["label"].append([Js-start, [0, 0, ref_values]])
            d["label"].append([Je-start, [0, ref_values, 0]])
            d["mutlabel"].append([Js-start, [0, 0, deltaPsi]])
            d["mutlabel"].append([Je-start, [0, deltaPsi, 0]])

            d_vitro["label"].append([Js-start, [0, 0, Psi]])
            d_vitro["label"].append([Je-start, [0, Psi, 0]])

        else:
            assert Strand == "-"
            # print(Annotation[species][Chr]["-"][str(Je)])
            try:
                assert float(Annotation[species][Chr]["-"][str(Js)][1]) > 0 or np.isnan(
                    float(Annotation[species][Chr]["-"][str(Js)][1]))
                assert float(Annotation[species][Chr]["-"][str(Je)][2]) > 0 or np.isnan(
                    float(Annotation[species][Chr]["-"][str(Je)][2]))
            except:
                numNonexp += 1
                logger.info(
                    "Encounter nonexpression sites of {}/{}, will not pass ...".format(numNonexp, NumTotal))
            d["label"].append([Js-start, [0, ref_values, 0]])
            d["label"].append([Je-start, [0, 0, ref_values]])
            d["mutlabel"].append([Js-start, [0, deltaPsi, 0]])
            d["mutlabel"].append([Je-start, [0, 0, deltaPsi]])

            d_vitro["label"].append([Js-start, [0, Psi, 0]])
            d_vitro["label"].append([Je-start, [0, 0, Psi]])

        SaveData.append(d)
        SaveData_vitro.append(d_vitro)
    logger.info("finish handling {}".format(len(SaveData)))
    with open(os.path.join(SavePath, "data.json"), "w") as f:
        f.write(json.dumps(SaveData))
    SaveDataLabeled = refine_with_info(
        {x["name"]: x for x in SaveData}, infoFile)
    with open(os.path.join(SavePath, "data_withlabel.json"), "w") as f:
        f.write(json.dumps(SaveDataLabeled))


if __name__ == "__main__":
    import sys
    if not os.path.exists("data/IKBKAP"):
        os.mkdir("data/IKBKAP")
    if not os.path.exists("data/SMN2"):
        os.mkdir("data/SMN2")
    main(
        "data/ASO/IKBKAP_exon20_exon500bp.ext5000_IVS20.6TC_ASO.targets_dataset.txt",
        "data/ASO/IKBKAP_exon20_exon500bp.ext5000_IVS20.6TC_ASO.targets_seq.txt",
        "data/ASO/IKBKAP_exon20_measured.ASO_combined_info.txt", "hg19",
        "data/IKBKAP", os.path.join(Fapath, "hg19.fa"), 0.265 )

    main(
        "data/ASO/SMN2_exon7_exon500bp.ext5000_ASO.targets_dataset.txt",
        "",
        "data/ASO/SMN2_exon7_measured.ASO_combined_info.txt", "hg19",
        "data/SMN2", os.path.join(Fapath, "hg19.fa"), 0.377 )
