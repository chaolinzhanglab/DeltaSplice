import os
import numpy as np
from pyfasta import Fasta
import json
from loguru import logger
from config import Fapath, AnnoPath, NpyPath,  CL, EL, TableFile, Train_Chromes, Valid_Chromes, Test_Chromes
from utils import parse_bed_line
import random
from bisect import bisect_left


def main(Data, Path):
    logger.info("handling {} lines for {}".format(len(Data), Path))
    Fafiles = {}
    SaveData = []
    padsize = EL//2
    PASTSITES = set()

    with open(os.path.join(AnnoPath, "data.json"), "r") as f:
        Annotation = json.load(f)
    SortedKeys = {}
    for sp in Annotation:
        SortedKeys[sp] = {}
        for chr in Annotation[sp]:
            SortedKeys[sp][chr] = {}
            SortedKeys[sp][chr]["+"] = sorted([int(_) for _ in Annotation[sp][chr]["+"].keys()])
            SortedKeys[sp][chr]["-"] = sorted([int(_) for _ in Annotation[sp][chr]["-"].keys()])

    for idx, line in enumerate(Data):
        species, chrom, pos, strand, name = line
        key = "{}_{}_{}".format(chrom, pos, strand)
        if key not in PASTSITES:
            PASTSITES.add(key)
        else:
            continue
        # prepare data
        if species not in Fafiles:
            logger.info("loading fa file of {}".format(species))
            if len(Fafiles) > 3:
                Fafiles.pop(random.choice(list(Fafiles.keys())))
            Fafiles[species] = Fasta(os.path.join(Fapath, species+".fa"))
        # write

        start, end = pos-padsize, pos+padsize+CL
        end = min(end, len(Fafiles[species][chrom])-1)
        start = end-CL-EL
        if start < 0 or end >= len(Fafiles[species][chrom]):
            print("exceed", start, end, len(Fafiles[species][chrom]))
            continue
        d = {"species": species, "chrom": chrom, "start": start, "end": end, "strand": strand, "label": [], "name": name}

        posidx, startidx, endidx = bisect_left(SortedKeys[species][chrom][strand], pos), bisect_left(SortedKeys[species][chrom][strand], start), bisect_left(SortedKeys[species][chrom][strand], end)
        assert posidx < len(SortedKeys[species][chrom][strand]) and SortedKeys[species][chrom][strand][posidx] == pos
        assert str(pos) in Annotation[species][chrom][strand]
        for v in SortedKeys[species][chrom][strand][startidx:endidx]:
            d["label"].append([v-start, Annotation[species][chrom][strand][str(v)]])
            if v-start >= EL//2 and v-start < end-start-EL//2:
                PASTSITES.add("{}_{}_{}".format(chrom, v, strand))

        SaveData.append(d)
    logger.info("finish handling {}".format(len(SaveData)))
    with open(os.path.join(Path, "data.json"), "w") as f:
        f.write(json.dumps(SaveData))
    with open(os.path.join(Path, "human.json"), "w") as f:
        f.write(json.dumps([x for x in SaveData if x["species"] == "hg19"]))


if __name__ == "__main__":
    if not os.path.exists(NpyPath):
        os.mkdir(NpyPath)
    trainpath, validpath, testpath = os.path.join(NpyPath, "train"), os.path.join(NpyPath, "valid"), os.path.join(NpyPath, "test")
    if not os.path.exists(trainpath):
        os.mkdir(trainpath)
    if not os.path.exists(validpath):
        os.mkdir(validpath)
    if not os.path.exists(testpath):
        os.mkdir(testpath)

    Trainidxs, Valididxs, Testidxs = [], [], []
    loginfo={"train":{}, "test":{}, "valid":{}}
    with open(TableFile, "r") as f:
        expinfo = f.readlines()
    parsedexpinfo = [parse_bed_line(_) for _ in expinfo]
    def add_line_to_loginfo(line, group):
        fa, chrom, name, strand, _, _, ss3s, ss5s, ss3vs, ss5vs, hgchrom = line
        if fa not in loginfo[group]:
            loginfo[group][fa]={"gene":set(), "num_non_nan_3":0, "num_non_nan_5":0, "num_nan_3":0, "num_nan_5":0} 
        loginfo[group][fa]["gene"].add(name)
        for idx, value in zip(ss3s, ss3vs):
            if np.isnan(value) :
                loginfo[group][fa]["num_nan_3"]+=1
            else:
                loginfo[group][fa]["num_non_nan_3"]+=1
          
        for idx, value in zip(ss5s, ss5vs):
            if np.isnan(value) :
                loginfo[group][fa]["num_nan_5"]+=1
            else:
                loginfo[group][fa]["num_non_nan_5"]+=1
            
    for line in parsedexpinfo:
        fa, chrom, name, strand, _, _, ss3s, ss5s, ss3vs, ss5vs, hgchrom = line
        if fa == "hg19":
            if chrom in Train_Chromes:
                add_line_to_loginfo(line, "train")
                for idx, value in zip(ss3s+ss5s, ss3vs+ss5vs):
                    Trainidxs.append([fa, chrom, idx, strand, name])
            elif chrom in Valid_Chromes:
                add_line_to_loginfo(line, "valid")
                for idx, value in zip(ss3s+ss5s, ss3vs+ss5vs):
                    Valididxs.append([fa, chrom, idx, strand, name])
            else:
                assert chrom in Test_Chromes
                add_line_to_loginfo(line, "test")
                for idx, value in zip(ss3s+ss5s, ss3vs+ss5vs):
                    Testidxs.append([fa, chrom, idx, strand, name])
        else:
            if hgchrom in Train_Chromes:
                add_line_to_loginfo(line, "train")
                for idx, value in zip(ss3s+ss5s, ss3vs+ss5vs):
                    Trainidxs.append([fa, chrom, idx, strand, name])
            elif hgchrom in Valid_Chromes:
                add_line_to_loginfo(line, "valid")
                for idx, value in zip(ss3s+ss5s, ss3vs+ss5vs):
                    Valididxs.append([fa, chrom, idx, strand, name])
            elif hgchrom in Test_Chromes:
                add_line_to_loginfo(line, "test")
                for idx, value in zip(ss3s+ss5s, ss3vs+ss5vs):
                    Testidxs.append([fa, chrom, idx, strand, name])
            else:
                print(line)
    main(Trainidxs, trainpath)
    main(Valididxs, validpath)
    main(Testidxs, testpath)
    for group in loginfo:
        with open(group, "w") as f:
            for species in loginfo[group]:
                f.writelines(",".join([species, str(len(loginfo[group][species]["gene"])), str(loginfo[group][species]["num_nan_3"]), str(loginfo[group][species]["num_non_nan_3"]),str(loginfo[group][species]["num_nan_5"]), str(loginfo[group][species]["num_non_nan_5"])])+"\n")
            
