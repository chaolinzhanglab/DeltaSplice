from loguru import logger
import os
import sys
import json
from config import CL, AnnoPath,  EL, Fapath
from pyfasta import Fasta
import pandas
SavePath = "data/FAS"
# hgfile["chr10"][90770509:90770572].upper() matches


def main(File, species, SavePath):
    mutinfo = pandas.read_csv(File)
    fafile = Fasta(os.path.join(Fapath, species+".fa"))

    padsize = (EL+CL)//2
    NumTotal = 0
    Strand = "+"

    seq_start, seq_end = central_pos-padsize, central_pos+padsize
    oriseq = fafile[Chr][seq_start:seq_end].upper()
    oriseq = (oriseq[:exon_start-seq_start]+"GATCCAGATCTAACTTGCTGTGGTTGTGTCTCCTGCTTCTCCCGATTCTAGTAATTGTTTGGG" +
              oriseq[exon_end-seq_start:]).upper()
    assert len(oriseq) == EL+CL

    SaveData = []
    for idx, line in mutinfo.iterrows():
        NumTotal += 1
        Name = line["Mutation.IDs"]
        if not isinstance(Name, str):
            continue
        s = line["Sequence"]
        assert len(s) == exon_end-exon_start

        mutseq = (oriseq[:exon_start-seq_start]+s +
                  oriseq[exon_end-seq_start:]).upper()
        mut_num = sum([1 if a != b else 0 for a, b in zip(oriseq, mutseq)])
        assert mut_num == len(Name.split(
            ";")), f'{mut_num} {len(Name.split(";"))} {len(mutseq)} {len(oriseq)}'
        Psi = line["Mean"]/100.
        deltaPsi = Psi-0.96
        if deltaPsi == 0:
            deltaPsi = deltaPsi+1e-5

        d = {"species": species, "chrom": Chr, "start": seq_start, "end": seq_end, "strand": Strand,
             "mutseq": mutseq, "seq": oriseq, "label": [], "mutlabel": [], "name": Name}

        d["label"].append([exon_end-seq_start-1, [0, 0, 0.96]])
        d["label"].append([exon_start-seq_start, [0, 0.96, 0]])
        d["mutlabel"].append([exon_end-seq_start-1, [0, 0, deltaPsi]])
        d["mutlabel"].append([exon_start-seq_start, [0, deltaPsi, 0]])
        SaveData.append(d)

    logger.info("finish handling {}".format(len(SaveData)))
    with open(os.path.join(SavePath, "data.json"), "w") as f:
        f.write(json.dumps(SaveData))


if __name__ == "__main__":
    import sys
    File = "data/fas.csv"
    Chr = "chr10"
    exon_start = 90770509
    exon_end = 90770572
    central_pos = (exon_start+exon_end)//2
    if not os.path.exists(SavePath):
        os.mkdir(SavePath)

    main(File, "hg19", SavePath)
