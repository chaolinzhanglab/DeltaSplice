from pyliftover import LiftOver
from functools import partial
from loguru import logger
import json
import os
from pyfasta import Fasta
from config import Fapath

TABLE_FILE = "data/gene_dataset.tsu.txt"
AnnoPath = "data/annotation/"
SavePath = "data/Hg19VsOthers"
hg19fa = Fasta(Fapath+"/hg19.fa")
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
Test_Chromes = ["chr3", "chr5", "chr9", "chr7", "chr1"]
Species = ["susScr11", "mm10", "rheMac10", "rn6", "panTro5", "bosTau9"]
ContextLength = 35000//2
if not os.path.exists(SavePath):
    os.mkdir(SavePath)


def parse_bed_line(line):
    eps = 1e-3
    line = line.replace("\n", "").split("\t")
    fa, hgchrom, name, _, chrom, strand, start, end, ss3s, ss5s, ss3vs, ss5vs = line
    ss3s = [int(x)-1 for x in ss3s.split(",") if len(x) > 0]  # with 0 base
    ss5s = [int(x)-1 for x in ss5s.split(",") if len(x) > 0]
    ss3vs = [max(float(x), eps) for x in ss3vs.split(",") if len(x) > 0]
    ss5vs = [max(float(x), eps) for x in ss5vs.split(",") if len(x) > 0]
    return fa, chrom, name, strand, int(start)-1, int(end)-1, ss3s, ss5s, ss3vs, ss5vs, hgchrom


def load_annotation():
    with open(os.path.join(AnnoPath, "data.json"), "r") as f:
        Annotation = json.load(f)
    return Annotation


def load_table(startspecies, endspecies, LO_FUNC, startfa, endfa):
    with open(TABLE_FILE, "r") as f:
        content = f.readlines()
    Save = []
    prefa = "chr1"
    LO = LO_FUNC
    Annotation = load_annotation()
    start = 35000//2
    RANGE = 200
    logger.info("start handling {}".format(endspecies))
    for line in content:
        if line.startswith(startspecies):
            fa, chrom, name, strand, _, _, ss3s, ss5s, ss3vs, ss5vs, _ = parse_bed_line(
                line)
            if chrom != prefa:
                print(prefa, chrom, len(Save))
                prefa = chrom
                LO = LO_FUNC
            if chrom not in Test_Chromes:
                continue
            for i, v in zip(ss3s+ss5s, ss3vs+ss5vs):
                cidx = LO.convert_coordinate(chrom, i, strand)
                v = float(v)
                if len(cidx) > 0 and v > 0:
                    for t in cidx:
                        c, ti, s, _ = t
                        if c in Annotation[endspecies] and str(ti) in Annotation[endspecies][c][s]:
                            tv = sum(Annotation[endspecies][c][s][str(ti)])
                            if tv > 0 and i+start+1 < len(startfa[chrom]) and i-start >= 0 and ti-start >= 0 and ti+start+1 < len(endfa[c]):
                                save = {"{}_seq".format(startspecies): startfa[chrom][i-start:i+start+1], "{}_seq".format(endspecies): endfa[c][ti-start:ti+start+1],
                                        "{}_strand".format(startspecies): strand, "{}_strand".format(endspecies): s,
                                        "{}_label".format(startspecies): Annotation[startspecies][chrom][strand][str(i)], "{}_label".format(endspecies): Annotation[endspecies][c][s][str(ti)],
                                        "{}_name".format(startspecies): name, "{}_chrom".format(startspecies): chrom, "{}_idx".format(startspecies): i, "{}_chrom".format(endspecies): c,"{}_idx".format(endspecies): ti,}
                                assert len(
                                    save["{}_seq".format(startspecies)]) > 0
                                assert len(
                                    save["{}_seq".format(endspecies)]) > 0

                                transitions = []

                                hseq = startfa[chrom][i -
                                                      RANGE:i+RANGE+1].upper()
                                pseq = endfa[c][ti-RANGE:ti+RANGE+1].upper()
                                endseq = endfa[c][ti-start:ti+start+1]
                                if strand != s:
                                    pseq = "".join([repdict[_]
                                                   for _ in pseq][::-1])
                                    endseq = "".join(
                                        [repdict[_] for _ in endseq.upper()][::-1])

                                diff_num = sum(
                                    [1 for a, b in zip(hseq, pseq) if a != b])
                                preseq = startfa[chrom][i-start:i+start+1]
                                transitions.append(preseq)
                                if diff_num > 0 and abs(v-tv) > 0.8:
                                    for idx in range(i-RANGE, i+RANGE+1):
                                        if hseq[idx-i+RANGE] != pseq[idx-i+RANGE]:
                                            transitions.append(
                                                preseq[:idx-i+start]+pseq[idx-i+RANGE]+preseq[idx-i+start+1:])
                                    print(len(transitions), diff_num)
                                transitions.append(endseq)
                                save["transition"] = transitions
                                Save.append(save)
    return Save


def analysis(data):
    for d in data:
        hg19_seq = d["hg19_seq"]
        panTro5_seq = d["panTro5_seq"]
        assert len(hg19_seq) == len(panTro5_seq)
        num = 0
        for a, b in zip(hg19_seq, panTro5_seq):
            if a != b:
                num += 1
        print(d["hg19_name"], num, d["hg19_label"], d["panTro5_label"],
              d["hg19_strand"], d["panTro5_strand"])


for s in Species:
    LO_FUNC = LiftOver(
        "data/Chains/hg19To{}.over.chain".format(s[0].upper()+s[1:]))
    print("finish loading")
    fa = Fasta(Fapath+"/{}.fa".format(s))
    Save = load_table("hg19", s, LO_FUNC, hg19fa, fa)
    logger.info("Num of Matched Gene is {}".format(len(Save)))
    with open(os.path.join(SavePath, "{}.json".format(s)), "w") as f:
        f.write(json.dumps(Save))
