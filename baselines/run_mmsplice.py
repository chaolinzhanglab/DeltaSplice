# Import
from mmsplice import MMSplice, LINEAR_MODEL
import sys
from mmsplice.exon_dataloader import ExonDataset
from kipoi.data import DataLoader
from kipoi.data_utils import numpy_collate
import numpy as np
import pandas
import json
from pyfasta import Fasta
Modelpath = "baselines/MMSplice_paper/data/vexseq/scale_model.pkl"
repdict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", '-': '-'}
# example files
vcf, fasta, csv = sys.argv[1:4]
fafile = Fasta(fasta)
temp_vcf_file = "tmp_vcf/{}".format(vcf.replace("/", "_"))
labels = []
if vcf.endswith("json"):
    with open(vcf, "r") as f:
        content = json.load(f)
    with open(temp_vcf_file, "w") as f:
        f.writelines(
            "ID,seqnames,start,end,width,strand,hg19_variant_position,reference,variant\n")
        for d in content:
            if "mutpos" in d:
                species, chrom, start, end, strand, mutpos, oriv, mutv = d["species"], d["chrom"], int(
                    d["start"]), int(d["end"]), d["strand"], int(d["mutpos"]), d["oriv"].upper(), d["mutv"].upper()
                mutpos = mutpos+start+1
                assert fafile[chrom][mutpos-1].upper() == oriv

            else:
                chrom, start, end, strand, seq, mutseq = d["chrom"], int(
                    d["start"]), int(d["end"]), d["strand"], d["seq"].upper(), d["mutseq"].upper()
                try:
                    species = d["species"]
                except:
                    species = "hg19"
                for i in range(len(seq)):
                    if seq[i] != mutseq[i]:
                        mutpos = start+i+1

                        oriv = seq[i].upper()
                        mutv = mutseq[i].upper()
                        break
                else:
                    print(seq)
                    print(mutseq)
                    mutpos = start
                    oriv = seq[0]
                    mutv = mutseq[0]
            # print(fafile[chrom][mutpos-3:mutpos+3].upper(), oriv, mutpos)

            name = "{}:{}:{}>{}".format(chrom,mutpos, oriv, mutv)
            l1, l2 = d["label"]
            l1, l2 = l1[0]+start, l2[0]+start
            f.writelines("{},{},{},{},{},{},{},{},{}\n".format(name, chrom, min(
                l1, l2)+1, max(l1, l2)+1, max(l1, l2)-min(l1, l2), strand.replace(" ", ""), mutpos, oriv, mutv))
    dl = ExonDataset(temp_vcf_file,
                     fasta, split_seq=False, overhang=(100, 100))
    dl = DataLoader(dl, batch_size=dl.__len__(),
                    collate_fn=numpy_collate, shuffle=False)
    dt = next(iter(dl))

    csvf = pandas.read_csv(temp_vcf_file, sep=',')

    assert len(dt["inputs"]["seq"]) == len(content)

    for i in range(len(content)):
        d = content[i]
        labels.append(sum([sum(_[1]) for _ in d["label"]])/len(d["label"]))

        if "mutpos" in d:
            species, chrom, start, end, strand, mutpos, oriv, mutv = d["species"], d["chrom"], int(
                d["start"]), int(d["end"]), d["strand"], int(d["mutpos"]), d["oriv"].upper(), d["mutv"].upper()
            seq = fafile[chrom][start:end].upper()
            mutseq = (seq[:mutpos]+mutv+seq[mutpos+1:]).upper()
            templateseq = seq
            assert seq[mutpos] == oriv.upper()
        else:
            species, chrom, start, end, strand, seq, mutseq = d["species"], d["chrom"], int(
                d["start"]), int(d["end"]), d["strand"], d["seq"].upper(), d["mutseq"].upper()
            templateseq = fafile[chrom][start:end].upper()
        if strand == "-":
            seq = "".join([repdict[_] for _ in seq[::-1]])
            mutseq = "".join([repdict[_] for _ in mutseq[::-1]])
            templateseq = "".join([repdict[_] for _ in templateseq[::-1]])
        idx = templateseq.find(dt["inputs"]["seq"][i])

        assert idx > -1
        length = len(dt["inputs"]["seq"][i])

        # assert sum([1 if a != b else 0 for a, b in zip(
        #    seq[idx:idx+length], mutseq[idx:idx+length])]) == sum([1 if a != b else 0 for a, b in zip(
        #        seq, mutseq)])

        dt["inputs"]["seq"][i] = seq[idx:idx+length]
        dt["mut_inputs"]["seq"][i] = mutseq[idx:idx+length]
else:
    dl = ExonDataset(vcf,
                     fasta, split_seq=False, overhang=(100, 100))
    csvf = pandas.read_csv(vcf, sep=',')
    dl = DataLoader(dl, batch_size=dl.__len__(),
                    collate_fn=numpy_collate, shuffle=False)
    dt = next(iter(dl))

model = MMSplice(
    exon_cut_l=0,
    exon_cut_r=0,
    acceptor_intron_cut=6,
    donor_intron_cut=6,
    acceptor_intron_len=50,
    acceptor_exon_len=3,
    donor_exon_len=5,
    donor_intron_len=13)
huber = LINEAR_MODEL

print("start ref pred")
ref_pred = model.predict_on_unsplitted_batch(dt['inputs'])
print("start mut pred")
alt_pred = model.predict_on_unsplitted_batch(dt['mut_inputs'])

X = alt_pred-ref_pred

exon_overlap = np.logical_or(np.logical_and(
    X[:, 1] != 0, X[:, 2] != 0), np.logical_and(X[:, 2] != 0, X[:, 3] != 0))

acceptor_intron_overlap = np.logical_and(X[:, 0] != 0, X[:, 1] != 0)
donor_intron_overlap = np.logical_and(X[:, 3] != 0, X[:, 4] != 0)

X = np.hstack([X, (X[:, 2]*exon_overlap).reshape(-1, 1)])
X = np.hstack([X, (X[:, 4]*donor_intron_overlap).reshape(-1, 1)])
X = np.hstack([X, (X[:, 0]*acceptor_intron_overlap).reshape(-1, 1)])

delt_measured = huber.predict(X)
assert len(csvf["ID"]) == len(delt_measured)
eps = 1e-10


def expit(x):
    return 1. / (1. + np.exp(-x))


def clip(x, clip_threshold=eps):
    return np.clip(x, clip_threshold, 1 - clip_threshold)


def logit(x, clip_threshold=eps):
    x = clip(x, clip_threshold=clip_threshold)
    return np.log(x) - np.log(1 - x)


with open(csv, "w") as f:
    for i, a, b in zip(range(len(delt_measured)), csvf["ID"], delt_measured):
        if len(labels) > 0 and not np.isnan(sum(labels)):
            b = expit(b+logit(labels[i]))-labels[i]
        f.writelines("{}\t{}\n".format(a, b))
