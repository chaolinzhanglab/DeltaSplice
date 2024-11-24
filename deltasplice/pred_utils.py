from deltasplice.utils import  write_splice_sites, write_splice_site_file_header
import numpy as np
from loguru import logger
from pyfasta  import Fasta
from deltasplice.constant import default_model_paths, model, IN_MAP, default_anno_file,EL, CL, SeqTable, repdict
import copy
import torch
import json
import pandas as pd
from pkg_resources import resource_filename
from bisect import bisect_left


# annotation
def normalise_chrom(source, target):

    def has_prefix(x):
        return x.startswith('chr')

    if has_prefix(source) and not has_prefix(target):
        return source.strip('chr')
    elif not has_prefix(source) and has_prefix(target):
        return 'chr'+source

    return source

class Annotator:

    def __init__(self, ref_fasta, annotations):

        if annotations == 'grch37':
            annotations = resource_filename(__name__, 'annotations/grch37.txt')
        elif annotations == 'grch38':
            annotations = resource_filename(__name__, 'annotations/grch38.txt')

        df = pd.read_csv(annotations, sep='\t', dtype={'CHROM': object})
        self.genes = df['#NAME'].to_numpy()
        self.chroms = df['CHROM'].to_numpy()
        self.strands = df['STRAND'].to_numpy()
        self.tx_starts = df['TX_START'].to_numpy()+1
        self.tx_ends = df['TX_END'].to_numpy()
        self.exon_starts = [np.asarray([int(i) for i in c.split(',') if i])+1
                            for c in df['EXON_START'].to_numpy()]
        self.exon_ends = [np.asarray([int(i) for i in c.split(',') if i])
                            for c in df['EXON_END'].to_numpy()]
        
        self.ref_fasta = Fasta(ref_fasta)
      

        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        self.models = [copy.deepcopy(model) for _ in default_model_paths]
        [m.load_state_dict(torch.load(resource_filename(__name__, b))) for m,b in zip(self.models, default_model_paths)]
        with open(resource_filename(__name__, default_anno_file), "r") as f:
            default_anno_info=json.load(f)
            SortedKeys = {}
            for sp in default_anno_info:
                SortedKeys[sp] = {}
                for chr in default_anno_info[sp]:
                    SortedKeys[sp][chr] = {}
                    SortedKeys[sp][chr]["+"] = sorted([int(_) for _ in default_anno_info[sp][chr]["+"].keys()])
                    SortedKeys[sp][chr]["-"] = sorted([int(_) for _ in default_anno_info[sp][chr]["-"].keys()])
        self.SortedKeys=SortedKeys
        self.default_anno_info=default_anno_info

    def get_name_and_strand(self, chrom, pos):

        chrom = normalise_chrom(chrom, list(self.chroms)[0])
        idxs = np.intersect1d(np.nonzero(self.chroms == chrom)[0],
                              np.intersect1d(np.nonzero(self.tx_starts <= pos)[0],
                              np.nonzero(pos <= self.tx_ends)[0]))

        if len(idxs) >= 1:
            return self.genes[idxs], self.strands[idxs], idxs
        else:
            return [], [], []

    def get_pos_data(self, idx, pos):

        dist_tx_start = self.tx_starts[idx]-pos
        dist_tx_end = self.tx_ends[idx]-pos
        dist_exon_bdry = min(np.union1d(self.exon_starts[idx], self.exon_ends[idx])-pos, key=abs)
        dist_ann = (dist_tx_start, dist_tx_end, dist_exon_bdry)

        return dist_ann



def get_prediction(inp, models, use_ref, ys=None):
    models[0].predict(inp)
    pred=sum([m.predict(inp, use_ref=use_ref)["single_pred_psi"]
                            for m in models])/len(models)
    if ys is None:
        return pred
    if len(ys.shape)==4:
        ys=ys[:, : , 0]
    assert len(ys.shape)==len(pred.shape)
    ly=ys.shape[1]
    bias=pred.shape[1]-ly
    assert bias%2==0
    pred=pred[:, bias//2:bias//2+ly]
    assert ys.shape==pred.shape
    
    return pred

def get_delta_prediction(record, distance, reference_genome, use_ref, ann, mask):
    cov=2*distance+1
    CL=cov
    wid=10000+cov
    delta_scores=[]
    
    chrom = normalise_chrom(record.chrom, list(ann.ref_fasta.keys())[0])
    pos=record.pos
    ref=record.ref
    pos=pos-1
    (genes, strands, idxs) = ann.get_name_and_strand(record.chrom, record.pos)
    seq_start=pos-(EL+CL)//2
    seq_end=seq_start+EL+CL
    seq= ann.ref_fasta[chrom][max(seq_start, 0):min(seq_end, len(ann.ref_fasta[chrom]))].upper()
    if seq_start<0:
        seq="N"*abs(seq_start)+seq
    if seq_end>len(ann.ref_fasta[chrom]):
        seq=seq+"N"*abs(seq_start)
    
    for j in range(len(record.alts)):
        for i in range(len(idxs)):
            refmat=np.zeros((CL+EL, 3))
    
            if use_ref:
                species=reference_genome.split("/")[-1].replace(".fa", "")
                assert species in ann.SortedKeys, f"{reference_genome} not exists in default reference"
                posidx, startidx, endidx = bisect_left(ann.SortedKeys[species][chrom][strands[i]], pos), bisect_left(ann.SortedKeys[species][chrom][strands[i]], seq_start), bisect_left(ann.SortedKeys[species][chrom][strands[i]], seq_end)
                for v in ann.SortedKeys[species][chrom][strands[i]][startidx:endidx]:
                    refmat[v-seq_start]=ann.default_anno_info[species][chrom][strands[i]][str(v)]
                refmat[np.isnan(refmat)]=1e-3
           
            if '.' in record.alts[j] or '-' in record.alts[j] or '*' in record.alts[j]:
                continue

            if '<' in record.alts[j] or '>' in record.alts[j]:
                continue

            if len(record.ref) > 1 or len(record.alts[j]) > 1:
                delta_scores.append("{}|{}|.|.|.|.|.|.|.|.".format(record.alts[j], genes[i]))
                continue
            assert seq[pos-seq_start]==ref.upper(), f"{ref} {seq[pos-seq_start-1:pos-seq_start+2]} {pos}"
            assert seq[len(seq)//2]==ref.upper(), f"{seq[len(seq)//2]} {ref.upper()} {len(seq)}"
           
            x_ref = seq
            x_alt =seq[:pos-seq_start]+str(record.alts[j]).upper()+seq[pos-seq_start+1:]

            if strands[i]=="+":
                x_ref = IN_MAP[[SeqTable[_] for _ in x_ref]][:, :4]
                x_alt = IN_MAP[[SeqTable[_] for _ in x_alt]][:, :4]
            else:
                x_ref=[repdict[_] for _ in x_ref][::-1]
                x_alt=[repdict[_] for _ in x_alt][::-1]
                x_ref = IN_MAP[[SeqTable[_] for _ in x_ref]][:, :4]
                x_alt = IN_MAP[[SeqTable[_] for _ in x_alt]][:, :4]
                refmat=refmat[::-1]
            refmat[:, 0]=1-refmat[:, 1:].sum(-1)
            refmat=refmat[EL//2:EL//2+CL].copy()
    
            if use_ref:
                d={
                    "X":torch.tensor(x_ref)[None],
                    "mutX":torch.tensor(x_alt)[None],
                    "single_pred_psi":torch.tensor(refmat)[None]
                }
            else:
                d={
                    "X":torch.tensor(x_ref)[None],
                    "mutX":torch.tensor(x_alt)[None],

                }
            pred = [m.predict(d, use_ref=use_ref) for m in ann.models]
            pred_ref = sum([v["single_pred_psi"] for v in pred])/len(pred)
            pred_alt = sum([v["mutY"] for v in pred])/len(pred)
    
            assert cov==pred_ref.shape[1], f"{cov} {pred_ref.shape}"
            if strands[i] == '-':
                pred_ref = pred_ref[:, ::-1]
                pred_alt = pred_alt[:, ::-1]
            ref_len = len(record.ref)
            alt_len = len(record.alts[j])
            del_len = max(ref_len-alt_len, 0)
            dist_ann = ann.get_pos_data(idxs[i], record.pos)
           
            if ref_len > 1 and alt_len == 1:
                y_alt = np.concatenate([
                    y_alt[:, :cov//2+alt_len],
                    np.zeros((1, del_len, 3)),
                    y_alt[:, cov//2+alt_len:]],
                    axis=1)
            elif ref_len == 1 and alt_len > 1:
                y_alt = np.concatenate([
                    y_alt[:, :cov//2],
                    np.max(y_alt[:, cov//2:cov//2+alt_len], axis=1)[:, None, :],
                    y_alt[:, cov//2+alt_len:]],
                    axis=1)

            y = np.concatenate([pred_ref, pred_alt], 0)

            idx_pa = (y[1, :, 1]-y[0, :, 1]).argmax()
            idx_na = (y[0, :, 1]-y[1, :, 1]).argmax()
            idx_pd = (y[1, :, 2]-y[0, :, 2]).argmax()
            idx_nd = (y[0, :, 2]-y[1, :, 2]).argmax()
           
            mask_pa = np.logical_and((idx_pa-cov//2 == dist_ann[2]), mask)
            mask_na = np.logical_and((idx_na-cov//2 != dist_ann[2]), mask)
            mask_pd = np.logical_and((idx_pd-cov//2 == dist_ann[2]), mask)
            mask_nd = np.logical_and((idx_nd-cov//2 != dist_ann[2]), mask)

            delta_scores.append("{}|{}|{:.4f}|{:.4f}|{:.4f}|{:.4f}|{}|{}|{}|{}".format(
                                record.alts[j],
                                genes[i],
                                (y[1, idx_pa, 1]-y[0, idx_pa, 1])*(1-mask_pa),
                                (y[0, idx_na, 1]-y[1, idx_na, 1])*(1-mask_na),
                                (y[1, idx_pd, 2]-y[0, idx_pd, 2])*(1-mask_pd),
                                (y[0, idx_nd, 2]-y[1, idx_nd, 2])*(1-mask_nd),
                                idx_pa-cov//2,
                                idx_na-cov//2,
                                idx_pd-cov//2,
                                idx_nd-cov//2))

    return delta_scores

def eval_test_data(data, models, save_path, use_ref):
    save_paths=(open(save_path+"_acceptor", "w"),open(save_path+"_donor", "w"))
    [write_splice_site_file_header(_) for _ in save_paths]
    Y_true_1 = []
    Y_true_2 = []
    Y_pred_1 = []
    Y_pred_2 = []
    for i, d in enumerate(data):
        if i%100==0:
            print("finish eval ", i)
        X, Y = d["X"], d["single_pred_psi"].numpy()
        pred=get_prediction(d, models, False, Y)
        CHROM = d["chrom"]
        NAME = d["name"]
        STRAND = d["strand"]
        TX_START = d["txstart"].numpy()
        TX_END = d["txend"].numpy()
        SPECIES = d["species"]
        write_splice_sites(save_paths, CHROM, NAME, STRAND,
                                    TX_START, TX_END, SPECIES, pred, Y)
        
        pred_sites = np.nonzero(Y[:, :, 1:].sum(-1))

        Y_true_1.extend(Y[pred_sites][:, 1].flatten())
        Y_true_2.extend(Y[pred_sites][:, 2].flatten())
        Y_pred_1.extend(pred[pred_sites][:, 1].flatten())
        Y_pred_2.extend(pred[pred_sites][:, 2].flatten())
        
    [_.close() for _ in save_paths]
    
def eval_mut_data(data, models, save_path, use_ref):
    Pred_ref = []
    Pred_delta = []
    Gt_delta = []
    for i, d in enumerate(data):
        if i%100==0:
            print("finish eval ", i)
        gt_delta = d["mutY"]  # GT of delta usage
        pred = [m.predict(d, use_ref=use_ref) for m in models]
        pred_ref = sum([v["single_pred_psi"] for v in pred])/len(pred)
        pred_delta = sum([v["mutY"] for v in pred])/len(pred)-pred_ref

        if len(gt_delta.shape) == 4:
            gt_delta = gt_delta[:, :, 0]
        assert len(pred_delta.shape) == 3
        assert len(gt_delta.shape) == 3

        gt_delta = gt_delta[:, :, 1:].reshape(-1)
        pred_ref = pred_ref[:, :, 1:].reshape(-1)
        pred_delta = pred_delta[:, :, 1:].reshape(-1)
        position = np.nonzero(gt_delta)
        if len(position) == 0:
            logger.warning(
                "Encounting mutations without non-zero labels, will be skipped")
            continue
        gt_delta = gt_delta[position].mean()
        pred_ref = pred_ref[position].mean()
        pred_delta = pred_delta[position].mean()

        Pred_delta.append(pred_delta)
        Gt_delta.append(gt_delta)
        Pred_ref.append(pred_ref)
   
    with open(save_path, "w") as f:
        f.writelines("Pred_ref\tPred_delta\tGt_delta\n")
        for a, b, c in zip(Pred_delta, Gt_delta, Pred_ref):
            f.writelines("{}\t{}\t{}\n".format(c, a, b))