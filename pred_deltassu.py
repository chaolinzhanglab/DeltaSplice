import os
import argparse
import torch
import copy
import pandas as pd
from deltasplice.constant import default_model_paths, model, Fapath, EL, CL,SeqTable, repdict, IN_MAP, default_anno_file
from pyfasta import Fasta
import numpy as np
import json
from bisect import bisect_left

def main():
    # load config file
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_path", help="the path to input file")
    parser.add_argument("--save_path", help="the path to output file")
    parser.add_argument("--use_reference", help="whether use the default reference information", default=False, action="store_true")
    parser.add_argument("--window_size", type=int, default=200, help="when exon is not given, the size of the window around the predicted mutation.")
    parser.add_argument("--genome", help="reference genome")
    args = parser.parse_args()
    Models = [copy.deepcopy(model) for _ in default_model_paths]
    [m.load_state_dict(torch.load(b)) for m,b in zip(Models, default_model_paths)]
    
    input_file=pd.read_csv(args.data_path)
    reference_genome=Fasta(os.path.join(Fapath, args.genome+".fa"))
    save_file=open(args.save_path, "w")
    if args.use_reference:
        print("using default reference information ...")
        with open(default_anno_file, "r") as f:
            default_anno_info=json.load(f)
            SortedKeys = {}
            for sp in default_anno_info:
                SortedKeys[sp] = {}
                for chr in default_anno_info[sp]:
                    SortedKeys[sp][chr] = {}
                    SortedKeys[sp][chr]["+"] = sorted([int(_) for _ in default_anno_info[sp][chr]["+"].keys()])
                    SortedKeys[sp][chr]["-"] = sorted([int(_) for _ in default_anno_info[sp][chr]["-"].keys()])
    # prediction without given exons and ssu
    if "exon_start" not in input_file.keys():
        print("predicting without exon information ...")
        save_file.writelines("chrom,mut_position,strand,position,reference_acceptor_ssu,reference_donor_ssu,pred_ref_acceptor_ssu,pred_ref_donor_ssu,pred_acceptor_deltassu,pred_donor_deltassu\n") 
        for chrom, mut_pos,ref,alt, strand in zip(input_file["chrom"], input_file["mut_position"],input_file["ref"], input_file["alt"], input_file["strand"]):
            pos=mut_pos  
            seq_start=pos-(EL+CL)//2
            seq_end=seq_start+EL+CL
            
            seq=reference_genome[chrom][max(seq_start, 0):min(seq_end, len(reference_genome[chrom]))].upper()
            
            if seq_start<0:
                seq="N"*abs(seq_start)+seq
            if seq_end>len(reference_genome[chrom]):
                seq=seq+"N"*abs(seq_start)
            assert seq[mut_pos-seq_start]==ref.upper()
            assert seq[len(seq)//2]==ref.upper()
            mutseq=seq[:mut_pos-seq_start]+alt.upper()+seq[mut_pos-seq_start+1:]
            refmat=np.zeros((CL+EL, 3))
            if args.use_reference:
                species=args.genome.split("/")[-1].replace(".fa", "")
                assert species in SortedKeys, f"{args.genome} not exists in default reference"
                posidx, startidx, endidx = bisect_left(SortedKeys[species][chrom][strand], pos), bisect_left(SortedKeys[species][chrom][strand], seq_start), bisect_left(SortedKeys[species][chrom][strand], seq_end)
                for v in SortedKeys[species][chrom][strand][startidx:endidx]:
                    refmat[v-seq_start]=default_anno_info[species][chrom][strand][str(v)]
                refmat[np.isnan(refmat)]=1e-3
                    
            if strand=="-":
                seq=[repdict[_] for _ in seq][::-1]
                mutseq=[repdict[_] for _ in mutseq][::-1]
                refmat=refmat[::-1]
                
            seq=IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
            mutseq=IN_MAP[[SeqTable[_] for _ in mutseq]][:, :4]
            refmat[:, 0]=1-refmat[:, 1:].sum(-1)
            
            refmat=refmat[EL//2:EL//2+CL].copy()
            
            if args.use_reference:
                d={
                    "X":torch.tensor(seq)[None],
                    "mutX":torch.tensor(mutseq)[None],
                    "single_pred_psi":torch.tensor(refmat)[None]
                }
                use_ref=False
            else:
                d={
                    "X":torch.tensor(seq)[None],
                    "mutX":torch.tensor(mutseq)[None],
                }
                use_ref=True
            pred = [m.predict(d, use_ref=use_ref) for m in Models]
            pred_ref = sum([v["single_pred_psi"] for v in pred])/len(pred)
            pred_delta = sum([v["mutY"] for v in pred])/len(pred)-pred_ref
            refmat=refmat[None]
            assert refmat.shape==pred_ref.shape, f"{refmat.shape}{pred_ref.shape}"
            if strand=="-":
                pred_ref=np.flip(pred_ref, 1)
                pred_delta=np.flip(pred_delta, 1)
                refmat=np.flip(refmat, 1)
            
            write_window_start=pred_ref.shape[1]//2-args.window_size//2
            write_window_end=pred_ref.shape[1]//2+args.window_size//2
            
            for i in range(write_window_start, write_window_end+1):
                acceptor_ref=refmat[0, i, 1]
                donor_ref=refmat[0, i, 2]
                acceptor_ref_pred=pred_ref[0, i, 1]
                donor_ref_pred=pred_ref[0, i, 2]
                acceptor_delta=pred_delta[0, i, 1]
                donor_delta=pred_delta[0, i, 2]
                position=pos-args.window_size//2+i-write_window_start
                
                save_file.writelines(f"{chrom},{mut_pos},{strand},{position},{acceptor_ref},{donor_ref},{acceptor_ref_pred},{donor_ref_pred},{acceptor_delta},{donor_delta}\n")
        save_file.close()
    else:
        save_file.writelines("chrom,mut_position,strand,exon_start,exon_end,ssu,pred_ref,pred_deltassu\n")
        for chrom, mut_pos,ref,alt, strand,jn_start, jn_end, ssu in zip(input_file["chrom"], input_file["mut_position"],input_file["ref"], input_file["alt"], input_file["strand"], input_file["exon_end"], input_file["exon_start"], input_file["ssu"]):
            pos=(jn_start+jn_end)//2        
            seq_start=pos-(EL+CL)//2
            seq_end=seq_start+EL+CL
            ssu=max(ssu, 1e-3)
            
            seq=reference_genome[chrom][max(seq_start, 0):min(seq_end, len(reference_genome[chrom]))].upper()
            
            if seq_start<0:
                seq="N"*abs(seq_start)+seq
            if seq_end>len(reference_genome[chrom]):
                seq=seq+"N"*abs(seq_start)
            assert seq[mut_pos-seq_start]==ref.upper()
            mutseq=seq[:mut_pos-seq_start]+alt.upper()+seq[mut_pos-seq_start+1:]
            refmat=np.zeros((CL+EL, 3))
            if args.use_reference:
                #assert np.isnan(ssu), "when set use_reference, ssu must be nan"
                species=args.genome.split("/")[-1].replace(".fa", "")
                assert species in SortedKeys, f"{args.genome} not exists in default reference"
                if str(jn_end) in default_anno_info[species][chrom][strand]:
                    ssu_end=sum( default_anno_info[species][chrom][strand][str(jn_end)])
                    if np.isnan(ssu_end):
                        ssu_end=None
                else:
                    print(f"{jn_end} not in defalt reference ssu")
                    ssu_end=None
                    
                if str(jn_start) in default_anno_info[species][chrom][strand]:
                    ssu_start= sum(default_anno_info[species][chrom][strand][str(jn_start)])
                    if np.isnan(ssu_start):
                        ssu_start=None
                else:
                    print(f"{jn_start} not in defalt reference ssu")
                    ssu_start=None
                if ssu_end is not None and ssu_start is not None:
                    ssu=(ssu_start+ssu_end)/2
                elif ssu_end is not None:
                    ssu=ssu_end
                elif ssu_start is not None:
                    ssu=ssu_start
                else:
                    print(f"exon not found, 1e-3 will be used as ssu")
                    ssu=1e-3
                ssu=max(ssu, 1e-3)
                
            if strand=="-":
                seq=[repdict[_] for _ in seq][::-1]
                mutseq=[repdict[_] for _ in mutseq][::-1]
                
                refmat[jn_end-seq_start, 2]=ssu
                refmat[jn_start-seq_start, 1]=ssu
                refmat=refmat[::-1]
            else:
                refmat[jn_start-seq_start, 2]=ssu
                refmat[jn_end-seq_start, 1]=ssu
                
                
            seq=IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
            mutseq=IN_MAP[[SeqTable[_] for _ in mutseq]][:, :4]
            refmat[:, 0]=1-refmat[:, 1:].sum(-1)
            
            refmat=refmat[EL//2:EL//2+CL].copy()
            
            
            d={
                "X":torch.tensor(seq)[None],
                "mutX":torch.tensor(mutseq)[None],
                "single_pred_psi":torch.tensor(refmat)[None]
            }
            pred = [m.predict(d, use_ref=True) for m in Models]
            pred_ref = sum([v["single_pred_psi"] for v in pred])/len(pred)
            pred_delta = sum([v["mutY"] for v in pred])/len(pred)-pred_ref
            
            position=np.nonzero(refmat[:, 1:].reshape(-1))
            pred_ref=(pred_ref[:, :, 1:].reshape(-1)[position]).mean()
            pred_delta=(pred_delta[:, :, 1:].reshape(-1)[position]).mean()
            save_file.writelines(f"{chrom},{mut_pos},{ref},{alt},{strand},{jn_end},{jn_start},{ssu},{pred_ref},{pred_delta}\n")
        save_file.close()
            
        
    

if __name__ == "__main__":
    main()
