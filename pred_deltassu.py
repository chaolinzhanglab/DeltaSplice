import os
import argparse
import torch
import copy
import pandas as pd
from deltasplice.constant import default_model_paths, model, Fapath, EL, CL,SeqTable, repdict, IN_MAP, default_anno_file
from pyfasta import Fasta
import numpy as np
import json

def main():
    # load config file
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_path", help="the path to input file")
    parser.add_argument("--save_path", help="the path to output file")
    parser.add_argument("--use_reference", help="whether use the default reference information", type=bool, default=False, action="store_true")
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
    # prediction without given exons and ssu
    if "exon_start" not in input_file.keys():
        print("predicting without exon information ...")
        save_file.writelines("chrom,mut_position,strand,pred_ref_acceptor_ssu,pred_ref_donor_ssu,pred_acceptor_deltassu,pred_donor_deltassu\n") 
        for chrom, mut_pos,ref,alt, strand in zip(input_file["chrom"], input_file["mut_position"],input_file["ref"], input_file["alt"], input_file["strand"]):
            pos=mut_pos  
            seq_start=pos-(EL+CL)//2
            seq_end=seq_start+EL+CL
            
            seq=reference_genome[chrom][max(seq_start, 0):min(seq_end, len(reference_genome[chrom]))].upper()
            
            if seq_start<0:
                seq="N"*abs(seq_start)+seq
            if seq_end>len(reference_genome[chrom]):
                seq=seq+"N"*abs(seq_start)
            assert seq[mut_pos-seq_start]==ref.upper(), "the reference allele is not consistant with the reference genome"
            mutseq=seq[:mut_pos-seq_start]+alt.upper()+seq[mut_pos-seq_start+1:]
            
            if strand=="-":
                seq=[repdict[_] for _ in seq][::-1]
                mutseq=[repdict[_] for _ in mutseq][::-1]
                
                
            seq=IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
            mutseq=IN_MAP[[SeqTable[_] for _ in mutseq]][:, :4]
            
            
            d={
                "X":torch.tensor(seq)[None],
                "mutX":torch.tensor(mutseq)[None],
            }
            pred = [m.predict(d, use_ref=True) for m in Models]
            pred_ref = sum([v["single_pred_psi"] for v in pred])/len(pred)
            pred_delta = sum([v["mutY"] for v in pred])/len(pred)-pred_ref
            
            position=np.nonzero(refmat[:, 1:].reshape(-1))
            pred_ref=(pred_ref[:, :, 1:].reshape(-1)[position]).mean()
            pred_delta=(pred_delta[:, :, 1:].reshape(-1)[position]).mean()
            save_file.writelines(f"{chrom},{mut_pos},{ref},{alt},{strand},{jn_end},{jn_start},{ssu},{pred_ref},{pred_delta}\n")
    
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
