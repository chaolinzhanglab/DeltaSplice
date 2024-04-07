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
    parser.add_argument("--simple_output", help="whether only given one prediction value to each mutation", default=False, action="store_true")
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
        if not args.simple_output:
            save_file.writelines("chrom,mut_position,ref,alt,strand,position,reference_acceptor_ssu,reference_donor_ssu,pred_ref_acceptor_ssu,pred_ref_donor_ssu,pred_acceptor_deltassu,pred_donor_deltassu\n") 
        else:
            save_file.writelines("chrom,mut_position,strand,pred_acceptor_deltassu,pred_donor_deltassu\n") 
        for chrom, mut_pos,ref,alt, strand in zip(input_file["chrom"], input_file["mut_position"],input_file["ref"], input_file["alt"], input_file["strand"]):
            pos=mut_pos  
            seq_start=pos-(EL+CL)//2
            seq_end=seq_start+EL+CL
            if len(ref)>1 or len(alt)>1 or strand not in ["+", "-"]:
                if  args.simple_output:
                    save_file.writelines(f"{chrom},{mut_pos},{strand},Nan,Nan\n")
                
                continue
            seq=reference_genome[chrom][max(seq_start, 0):min(seq_end, len(reference_genome[chrom]))].upper()
            
            if seq_start<0:
                seq="N"*abs(seq_start)+seq
            if seq_end>len(reference_genome[chrom]):
                seq=seq+"N"*abs(seq_start)
            assert seq[mut_pos-seq_start]==ref.upper(), f"{ref} {seq[mut_pos-seq_start-1:mut_pos-seq_start+2]}"
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
            max_acceptor_impact=0
            max_donor_impact=0
            
            positions=[]
            acceptor_refs=[]
            donor_refs=[]
            acceptor_ref_preds=[]
            donor_ref_preds=[]
            acceptor_deltas=[]
            donor_deltas=[]
            for i in range(write_window_start, write_window_end+1):
                acceptor_ref=refmat[0, i, 1]
                donor_ref=refmat[0, i, 2]
                acceptor_ref_pred=pred_ref[0, i, 1]
                donor_ref_pred=pred_ref[0, i, 2]
                acceptor_delta=pred_delta[0, i, 1]
                donor_delta=pred_delta[0, i, 2]
                if abs(acceptor_delta)>abs(max_acceptor_impact):
                    max_acceptor_impact=acceptor_delta
                if abs(donor_delta)>abs(max_donor_impact):
                    max_donor_impact=donor_delta
                
                position=pos-args.window_size//2+i-write_window_start
                
                positions.append(str(position))
                acceptor_refs.append(str(acceptor_ref))
                donor_refs.append(str(donor_ref))
                acceptor_ref_preds.append(str(acceptor_ref_pred))
                donor_ref_preds.append(str(donor_ref_pred))
                acceptor_deltas.append(str(acceptor_delta))
                donor_deltas.append(str(donor_delta))
            
            if not args.simple_output:
                save_file.writelines(f'{chrom},{mut_pos},{ref},{alt},{strand},{";".join(positions)},{";".join(acceptor_refs)},{";".join(donor_refs)},{";".join(acceptor_ref_preds)},{";".join(donor_ref_preds)},{";".join(acceptor_deltas)},{";".join(donor_deltas)}\n')
            
            if args.simple_output:
                save_file.writelines(f"{chrom},{mut_pos},{ref},{alt},{strand},{max_acceptor_impact},{max_donor_impact}\n")
        save_file.close()
    else:
        save_file.writelines("chrom,mut_position,ref,alt,strand,reference_acceptor_ssu,reference_donor_ssu,pred_ref_acceptor_ssu,pred_ref_donor_ssu,pred_acceptor_deltassu,pred_donor_deltassu\n")
        for chrom, mut_pos,ref,alt, strand,jn_start, jn_end, acceptor_ssu, donor_ssu in zip(input_file["chrom"], input_file["mut_position"],input_file["ref"], input_file["alt"], input_file["strand"], input_file["exon_end"], input_file["exon_start"], input_file["acceptor_ssu"], input_file["donor_ssu"]):
            pos=(jn_start+jn_end)//2        
            seq_start=pos-(EL+CL)//2
            seq_end=seq_start+EL+CL
            acceptor_ssu=max(acceptor_ssu, 1e-3)
            donor_ssu=max(donor_ssu, 1e-3)
            
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
                    ssu_end=np.array(default_anno_info[species][chrom][strand][str(jn_end)])
                    ssu_end[np.isnan(ssu_end)]=1e-3
                else:
                    print(f"{jn_end} not in defalt reference ssu")
                    ssu_end=np.zeros(3)
                    
                if str(jn_start) in default_anno_info[species][chrom][strand]:
                    ssu_start= np.array(default_anno_info[species][chrom][strand][str(jn_start)])
                    ssu_start[np.isnan(ssu_start)]=1e-3
                else:
                    print(f"{jn_start} not in defalt reference ssu")
                    ssu_start=np.zeros(3)
            else:
                ssu_start=np.zeros(3)
                ssu_end=np.zeros(3)
                if strand=="-":
                    ssu_end[..., 2]=donor_ssu
                    ssu_start[..., 1]=acceptor_ssu
                else:
                    ssu_end[..., 1]=acceptor_ssu
                    ssu_start[..., 2]=donor_ssu
                
            if strand=="-":
                seq=[repdict[_] for _ in seq][::-1]
                mutseq=[repdict[_] for _ in mutseq][::-1]
                if ssu_end[..., 2]==0:
                    ssu_end[..., 2]=1e-3
                if ssu_start[..., 1]==0:
                    ssu_start[..., 1]=1e-3
                refmat[jn_end-seq_start]=ssu_end
                refmat[jn_start-seq_start]=ssu_start
                refmat=refmat[::-1]
            else:
                if ssu_end[..., 1]==0:
                    ssu_end[..., 1]=1e-3
                if ssu_start[..., 2]==0:
                    ssu_start[..., 2]=1e-3
                refmat[jn_start-seq_start]=ssu_start
                refmat[jn_end-seq_start]=ssu_end
                
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
            
            
            pred_acceptor_ref=(pred_ref[:, :, 1].reshape(-1)[refmat[:, 1]>0]).mean()
            pred_donor_ref=(pred_ref[:, :, 2].reshape(-1)[refmat[:, 2]>0]).mean()
            pred_acceptor_delta=(pred_delta[:, :, 1].reshape(-1)[refmat[:, 1]>0]).mean()
            pred_donor_delta=(pred_delta[:, :, 2].reshape(-1)[refmat[:, 2]>0]).mean()
            acceptor_ref=(refmat[..., 1].reshape(-1)[refmat[:, 1]>0]).mean()
            donor_ref=(refmat[..., 2].reshape(-1)[refmat[:, 2]>0]).mean()

            save_file.writelines(f"{chrom},{mut_pos},{ref},{alt},{strand},{acceptor_ref},{donor_ref},{pred_acceptor_ref},{pred_donor_ref},{pred_acceptor_delta},{pred_donor_delta}\n")
        save_file.close()
            
        
    

if __name__ == "__main__":
    main()
