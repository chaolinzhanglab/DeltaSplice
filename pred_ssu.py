from loguru import logger
import os
import argparse
import importlib
import random
import numpy as np
from deltasplice.utils import MutGenerator, GetSummaryStatisticsCallback, write_splice_site_file_header, write_splice_sites, get_top1_statistics, get_correlation, collect_predictions, density_scatter
import torch
from torch.utils.data import DataLoader
import copy
import pandas as pd
from deltasplice.constant import default_model_paths, model, Fapath, EL, SeqTable, repdict, IN_MAP
from pyfasta import Fasta

def main():
    # load config file
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_path", help="the path to input file")
    parser.add_argument("--save_path", help="the path to output file")
    parser.add_argument("--genome", help="reference genome")
    args = parser.parse_args()
    Models = [copy.deepcopy(model) for _ in default_model_paths]
    [m.load_state_dict(torch.load(b)) for m,b in zip(Models, default_model_paths)]
    
    input_file=pd.read_csv(args.data_path)
    reference_genome=Fasta(os.path.join(Fapath, args.genome+".fa"))
    save_file=open(args.save_path, "w")
    save_file.writelines("chrom,position,strand,accptor_ssu,donor_ssu\n")
    for chrom, pos, strand in zip(input_file["chrom"], input_file["position"], input_file["strand"]):
        seq_start=pos-EL//2
        seq_end=seq_start+EL+1
        seq=reference_genome[chrom][max(seq_start, 0):min(seq_end, len(reference_genome[chrom]))]
        if seq_start<0:
            seq="N"*abs(seq_start)+seq
        if seq_end>len(reference_genome[chrom]):
            seq=seq+"N"*abs(seq_start)
        seq=seq.upper()
        if strand=="-":
            seq=[repdict[_] for _ in seq][::-1]
        seq=IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
        pred=0
        for m in Models:
            pred+=m.predict({"X":torch.tensor(seq)[None]}, use_ref=False)["single_pred_psi"]
        pred=(pred/len(Models))[0]
        
        pred=pred[pred.shape[0]//2]
        save_file.writelines(f"{chrom},{pos},{strand},{pred[1]},{pred[2]}\n")
    save_file.close()
            
        
    

if __name__ == "__main__":
    main()
