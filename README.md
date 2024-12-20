# DeltaSplice 

A neural network model to predict splice site usage and splicing-altering mutations

Citation:

Xu, C., Bao, S., Wang, Y., Li, W., Chen, H., Shen, Y., ... & Zhang, C. (2024). Reference-informed prediction of alternative splicing and splicing-altering mutations from sequences. *Genome Research*, 34(7), 1052-1065.

## Update
Now we align the input format with [SpliceAI](https://github.com/Illumina/SpliceAI/tree/master). Installation:

>>>
    cd deltasplice/data/anno
    tar -xzvf data.tar.gz
    cd ../../..
    bash Tools/download_files.sh
    pip install .
>>>

Usage:
>>>
    deltasplice -I test.vcf -O temp.vcf -R data/genomes/hg19.fa  -A grch37 -U 1 -M 0
>>>

Required parameters:

    -I: Input VCF with variants of interest.
    -O: Output VCF with deltasplice predictions ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL included in the INFO column (see table below for details). Only SNVs and simple INDELs (REF or ALT is a single base) within genes are annotated. Variants in multiple genes have separate predictions for each gene.
    -R: Reference genome fasta file. Can be downloaded from GRCh37/hg19 or GRCh38/hg38.
    -A: Gene annotation file. Can instead provide grch37 or grch38 to use GENCODE V24 canonical annotation files included with the package. To create custom gene annotation files, use deltasplice/annotations/grch37.txt in repository as template.
Optional parameters:

    -D: Maximum distance between the variant and gained/lost splice site (default: 50).
    -U: Whether to use default reference sequences and reference usages (default: 1).
    -M: Mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss (default: 0).

## Installation
[Anaconda](https://www.anaconda.com/download)  is recommended for installation. 
>>>
    git clone https://github.com/chaolinzhanglab/DeltaSplice.git
    cd DeltaSplice

    # init conda environment 
    conda create -n deltasplice python=3.9 bioconda::pyfasta
    conda activate deltasplice

    # install torch
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

    # install dependents
    pip install -r requirements.txt
>>>


## Data preparation

Download genome reference and liftOver files from UCSC and save them to `data/genomes`.


>>>
    bash Tools/download_files.sh
>>>

## Quick start with pretrained model
Currently DeltaSplice support the prediction of SSU and delta-SSU for mutations. Example data are provided under `data/` and pretrained models are under `pretrained_models`. The file `deltasplice/constant.py` contains the default path to pretrained models and reference genomes. The average prediction of five models under `pretrained_models/DeltaSplice_models/` is used as the final prediction for SSU and delta-SSU. Here, splice acceptor refers to the second nucleotide of the corresponding exon, and splice donor refers to the second-to-last nucluotide of the corresponding exon.

The description of the header in the output file of DeltaSplice is as follows:

    |           ID           |                        Description                        |
    |------------------------|-----------------------------------------------------------|
    |         chrom          | id of chromosome                                          |
    |        position        | zero-based coordinates                                    |
    |         strand         | strand                                                    |
    |      acceptor_ssu      | predicted SSU (acceptor)                                  |
    |       donor_ssu        | predicted SSU (donor)                                     |
    |      mut_position      | zero-based coordinates of mutation sites                  |
    |           ref          | reference base before mutation                            |
    |           alt          | alternative base after mutation                           |
    | reference_acceptor_ssu | reference SSUs used in mutation prediction (acceptor)     |
    |  reference_donor_ssu   | reference SSUs used in mutation prediction (donor)        |
    |  pred_ref_acceptor_ssu | predicted SSUs for the reference allele (acceptor)        |
    |   pred_ref_donor_ssu   | predicted SSUs for the reference allele (donor)           |
    | pred_acceptor_deltassu | predicted delta-SSUs for mutations (acceptor)             |
    |   pred_donor_deltassu  | predicted delta-SSUs for mutations (donor)                |

It is worth noting that when predicting SSU, we given both the predicted acceptor SSU and donor SSU. However, the predicted score holds significance only when it exceeds 1e-3. When predicting delta-SSU values for mutations, if specific exons are provided, DeltaSplice exclusively predicts the delta-SSU for these exons. In cases where exons are not specified, DeltaSplice predicts both acceptor and donor delta-SSU values.

### SSU prediction

For the prediction of SSU, the input file should be in the csv format with chrom, zero-based position and strand, as follows,

    | chrom   | position | strand |
    |---------|----------|--------|
    | chr1    | 151003847| +      |
    | chr1    | 176015316| -      |

The output file contains chrom, position, strand, acceptor_ssu, and donor_ssu columns.
#### Usage:

Run following code to generate prediction results, in which the reference genome is used to extract input sequences.

>>>
    python pred_ssu.py --data_path /path/to/input.csv --save_path /path/to/output.csv --genome reference_genome
>>>

Required parameters:

 - ```--data_path```: Input CSV file with coordinates and strands of sites to be predicted. Please refer to `data/example_pred_ssu.csv`.
 - ```--save_path```: Output CSV file with prediction results. The output file contains five columns, i.e. chrom, position, strand,acceptor_ssu and donor_ssu, where acceptor_ssu and donor_ssu are predicted SSUs for each site when it is used as acceptor or donor, respectively. Sites with SSU predicted values lower than 1e-3 are not considered as splicing sites.
 - ```--genome```   : Which reference genome to use, for example, hg19, hg38 or other reference genomes. Note that the default path for reference genome is `data/genomes`.

#### Example:

>>>
    python pred_ssu.py --data_path data/example_pred_ssu.csv --save_path data/example_ssu_pred_out.csv  --genome hg19 
>>>

### Delta-SSU prediction without exon information
For the prediciton of delta-SSU for mutations without given exon information, the input file should be in csv format and contain chrom, mut_position, ref, alt and strand columns, as shown in the following table. The mut_position should be zero-based, and ref/alt are bases on the positive strand.

    | chrom   | mut_position | ref | alt | strand |
    |---------|--------------|-----|-----|--------|
    | chr1    | 114161115    | G   | A   | +      |
    | chr1    | 119584866    | G   | C   | -      |

The output file contains chrom, mut_position, ref, alt, strand, position, reference_acceptor_ssu, reference_donor_ssu, pred_ref_acceptor_ssu, pred_ref_donor_ssu, pred_acceptor_deltassu and pred_donor_deltassu columns.

#### Usage:

  Run following code to generate prediction results
>>>
    # Exon positions are unspecified and do not use reference information
    python pred_deltassu.py --data_path /path/to/data --save_path /path/to/save  --window_size 200  --genome reference_genome 
    #  Exon positions are unspecified and use SSU values of brain tissues from RNA-seq data as reference
    python pred_deltassu.py --data_path /path/to/data --save_path /path/to/save  --window_size 200  --genome reference_genome --use_reference
>>>

Required parameters:

 - ```--data_path```: Input CSV file with coordinates, ref/alt bases, strands and exon positions. Please refer to `data/vexseq_out.csv`.
 - ```--save_path```: Output CSV file with prediction results. The output file contains eight columns, i.e. chrom, mut_position, strand, exon_start, exon_end, ssu, pred_ref, pred_deltassu, where pred_ref is the predicted SSU for the sequence before mutation, and pred_deltassu is the predicted delta-SSU for current mutation.
 - ```--window_size```: Predicted window size around mutation sites, the default value is 200.
 - ```--use_reference```: Whether use the default SSU information of brain tissues from RNA-seq data. The default value is False.
 - ```--genome```   : Which reference genome to use, for example, hg19, hg38 or other reference genomes. The default path for reference genome is `data/genomes`.
#### Example:

>>>
    # Exon positions are unspecified and do not use reference information
    python pred_deltassu.py --data_path data/vexseq_wo_exon.csv  --save_path data/vexseq_woexon_out.csv --window_size 200  --genome hg19 
    # Exon positions are unspecified and use SSU values of brain tissues from RNA-seq data as reference
    python pred_deltassu.py --data_path data/vexseq_wo_exon.csv  --save_path data/vexseq_woexon_wref_out.csv --window_size 200  --genome hg19 --use_reference
>>>

### Delta-SSU prediction with exon information

For the prediction of delta-SSU for mutations with given exons, the input file should be in csv format and contain the following columns, in which if there's no SSU information, set acceptor_ssu and donor_ssu as Nan. Note that all positions should be zero-based. Here ref/alt are bases on the positive strand.

    | chrom   | mut_position | ref | alt | strand | exon_start | exon_end | acceptor_ssu|donor_ssu|
    |---------|--------------|-----|-----|--------|------------|----------|-------------|---------|
    | chr1    | 114161115    | G   | A   | +      | 114161153  |114161227 |   0.4888    | 0.4888  |
    | chr1    | 119584866    | G   | C   | -      | 119584886  |119584971 |   0.8859    | 0.8859  |

The output file contains chrom, mut_position, ref, alt, strand, reference_acceptor_ssu, reference_donor_ssu, pred_ref_acceptor_ssu, pred_ref_donor_ssu, pred_acceptor_deltassu and pred_donor_deltassu columns.
#### Usage:
Note that 5' splice site is represented by the upstream nucleotide (junction start) and 3' splice site is represented by the downstream nucleotide (junction end).
  Run following code to generate prediction results
>>>
    # do not use SSU values in brain tissues from RNA-seq data as reference information
    python pred_deltassu.py --data_path /path/to/data --save_path /path/to/save  --genome reference_genome 
    # use SSU values in brain tissues from RNA-seq data as reference information
    python pred_deltassu.py --data_path /path/to/data --save_path /path/to/save  --genome reference_genome --use_reference
>>>

Required parameters:

 - ```--data_path```: Input CSV file with coordinates, ref/alt bases, strands and exon positions. Please refer to `data/vexseq_out.csv`.
 - ```--save_path```: Output CSV file with prediction results. The output file contains eight columns, i.e. chrom, mut_position, strand, exon_start, exon_end, ssu, pred_ref, pred_deltassu, where pred_ref is the predicted SSU for the sequence before mutation, and pred_deltassu is the predicted delta-SSU for current mutation.
 - ```--use_reference```: Whether use the default usage information, the default value is False.
 - ```--genome```   : Which reference genome to use, for example, hg19, hg38 or other reference genomes. The default path for reference genome is `data/genomes`.
#### Example:

>>>
    # do not use SSU values in brain tissues from RNA-seq data as reference information
    python pred_deltassu.py --data_path data/vexseq.csv  --save_path data/vexseq_out.csv --genome hg19 
    # use SSU values in brain tissues from RNA-seq data as reference information
    python pred_deltassu.py --data_path data/vexseq.csv  --save_path data/vexseq_out_wref.csv --genome hg19 --use_reference
>>>

## Train models from scratch
We provided data and scripts for users to train the model from scratch, and evaluate the performance of the obtained model. 

### Prepare train/test/valid data from gene annotation file
- `data/gene_dataset.tsu.txt` contains splice site usage in the adult brains of eight mammalian species.
- Run the following code to generate necessary data for model training and evaluation.
>>>
    #Generate gene annotations on the genome
    python -m Tools.annotate_gene --save_path data/anno/ --input_file data/gene_dataset.tsu.txt

    #Generate data for training/testing/validation
    python -m Tools.generate_data --save_path data/train_val_test --anno_path data/anno --input_file data/gene_dataset.tsu.txt
>>>

### Run model training/evaluation
The script to reproduce model training is experiments/model_train/run.sh, which will train 5 models with different random seeds with the same data. You can run it directly by `bash experiments/model_train/run.sh`. 

The python command lines to train models are as follows,
>>>
    python main.py --save_path experiments/model_train/DeltaSplice_rep0/ --is_train=True --train_data_path=data/train_val_test/train/data.json --valid_data_path=data/train_val_test/valid/data.json --seed=321

    python main.py --save_path experiments/model_train/DeltaSplice_rep1/ --is_train=True --train_data_path=data/train_val_test/train/data.json --valid_data_path=data/train_val_test/valid/data.json --seed=49

    python main.py --save_path experiments/model_train/DeltaSplice_rep2/ --is_train=True --train_data_path=data/train_val_test/train/data.json --valid_data_path=data/train_val_test/valid/data.json --seed=976

    python main.py --save_path experiments/model_train/DeltaSplice_rep3/ --is_train=True --train_data_path=data/train_val_test/train/data.json --valid_data_path=data/train_val_test/valid/data.json --seed=753

    python main.py --save_path experiments/model_train/DeltaSplice_rep4/ --is_train=True --train_data_path=data/train_val_test/train/data.json --valid_data_path=data/train_val_test/valid/data.json --seed=491
>>>
where 

 - ```--save_path```: The path to save the generated files in the training process.
 - ```--is_train```: Set the state as training mode.
 - ```--train_data_path```: The path to training data, which is generated by `Tools.generate_data`.
 - ```--valid_data_path```: The path to validation data, which is generated by `Tools.generate_data`.
- ```--seed```: The random seed.

An example command line to evaluate the performance of the trained models to predict SSU is as follows,
>>>
    # example
    python main.py --save_path experiments/evaluate_on_test_and_val --test_data_path data/train_val_test/test/data.json  data/train_val_test/test/human.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4    
>>>
Similarly, where 

 - ```--save_path```: The path to save the generated files in the testing process.
 - ```--test_data_path```: One or more paths to the testing data, which should be in the json format, as generated by `Tools.generate_data`.
 - ```--load_model_path```: One or more paths to load saved model.

## Analysis of additional datasets

 Scripts used to reproduce experiments described in our paper are under `experiments/`, where the folder `experiments/model_train` contains the script to train models, and other folders contain scripts for model evaluations (the default pretrained models under `pretrained_models/DeltaSplice_models` are used). In detail, 
 - `experiments/eval_on_test_and_val` contains script to run evaluations in testing data and validation data.
 - `experiments/eval_mut` contains script to run evaluations on MFASS dataset and VexSeq dataset.
 - `experiments/eval_sqtls` constains script to run evaluations on sQTLs data.
 - `experiments/eval_fas` contains script to run evaluations on FAS gene.
 - `experiments/eval_autism` contains script to run evaluations on Autism data.