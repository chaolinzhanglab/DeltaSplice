# DeltaSplice 

A neural network model to predict splice site usage and splicing-altering mutations

Citation:

Xu, C., Bao, S., Chen, H., Jiang, T., Zhang, C. "Reference-informed prediction of alternative splicing and splicing-altering mutations from sequences." *In submission*.

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
Currently DeltaSplice support the prediction of SSU and delta-SSU for mutations. Example data are provided under `data/` and pretrained models are under `pretrained_models`. The file `deltasplice/constant.py` contains the default path to pretrained models and reference genomes. The average prediction of five models under `pretrained_models/DeltaSplice_models/` is used as the final prediction for SSU and delta-SSU.

### SSU prediction

For the prediction of SSU, the input file should be in the csv format with chrom, zero-based position and strand, as follows,

    | chrom   | position | strand |
    |---------|----------|--------|
    | chr1    | 151003847| +      |
    | chr1    | 176015316| -      |


#### Usage:

Run following code to generate prediction results, in which the reference genome is used to extract input sequences.

>>>
    python pred_ssu.py --data_path /path/to/input.csv --save_path /path/to/output.csv --genome reference_genome
>>>

Required parameters:

 - ```--data_path```: Input CSV file with coordinates and strands of sites to be predicted, as mentioned before.
 - ```--save_path```: Output CSV file with prediction results. The output file contains five columns, i.e. chrom, position, strand,acceptor_ssu and donor_ssu, where acceptor_ssu and donor_ssu are predicted SSUs for each site when it is used as acceptor or donor, respectively. Sites with SSU predicted values lower than 1e-3 are not considered as splicing sites.
 - ```--genome```   : Which reference genome to use, for example, hg19, hg38 or other reference genomes. Note that the default path for reference genome is `data/genomes`.

#### Example:

>>>
    python pred_ssu.py --data_path data/example_pred_ssu.csv --save_path data/example_ssu_pred_out.csv --genome hg19 
>>>

### Delta-SSU prediction

For the prediction of delta-SSU for mutations, the input file should be in csv format and contain the following columns, in which if there's no psi information, set psi as Nan. Note that all positions should be zero-based. Here psi means psi of the reference allele, and ref/alt are bases on the positive strand.

    | chrom   | mut_position | ref | alt | strand | exon_start | exon_end | psi   |
    |---------|--------------|-----|-----|--------|------------|----------|-------|
    | chr1    | 114161115    | G   | A   | +      | 114161153  |114161227 | 0.4888|
    | chr1    | 119584866    | G   | C   | -      | 119584886  |119584971 | 0.8859|

#### Usage:

  Run following code to generate prediction results
>>>
    python pred_deltassu.py --data_path /path/to/data --save_path /path/to/save --genome reference_genome
>>>

Required parameters:

 - ```--data_path```: Input CSV file with coordinates, ref/alt bases, strands and exon positions, as mentioned before.
 - ```--save_path```: Output CSV file with prediction results. The output file contains eight columns, i.e. chrom, mut_position, strand, exon_start, exon_end, psi, pred_ref, pred_deltassu, where pred_ref is the predicted SSU for the sequence before mutation, and pred_deltassu is the predicted delta-SSU for current mutation.
 - ```--genome```   : Which reference genome to use, for example, hg19, hg38 or other reference genomes. The default path for reference genome is `data/genomes`.
#### Example:

>>>
    python pred_deltassu.py --data_path data/vexseq.csv  --save_path data/vexseq_out.csv --genome hg19 
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

 Scripts used to reproduce experiments described in our paper are under `experiments/`.
