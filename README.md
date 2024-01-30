# DeltaSplice
Code for the paper "Reference-informed prediction of alternative splicing and splicing-altering mutations from sequences". All experiments mentioned in the manuscript can be reproduced with the scripts under the experiments/ folder.

## Data preparation
1. models for MMSplice: 
>>>
    cd baselines
    git clone https://github.com/gagneurlab/MMSplice_paper.git
    cd ..
>>>
2. models for pangolin:
>>>
    cd baselines
    git clone https://github.com/tkzeng/Pangolin.git
    mv Pangolin/pangolin/models/ pangolin_models
    cd ..
>>>
3. models for spliceai
>>>
    cd baselines
    git clone https://github.com/Illumina/SpliceAI.git
    mv SpliceAI/spliceai/models spliceai_models
    cd ..
>>>
4. data used in this work is avaiable at xxxx

## Generate train/test/valid data from bed file
1. Please refer to gene_dataset.tsu.txt for the format of bed file.
2. Change the custom path in config.py
3. Run
>>>
    python -m Tools.annotate_gene
    python -m Tools.generate_data
>>>

## Run model training/evaluation
1. Please refer to configs under folder tasks for the format of config/test_config/mut_config file
2. Run
>>>
    python main.py -c path/to/config
    ## examples
    # train a model: python main.py -c tasks/Pretrain_withcons/config
    # test a model: python main.py -c tasks/Pretrain_withcons/test_config
>>>

## Reproduce experiments mentioned in the manuscript
1. All scripts for experiments are under the folder experiments. In each folder, run.sh contains all the command lines. Directly run "bash run.sh" can generate all the results.
2. experiments/9_plot_and_merge.ipynb is used to summarize all the generated results.

## Quick start with pretrained model
1. Please refer to VexSeq_snps2exon_ref_dataset.txt for the format of input data.
2. Write the config file following experiments/2_eval_mut/RefSplice_mut_config.py.
3. Generate file with:
>>>
    python -m Tools.generate_mutdata /path/to/data /path/to/save reference genome
    # example
    # python -m Tools.generate_mutdata data/VexSeq_snps2exon_ref_dataset.txt data/vexseq/ hg19 
>>>
4. Run the valuation with:
>>>
    python main_mut.py -c /the/path/to/mut_config
>>>
