# eSplice
Pytorch version of eSplice

## Generate train/test/valid data from bed file
1. Please refer to gene_dataset.tsu.txt for the format of bed file.
2. Change the custom path in config.py
3. Run
>>>
    python -m Tools.annotate_gene
    python -m Tools.generate_data
>>>

## Run model training/evaluation
1. Please refer to example/ for the format of config/test_config/mut_config file
2. Run
>>>
    python main.py -c path/to/config
    ## examples
    # train a model: python main.py -c tasks/Pretrain_withcons/config
    # test a model: python main.py -c tasks/Pretrain_withcons/test_config
>>>

## Evaluate on mutation data
1. Please refer to VexSeq_snps2exon_ref_dataset.txt for the format of input data.
2. Generate data by
>>>
    python -m Tools.generate_mutdata /the/path/to/raw/mutation/data /the/path/to/save/processed/data the_name_of_reference_genome
    ## example
    #python -m Tools.generate_mutdata /temp/xuchencheng/eSplicedata/Mutations/VexSeq_snps2exon_ref_dataset.txt /temp/xuchencheng/eSplicedata/Species_37_ProcessedData/VexSeq/ hg19
>>>
3. Set the path in mut_config and run evaluation as
>>>
    python main_mut.py -c /the/path/to/mut_config
>>>

## Quick start with pretrained model
1. Prepare data
>>>
    bash example/preparedata.sh
>>>
2. Run evaluation with conservation score
>>>
    python main.py -c example/test_config_wcons
    python main_mut.py -c example/mut_config_wcons
    # check data under TestResultWCons/
>>>
3. Run evaluation without conservation score
>>>
    python main.py -c example/test_config_wocons
    python main_mut.py -c example/mut_config_wocons
    # check data under TestResultWOCons/
>>>