## DeltaSplice 

A neural network model to predict splice site usage and splicing-altering mutations

Citation:

Xu, C., Bao, S., Chen, H., Jiang, T., Zhang, C. "Reference-informed prediction of alternative splicing and splicing-altering mutations from sequences." *In submission*.

### Installation

>>>
    git clone https://github.com/chaolinzhanglab/DeltaSplice.git
    cd DeltaSplice
>>>


### Data preparation

Download genome reference and liftOver files from UCSC.
>>>
    bash Tools/download_files.sh
>>>

### Quick start with pretrained model
Currently DeltaSplice support the prediction of ssu for splice sites and delta-ssu for mutations
- For the prediction of ssu for splice sites, the input file should be in the csv format with chrom, zero-based position and strand, as follows,

    | chrom   | position | strand |
    |---------|----------|--------|
    | chr1    | 100      | +      |
    | chr2    | 250      | -      |
    | chr3    | 500      | +      |
    | chr4    | 750      | -      |


  Run following code to generate prediction results
>>>
    python pred_ssu.py --data_path /path/to/data --save_path /path/to/save --genome reference_genome
    # example
    # python pred_ssu.py --data_path data/example_pred_ssu.csv --save_path temp.csv --genome hg19 
>>>

- For the prediction of delta-ssu for mutations, the input file should be in csv format and contain the following columns, in which if there's no psi information, set psi as Nan. Note that all positions should be zero-based.
    | chrom   | mut_position | ref | alt | strand | jn_start | jn_end | psi  |
    |---------|--------------|-----|-----|--------|----------|--------|------|
    | chr1    | 100          | A   | T   | +      | 50       | 150    | 0.8  |
    | chr2    | 250          | C   | G   | -      | 200      | 300    | 0.5  |
    | chr3    | 500          | G   | A   | +      | 450      | 550    | 0.6  |
    | chr4    | 750          | T   | C   | -      | 700      | 800    | 0.9  |


  Run following code to generate prediction results
>>>
    python pred_deltassu.py --data_path /path/to/data --save_path /path/to/save --genome reference_genome
    # example
    # python pred_deltassu.py --data_path data/example_pred_deltassu.csv  --save_path temp.csv --genome hg19 
>>>


### Retrain the model using gene annotations

#### Prepare train/test/valid data from gene annotation file

- `gene_dataset.tsu.txt` contains splice site usage in the adult brains of eight mammalian species.
- Run
>>>
    #Generate gene annotations on the genome
    python -m Tools.annotate_gene --save_path data/anno/ --input_file data/gene_dataset.tsu.txt

    #Generate data for training/testing/validation
    python -m Tools.generate_data --save_path data/train_val_test --anno_path data/anno --input_file data/gene_dataset.tsu.txt
>>>

#### Run model training/evaluation
The script for model training is experiments/model_train/run.sh. In detail, to train a model:
>>>
    # train a model: 
    # example
    python main.py --save_path experiments/model_train/DeltaSplice_rep0/ --is_train=True --train_data_path=data/train_val_test/train/data.json --valid_data_path=data/train_val_test/valid/data.json --seed=321
>>>

To evaluate the performance of a model to predict ssu:
>>>
    # test a model: 
    # example
    python main.py --save_path experiments/evaluate_on_test_and_val --test_data_path data/train_val_test/test/data.json  data/train_val_test/test/human.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4    
>>>

To evaluate the performance of a model to predict delta-ssu

>>>
    # test a model: 
    # example
    python main.py --save_path experiments/eval_mut --mut_data_path data/vexseq/data.json  data/mfass/data.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True
>>>

### Reproduce experiments described in the manuscript

All experiments mentioned in the manuscript can be reproduced with the scripts under the `experiments/` folder.

In each folder, `run.sh` contains all the command lines. Directly run `bash run.sh` can generate all the results.




