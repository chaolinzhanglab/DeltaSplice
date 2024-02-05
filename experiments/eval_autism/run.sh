python -m experiments.eval_autism.generate_autism_exome
python -m experiments.eval_autism.generate_autism_genome
python main_mut.py python main.py --save_path experiments/eval_sqtls --mut_data_path data/autism_exome/data.json  data/autism_genome/data.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True

python baselines/pangolin.py data/autism_exome/data.json experiments/eval_autism/test_results/exome_pangolin.txt
python baselines/pangolin.py data/autism_genome/data.json experiments/eval_autism/test_results/genome_pangolin.txt


python baselines/spliceai.py data/autism_exome/data.json  experiments/eval_autism/test_results/exome_spliceai.txt
python baselines/spliceai.py data/autism_genome/data.json  experiments/eval_autism/test_results/genome_spliceai.txt

python baselines/run_mmsplice.py  data/autism_exome/data.json py36hgfile/hg19.fa experiments/eval_autism/test_results/exome_mmsplice.txt
python baselines/run_mmsplice.py  data/autism_genome/data.json py36hgfile/hg19.fa experiments/eval_autism/test_results/genome_mmsplice.txt
