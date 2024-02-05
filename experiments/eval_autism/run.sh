python -m experiments.eval_autism.generate_autism_exome
python -m experiments.eval_autism.generate_autism_genome
python main_mut.py -c experiments/eval_autism/refsplice_config 
python main_mut_single.py -c experiments/eval_autism/refsplice_config 
python baselines/pangolin.py data/autism_exome/data.json experiments/eval_autism/test_results/exome_pangolin.txt
python baselines/pangolin.py data/autism_genome/data.json experiments/eval_autism/test_results/genome_pangolin.txt


python baselines/spliceai.py data/autism_exome/data.json  experiments/eval_autism/test_results/exome_spliceai.txt
python baselines/spliceai.py data/autism_genome/data.json  experiments/eval_autism/test_results/genome_spliceai.txt

python baselines/run_mmsplice.py  data/autism_exome/data.json py36hgfile/hg19.fa experiments/eval_autism/test_results/exome_mmsplice.txt
python baselines/run_mmsplice.py  data/autism_genome/data.json py36hgfile/hg19.fa experiments/eval_autism/test_results/genome_mmsplice.txt
