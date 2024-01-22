rm -r experiments/5_eval_autism/test_results
source /home/xuc0d/anaconda3/bin/activate openfold
python -m experiments.5_eval_autism.generate_autism_exome
python -m experiments.5_eval_autism.generate_autism_genome
CUDA_VISIBLE_DEVICES=1  python main_mut.py -c experiments/5_eval_autism/refsplice_config 
CUDA_VISIBLE_DEVICES=0 python main_mut_single.py -c experiments/5_eval_autism/refsplice_config 
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/autism_exome/data.json experiments/5_eval_autism/test_results/exome_pangolin.txt
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/autism_genome/data.json experiments/5_eval_autism/test_results/genome_pangolin.txt


source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/autism_exome/data.json  experiments/5_eval_autism/test_results/exome_spliceai.txt
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/autism_genome/data.json  experiments/5_eval_autism/test_results/genome_spliceai.txt

source /home/xuc0d/anaconda3/bin/activate mmsplice
python baselines/run_mmsplice.py  data/autism_exome/data.json py36hgfile/hg19.fa experiments/5_eval_autism/test_results/exome_mmsplice.txt
python baselines/run_mmsplice.py  data/autism_genome/data.json py36hgfile/hg19.fa experiments/5_eval_autism/test_results/genome_mmsplice.txt
