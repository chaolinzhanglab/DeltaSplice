rm -r experiments/3_eval_fas/test_results
source /home/xuc0d/anaconda3/bin/activate openfold
python -m experiments.3_eval_fas.generate_data
CUDA_VISIBLE_DEVICES=1 python main_mut.py -c experiments/3_eval_fas/RefSplice_mut_config
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/FAS/data.json experiments/3_eval_fas/test_results/fas_pangolin.txt

source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/FAS/data.json experiments/3_eval_fas/test_results/fas_spliceai.txt

source /home/xuc0d/anaconda3/bin/activate mmsplice
python baselines/run_mmsplice.py  data/FAS/data.json py36hgfile/hg19.fa experiments/3_eval_fas/test_results/fas_mmsplice.txt

