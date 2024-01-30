python -m experiments.6_eval_clinvar.generate_clinvar
CUDA_VISIBLE_DEVICES=0 python main_mut.py -c experiments/6_eval_clinvar/refsplice_config
CUDA_VISIBLE_DEVICES=0 python baselines/pangolin.py data/clinvar/data.json experiments/6_eval_clinvar/test_results/pangolin.txt

source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=0 python baselines/spliceai.py data/clinvar/data.json  experiments/6_eval_clinvar/test_results/spliceai.txt

source /home/xuc0d/anaconda3/bin/activate mmsplice
python baselines/run_mmsplice.py  data/clinvar/data.json py36hgfile/hg19.fa experiments/6_eval_clinvar/test_results/mmsplice.txt

