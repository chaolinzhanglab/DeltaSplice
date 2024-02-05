rm -r experiments/eval_sqtls/test_results
source /home/xuc0d/anaconda3/bin/activate openfold
python -m experiments.eval_sqtls.generate_sqtl
python -m experiments.eval_sqtls.generate_neg
CUDA_VISIBLE_DEVICES=0 python main_mut.py -c experiments/eval_sqtls/refsplice_config
CUDA_VISIBLE_DEVICES=1 python main_mut_single.py -c experiments/eval_sqtls/refsplice_config
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/sQTLs/data.json experiments/eval_sqtls/test_results/sqtl_pangolin.txt
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/sQTLs_neg/data.json experiments/eval_sqtls/test_results/neg_pangolin.txt

module load cudnn
source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=0 python baselines/spliceai.py data/sQTLs/data.json experiments/eval_sqtls/test_results/sqtl_spliceai.txt
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/sQTLs_neg/data.json experiments/eval_sqtls/test_results/neg_spliceai.txt

source /home/xuc0d/anaconda3/bin/activate mmsplice
python baselines/run_mmsplice.py  data/sQTLs/data.json py36hgfile/hg38.fa experiments/eval_sqtls/test_results/sqtl_mmsplice.txt
python baselines/run_mmsplice.py  data/sQTLs_neg/data.json py36hgfile/hg38.fa experiments/eval_sqtls/test_results/neg_mmsplice.txt
