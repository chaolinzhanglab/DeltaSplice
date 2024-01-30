source /home/xuc0d/anaconda3/bin/activate openfold
python -m experiments.7_eval_aso.generate_aso
CUDA_VISIBLE_DEVICES=1 python main_mut.py -c experiments/7_eval_aso/refsplice_config

CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/IKBKAP/data_withlabel.json experiments/7_eval_aso/test_results/IKBKAP_pangolin_withlabel.txt
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/SMN2/data_withlabel.json experiments/7_eval_aso/test_results/SMN2_pangolin_withlabel.txt
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/IKBKAP/data.json experiments/7_eval_aso/test_results/IKBKAP_pangolin.txt
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/SMN2/data.json experiments/7_eval_aso/test_results/SMN2_pangolin.txt

module load cudnn
source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/IKBKAP/data_withlabel.json  experiments/7_eval_aso/test_results/IKBKAP_spliceai_withlabel.txt
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/SMN2/data_withlabel.json  experiments/7_eval_aso/test_results/SMN2_spliceai_withlabel.txt
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/IKBKAP/data.json  experiments/7_eval_aso/test_results/IKBKAP_spliceai.txt
CUDA_VISIBLE_DEVICES=0 python baselines/spliceai.py data/SMN2/data.json  experiments/7_eval_aso/test_results/SMN2_spliceai.txt

source /home/xuc0d/anaconda3/bin/activate mmsplice
python baselines/run_mmsplice.py  data/IKBKAP/data_withlabel.json py36hgfile/hg19.fa experiments/7_eval_aso/test_results/IKBKAP_mmsplice_withlabel.txt
python baselines/run_mmsplice.py  data/SMN2/data_withlabel.json py36hgfile/hg19.fa experiments/7_eval_aso/test_results/SMN2_mmsplice_withlabel.txt
python baselines/run_mmsplice.py  data/IKBKAP/data.json py36hgfile/hg19.fa experiments/7_eval_aso/test_results/IKBKAP_mmsplice.txt
python baselines/run_mmsplice.py  data/SMN2/data.json py36hgfile/hg19.fa experiments/7_eval_aso/test_results/SMN2_mmsplice.txt

