source /home/xuc0d/anaconda3/bin/activate openfold
python -m experiments.0_eval_on_multiple_species.generate_data
CUDA_VISIBLE_DEVICES=0 python -m experiments.0_eval_on_multiple_species.main
#CUDA_VISIBLE_DEVICES=0 python -m experiments.0_eval_on_multiple_species.main_human  # this experiment is not necessary
CUDA_VISIBLE_DEVICES=0 python -m experiments.0_eval_on_multiple_species.main_aba
CUDA_VISIBLE_DEVICES=0 python -m experiments.0_eval_on_multiple_species.pangolin_human
module load cudnn
source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=0 python baselines/spliceai_s.py data/Hg19VsOthers/ experiments/0_eval_on_multiple_species/test_results