rm -r experiments/1_evaluate_on_test_and_val/test_results
source /home/xuc0d/anaconda3/bin/activate openfold
python -m Tools.annotate_gene
python -m Tools.generate_data
CUDA_VISIBLE_DEVICES=0 python main.py -c experiments.1_evaluate_on_test_and_val.RefSplice_test_config
CUDA_VISIBLE_DEVICES=0 python main.py -c experiments.1_evaluate_on_test_and_val.pangolin_test_config
CUDA_VISIBLE_DEVICES=0 python main.py -c experiments.1_evaluate_on_test_and_val.RefSplice_human_test_config
module load cudnn
source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=0 python baselines/spliceai_eval.py data/Npy/test/data.json  experiments/1_evaluate_on_test_and_val/test_results/test
CUDA_VISIBLE_DEVICES=0 python baselines/spliceai_eval.py data/Npy/valid/data.json  experiments/1_evaluate_on_test_and_val/test_results/valid

