rm -r experiments/evaluate_on_test_and_val/test_results
python -m Tools.annotate_gene
python -m Tools.generate_data
python main.py -c experiments.evaluate_on_test_and_val.RefSplice_test_config
python main.py -c experiments.evaluate_on_test_and_val.pangolin_test_config
python main.py -c experiments.evaluate_on_test_and_val.RefSplice_human_test_config

python baselines/spliceai_eval.py data/Npy/test/data.json  experiments/evaluate_on_test_and_val/test_results/test
python baselines/spliceai_eval.py data/Npy/valid/data.json  experiments/evaluate_on_test_and_val/test_results/valid

