
python -m experiments.eval_on_multiple_species.generate_data
python -m experiments.eval_on_multiple_species.main
python -m experiments.eval_on_multiple_species.main_aba

python baselines/spliceai_s.py data/Hg19VsOthers/ experiments/eval_on_multiple_species/test_results