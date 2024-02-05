python -m experiments.eval_fas.generate_data
python main_mut.py -c experiments/eval_fas/RefSplice_mut_config
python baselines/pangolin.py data/FAS/data.json experiments/eval_fas/test_results/fas_pangolin.txt

python baselines/spliceai.py data/FAS/data.json experiments/eval_fas/test_results/fas_spliceai.txt

python baselines/run_mmsplice.py  data/FAS/data.json py36hgfile/hg19.fa experiments/eval_fas/test_results/fas_mmsplice.txt

