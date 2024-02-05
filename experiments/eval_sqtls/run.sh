rm -r experiments/eval_sqtls/test_results

python -m experiments.eval_sqtls.generate_sqtl
python -m experiments.eval_sqtls.generate_neg
python main_mut.py -c experiments/eval_sqtls/refsplice_config
python main_mut_single.py -c experiments/eval_sqtls/refsplice_config
python baselines/pangolin.py data/sQTLs/data.json experiments/eval_sqtls/test_results/sqtl_pangolin.txt
python baselines/pangolin.py data/sQTLs_neg/data.json experiments/eval_sqtls/test_results/neg_pangolin.txt


python baselines/spliceai.py data/sQTLs/data.json experiments/eval_sqtls/test_results/sqtl_spliceai.txt
python baselines/spliceai.py data/sQTLs_neg/data.json experiments/eval_sqtls/test_results/neg_spliceai.txt

python baselines/run_mmsplice.py  data/sQTLs/data.json py36hgfile/hg38.fa experiments/eval_sqtls/test_results/sqtl_mmsplice.txt
python baselines/run_mmsplice.py  data/sQTLs_neg/data.json py36hgfile/hg38.fa experiments/eval_sqtls/test_results/neg_mmsplice.txt
