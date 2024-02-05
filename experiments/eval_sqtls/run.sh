rm -r experiments/eval_sqtls/test_results

python -m experiments.eval_sqtls.generate_sqtl
python -m experiments.eval_sqtls.generate_neg
python main_mut.py python main.py --save_path experiments/eval_sqtls --mut_data_path data/sQTLs/data.json  data/sQTLs_neg/data.json  --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True

python baselines/pangolin.py data/sQTLs/data.json experiments/eval_sqtls/test_results/sqtl_pangolin.txt
python baselines/pangolin.py data/sQTLs_neg/data.json experiments/eval_sqtls/test_results/neg_pangolin.txt


python baselines/spliceai.py data/sQTLs/data.json experiments/eval_sqtls/test_results/sqtl_spliceai.txt
python baselines/spliceai.py data/sQTLs_neg/data.json experiments/eval_sqtls/test_results/neg_spliceai.txt

python baselines/run_mmsplice.py  data/sQTLs/data.json py36hgfile/hg38.fa experiments/eval_sqtls/test_results/sqtl_mmsplice.txt
python baselines/run_mmsplice.py  data/sQTLs_neg/data.json py36hgfile/hg38.fa experiments/eval_sqtls/test_results/neg_mmsplice.txt
