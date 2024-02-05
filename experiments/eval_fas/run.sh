python -m experiments.eval_fas.generate_data
python main_mut.py python main.py --save_path experiments/eval_fas --mut_data_path data/FAS/data.json  --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True
python baselines/pangolin.py data/FAS/data.json experiments/eval_fas/test_results/fas_pangolin.txt

python baselines/spliceai.py data/FAS/data.json experiments/eval_fas/test_results/fas_spliceai.txt

python baselines/run_mmsplice.py  data/FAS/data.json py36hgfile/hg19.fa experiments/eval_fas/test_results/fas_mmsplice.txt

