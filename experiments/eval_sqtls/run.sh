rm -r experiments/eval_sqtls/test_results

python -m experiments.eval_sqtls.generate_sqtl
python -m experiments.eval_sqtls.generate_neg
python main.py --save_path experiments/eval_sqtls --mut_data_path data/sQTLs/data.json  data/sQTLs_neg/data.json  --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True
