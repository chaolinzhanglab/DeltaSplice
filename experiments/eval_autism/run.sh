python -m experiments.eval_autism.generate_autism_exome
python -m experiments.eval_autism.generate_autism_genome
python main.py --save_path experiments/eval_sqtls --mut_data_path data/autism_exome/data.json  data/autism_genome/data.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True
