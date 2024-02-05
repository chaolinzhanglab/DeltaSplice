python main.py --save_path experiments/evaluate_on_test_and_val --test_data_path data/train_val_test/test/data.json  data/train_val_test/valid/data.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4    
python main.py --save_path experiments/evaluate_on_test_and_val --test_data_path data/train_val_test/test/data.json  data/train_val_test/valid/data.json --load_model_path pretrained_models/DeltaSplice_human/model.ckpt-0 pretrained_models/DeltaSplice_human/model.ckpt-1 pretrained_models/DeltaSplice_human/model.ckpt-2 pretrained_models/DeltaSplice_human/model.ckpt-3 pretrained_models/DeltaSplice_human/model.ckpt-4    


python experiments/baselines/spliceai_eval.py data/train_val_test/test/data.json  experiments/evaluate_on_test_and_val/test_results/test
python experiments/baselines/spliceai_eval.py data/train_val_test/valid/data.json  experiments/evaluate_on_test_and_val/test_results/valid

