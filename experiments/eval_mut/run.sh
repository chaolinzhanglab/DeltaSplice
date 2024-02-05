python -m Tools.generate_mutdata experiments/eval_mut/VexSeq_snps2exon_ref_dataset.txt data/vexseq/ hg19  #File, SavePath, species =
python -m Tools.generate_mutdata experiments/eval_mut/MFASS.txt data/mfass hg19
python main.py --save_path experiments/eval_mut --mut_data_path data/vexseq/data.json  data/mfass/data.json --load_model_path pretrained_models/DeltaSplice_models/model.ckpt-0 pretrained_models/DeltaSplice_models/model.ckpt-1 pretrained_models/DeltaSplice_models/model.ckpt-2 pretrained_models/DeltaSplice_models/model.ckpt-3 pretrained_models/DeltaSplice_models/model.ckpt-4   --use_reference=True
python baselines/pangolin.py data/vexseq/data.json experiments/eval_mut/test_results/vexseq_pangolin.txt
python baselines/pangolin.py data/mfass/data.json experiments/eval_mut/test_results/mfass_pangolin.txt

python baselines/spliceai.py data/vexseq/data.json experiments/eval_mut/test_results/vexseq_spliceai.txt
python baselines/spliceai.py data/mfass/data.json experiments/eval_mut/test_results/mfass_spliceai.txt

python baselines/run_mmsplice.py  data/vexseq/data.json py36hgfile/hg19.fa experiments/eval_mut/test_results/vexseq_mmsplice.txt
python baselines/run_mmsplice.py  data/mfass/data.json py36hgfile/hg19.fa experiments/eval_mut/test_results/mfass_mmsplice.txt


