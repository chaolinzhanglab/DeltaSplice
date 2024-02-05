python -m Tools.generate_mutdata data/VexSeq_snps2exon_ref_dataset.txt data/vexseq/ hg19  #File, SavePath, species =
python -m Tools.generate_mutdata data/MFASS.bed data/mfass hg19
python main_mut.py -c experiments/eval_mut/RefSplice_mut_config
python baselines/pangolin.py data/vexseq/data.json experiments/eval_mut/test_results/vexseq_pangolin.txt
python baselines/pangolin.py data/mfass/data.json experiments/eval_mut/test_results/mfass_pangolin.txt

python baselines/spliceai.py data/vexseq/data.json experiments/eval_mut/test_results/vexseq_spliceai.txt
python baselines/spliceai.py data/mfass/data.json experiments/eval_mut/test_results/mfass_spliceai.txt

python baselines/run_mmsplice.py  data/vexseq/data.json py36hgfile/hg19.fa experiments/eval_mut/test_results/vexseq_mmsplice.txt
python baselines/run_mmsplice.py  data/mfass/data.json py36hgfile/hg19.fa experiments/eval_mut/test_results/mfass_mmsplice.txt


