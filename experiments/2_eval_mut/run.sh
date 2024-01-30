rm -r experiments/2_eval_mut/test_results
source /home/xuc0d/anaconda3/bin/activate openfold
python -m Tools.generate_mutdata data/VexSeq_snps2exon_ref_dataset.txt data/vexseq/ hg19  #File, SavePath, species =
python -m Tools.generate_mutdata data/MFASS.bed data/mfass hg19
CUDA_VISIBLE_DEVICES=1 python main_mut.py -c experiments/2_eval_mut/RefSplice_mut_config
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/vexseq/data.json experiments/2_eval_mut/test_results/vexseq_pangolin.txt
CUDA_VISIBLE_DEVICES=1 python baselines/pangolin.py data/mfass/data.json experiments/2_eval_mut/test_results/mfass_pangolin.txt

source /home/xuc0d/anaconda3/bin/activate spliceai
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/vexseq/data.json experiments/2_eval_mut/test_results/vexseq_spliceai.txt
CUDA_VISIBLE_DEVICES=1 python baselines/spliceai.py data/mfass/data.json experiments/2_eval_mut/test_results/mfass_spliceai.txt

source /home/xuc0d/anaconda3/bin/activate mmsplice
mkdir tmp_vcf
python baselines/run_mmsplice.py  data/vexseq/data.json py36hgfile/hg19.fa experiments/2_eval_mut/test_results/vexseq_mmsplice.txt
python baselines/run_mmsplice.py  data/mfass/data.json py36hgfile/hg19.fa experiments/2_eval_mut/test_results/mfass_mmsplice.txt


