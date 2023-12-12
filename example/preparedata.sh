mkdir data
mkdir data/Fa
cd data/Fa
for name in hg19 susScr11 mm10 rheMac10  rn6 panTro5 bosTau9
do
wget http://hgdownload.cse.ucsc.edu/goldenpath/${name}/bigZips/${name}.fa.gz
done
gzip -d *

cd ..
mkdir WGFile
cd WGFile
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals.phyloP46way.bw
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phyloP60way/mm10.60way.phyloP60way.bw

cd ..
mkdir PretrainedModel
mkdir TestResultWCons
mkdir TestResultWOCons
cd PretrainedModel
## todo download pretrained model
cd ..
cd example
## todo download vexseq/table file
cd ..
cp example/config.py config.py
python -m Tools.annotate_gene
python -m Tools.generate_data
python -m Tools.generate_mutdata example/VexSeq_snps2exon_ref_dataset.txt data/VexSeq/ hg19

