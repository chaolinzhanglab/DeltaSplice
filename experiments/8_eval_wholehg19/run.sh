source /home/xuc0d/anaconda3/bin/activate openfold
for n in hg19.fa bosTau9.fa calJac4.fa mm10.fa panTro5.fa rheMac10.fa rn6.fa susScr11.fa
do 
    CUDA_VISIBLE_DEVICES=0 python -m experiments.8_eval_wholehg19.main ${n}
done