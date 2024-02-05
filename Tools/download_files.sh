mkdir fafiles
cd fafiles
for name in hg19  hg38 mm10 panTro5 rheMac10 rn6 susScr11 bosTau9 calJac4
do
wget https://hgdownload.cse.ucsc.edu/goldenpath/${name}/bigZips/${name}.fa.gz
done
gzip -d *.gz
cd ..
cp -r fafiles py36hgfile
mkdir data/Chains
cd data/Chains
for name in hg19ToBosTau9.over.chain hg19ToMm10.over.chain hg19ToPanTro5.over.chain  hg19ToRheMac10.over.chain hg19ToRn6.over.chain hg19ToSusScr11.over.chain
do 
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/${name}.gz
done
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gzip -d *.gz

cd ../..