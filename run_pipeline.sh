#! /bin/sh
gff3=$1
fasta=$2
query=$3
length=$4

rm -rf ./tmp/*
python3 --version
python3 query_preProc.py $query

for fl in $(ls ./tmp/*txt | sed 's/.txt//g')
do
    prefix=$(basename $fl)
    python3 GFF3seqExtra.py $gff3 ${fl}.txt $length
    bedtools getfasta -fi $fasta -bed geneAppend.bed -name -fo ./res/${prefix}_length_${length}.fasta
    rm -rf geneAppend.bed
done
rm -rf geneAppend.bed
#rm -rf ./tmp/*
