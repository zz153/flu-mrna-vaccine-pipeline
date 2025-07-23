#concatenate the sequences
cat ../data/raw/HA_H1.fa ../data/raw/HA_H3.fa ../data/raw/HA_H5.fa ../data/raw/HA_H7.fa > ../data/processed/combined_HA.fa

#count the number of sequences
grep "^>" ../data/processed/combined_HA.fa | wc -l

#run msa using mafft tool
mafft --auto --thread -1 ../data/processed/combined_HA.fa > ../data/processed/combined_HA_aligned.fa
