
How two repoduce analysis of [state-of-the-art long-read overlapping tools comparative analysis](http://blog.pierre.marijon.fr/2018/04/04/State-of-the-art-long-reads-overlaper-compare.html)

# Acknowledgment

- @tctianchi for [pyvenn](https://github.com/tctianchi/pyvenn) script 

# Requirement

- matplotlib
- begins
- BioPython

# Dataset

## Nanopore 

Download [dataset](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/):
```
mkdir e_coli_nanopore
cd e_coli_nanopore
wget https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta
wget https://nanopore.s3.climb.ac.uk/MAP006-2_2D_pass.fasta
wget https://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta
wget https://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta
```

Merge reads in one file:
```
cat *pass.fasta > tmp.fasta
```

Set readname:
```
awk '/^>/{print ">" ++i; next}{print}' < tmp.fasta > merged.fasta
```

Clean:
```
rm *pass.fasta tmp.fasta
cd ..
```

## Pacbio

Download [dataset](https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly):
```
mkdir e_coli_pacbio
cd e_coli_pacbio
wget https://s3.amazonaws.com/files.pacb.com/datasets/secondary-analysis/e-coli-k12-P6C4/p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz
tar xvfz p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz
dextract E01_1/Analysis_Results/m141013_011508_sherri_c100709962550000001823135904221533_s1_p0.*.bax.h5
```


Merge reads:
```
cat E01_1/Analysis_Results/m141013_011508_sherri_c100709962550000001823135904221533_s1_p0*.fasta > tmp.fasta
```

Set readname:
```
awk '/^>/{print ">" ++i; next}{print}' < tmp.fasta > merged.fasta
```

Clean:
```
rm -rf E01_1
rm tmp.fasta
cd ..
```

## Synthetic

Create a random sequence (20.000, with equiprobability) and store it in `synthetic/genome.fasta`.
Create perfect read (length 10.000) begin at each 1000 base :
```
cd synthetic
./perfect_read.py genome.fasta merged.fasta 1000 10000
```

# Run mapper

## Nanopore

```
graphmap owler -r e_coli_nanopore/merged.fasta -d e_coli_nanopore/merged.fasta -o graphmap_nanopore.mhap
hisea --ref e_coli_nanopore/merged.fasta --self --threads 10 > hisea_nanopore.hisea
mhap -Xmx96g --num-threads 15 -s e_coli_nanopore/merged.fasta -k 14 --num-hashes 4844 --num-min-matches 3 --threshold 0.02 > mhap_nanopore.mhap
minimap -x ava10k -k 14 e_coli_nanopore/merged.fasta e_coli_nanopore/merged.fasta > minimap_nanopore.paf
minimap2 -x ava-ont -k 14 e_coli_nanopore/merged.fasta e_coli_nanopore/merged.fasta > minimap2_nanopore.paf
```

## Pacbio

```
graphmap owler -r e_coli_pacbio/merged.fasta -d e_coli_pacbio/merged.fasta -o graphmap_pacbio.mhap
hisea --ref e_coli_pacbio/merged.fasta --self --threads 10 > hisea_pacbio.hisea
mhap -Xmx96g --num-threads 15 -s e_coli_pacbio/merged.fasta -k 14 --num-hashes 4844 --num-min-matches 3 --threshold 0.02 > mhap_pacbio.mhap
minimap -x ava10k -k 14 e_coli_pacbio/merged.fasta e_coli_pacbio/merged.fasta > minimap_pacbio.paf
minimap2 -x ava-pb -k 14 e_coli_pacbio/merged.fasta e_coli_pacbio/merged.fasta > minimap2_pacbio.paf
```

## Synthetic

```
graphmap owler -r synthetic/merged.fasta -d synthetic/merged.fasta -o graphmap_synthetic.mhap
hisea --ref synthetic/merged.fasta --self --threads 10 > hisea_synthetic.hisea
mhap -Xmx96g --num-threads 15 -s synthetic/merged.fasta -k 14 --num-hashes 4844 --num-min-matches 3 --threshold 0.02 > mhap_synthetic.mhap
minimap -x ava10k -k 14 synthetic/merged.fasta synthetic/merged.fasta > minimap_synthetic.paf
minimap2 -x ava-ont -k 14 synthetic/merged.fasta synthetic/merged.fasta > minimap2_synthetic.paf
```

# Generate figure

## run

```
./venn_script.py figure.png minimap2_pacbio.paf hisea_pacbio.hisea
```
