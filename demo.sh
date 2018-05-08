#!/bin/bash

./lncrna_pairs.py  --sgrna data/sgrnalist.chr1.txt --gene-annotation data/Homo_sapiens.CRCh37.75.chr1.gtf --gene-info data/scan_v8.transcriptinfo.txt  promoter

./lncrna_pairs.py  --sgrna data/sgrnalist.chr1.txt --gene-annotation data/Homo_sapiens.CRCh37.75.chr1.gtf --gene-info data/scan_v8.transcriptinfo.txt  gene

# generating AAVS1 controls

./lncrna_pairs.py  --sgrna data/AAVS1_from_wglibrary.txt --gene-annotation data/Homo_sapiens.CRCh37.75.chr1.gtf --gene-info data/scan_v8.transcriptinfo.txt  aavs1

# generating non-target controls
./lncrna_pairs.py  --sgrna data/nontargeting.txt --gene-annotation data/Homo_sapiens.CRCh37.75.chr1.gtf --gene-info data/scan_v8.transcriptinfo.txt  nontarget


# merging results
cat results/*.design.txt > results/merged.txt
cat results/*.design.bed > results/merged.bed

# getting unique pairs and generating barcodes
./getuniqpair.py  results/merged.txt


