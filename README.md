# pgRNADesign: Designing paired gRNAs to knockout long non-coding RNAs (lncRNAs)

## Introduction

pgRNADesign is a python-based algorithm to design paired gRNAs for knocking out long non-coding RNAs (lncRNAs). pgRNADesign can be used when researchers need to design two gRNAs (pgRNAs) to delete large genomic regions of promoters or promoters+gene bodies for non-coding RNAs.

pgRNADesign can:

* design pgRNAs to target promoters (5kbp upstream of TSS) or promoters+lncRNA exons;
* assign a "barcode", or the identity of the first sgRNA in a pgRNA, to facilitate PCR and sequencing analysis;
* design pgRNAs to avoid the promoters and exons of coding genes;
* avoid some "blackout" regions (provided by the user) such that none of the pgRNAs cover these regions;

pgRNADesign is developed by Wei Li and Shiyou Zhu from X. Shirley Liu lab (Dana-Farber Cancer Institute and Harvard T.H. Chan School of Public Health) and Wensheng Wei lab (Peking University, China).


## Demo

There is a demo script (*demo.sh*) in the directory to illustrate the basic function of pgRNADesign. There is also a mini example (in *data* directory), providing all necessary files to run pgRNADesign. Simply run it on a terminal to get a basic idea of how it works:

```bash
./demo.sh
```

## Usage

The core program to generate gRNA pairs is *lncrna_pairs.py*.

lncrna_pairs.py: get sgRNA pairs (pgRNAs) for knocking out lncRNA promoters and bodies.

```
usage: ./lncrna_pairs.py 
       --sgrna SGRNA [--debug]
       --gene-annotation GENE_ANNOTATION 
       --gene-info GENE_INFO 
       [-h] [-o OUTPUT_FILE] [-n MAX_PAIRS] 
       [--blackout-region BLACKOUT_REGION]
       [--chromosome CHROMOSOME] 
       {promoter,gene,aavs1,nontarget}

```

**positional arguments:**

|{promoter,gene,aavs1,nontarget} | Types of pgRNAs: promoter (lncRNA promoters), gene (lncRNA promoter+gene body), aavs1 (aavs1 control), nontarget (nontarget control).

**optional arguments:**

|  -h, --help    |        show this help message and exit |
| -o OUTPUT_FILE, --output-file OUTPUT_FILE | Prefix of the output file, default results/sample1. |
| -n MAX_PAIRS, --max-pairs MAX_PAIRS | Maximum number of sgRNA pairs per gene, default -1 (no limits) | 
| --sgrna SGRNA  |      sgRNA sequence file (required). See *input and output files* section for more details. | 
| --debug       | Debug mode. |       
| --chromosome CHROMOSOME | chromosome to search. Default is to use all chromosomes. |
| --gene-annotation GENE_ANNOTATION | The annotation of genes and lncRNAs in GTF format (required). See *input and output files* section for more details. |
|  --gene-info GENE_INFO | The list of target lncRNAs (required). See *input and output files* section for more details. | 
| --blackout-region BLACKOUT_REGION | The list of blackout region (in bed file). If specified, none of the pgRNAs will span over the blackout region. |



There is another script, *getuniqpair.py*, that allows you to generate the maximum number of pgRNAs with barcodes. The barcode of a pgRNA is the first pair of its sgRNA pairs, which is unique across all pgRNAs in the library.

```
usage: ./getuniqpair.py {pgRNA list}
``` 

The script will generate a csv file containing all pgRNAs with barcodes. The barcode is the first sgRNA sequence of each pgRNA.

## Input and output files

### sgRNA list

The sgRNA list file is required by the *--sgrna* option, and contains the list of all putative sgRNAs. The file *data/sgrnalist.chr1.txt*  contains an example of how this file should look like:

```
chr     start   end     gene_symbol     strand  hg19_19_0       hg19_19_1       hg19_19_01      hg19_19_10      hg19_19_2       hg19_19_02      hg19_19_11      hg19_19_20      short_seq       long_seq        score
chr1    131652  131671  ENSG00000241860 -       5       2       1       1       2       1       1       0       GTCCTTGTCTGATTTGTTG     CACCGGATGTGGTCCTTGTCTGATTTGTTGGGGCTGCGTC
        0.026571
```

The first line is the header line, and the following lines are the sgRNA records, each line for one sgRNA. There are several required columns in this file, and you have to use the exact header to specify these columns (the order of the column doesn't matter).

**Required columns:**

| Header | Meaning |
|--------|---------|
|chr| The chromosome of the sgRNA|
|start| The start coordinate|
|end| The end coordinate|
|gene_symbol|The targeting lncRNA/gene IDs of the sgRNA. Note that if an sgRNA can target multiple lncRNAs/genes, use underscore `_` to separate multiple IDs. For example, `ENSG00000228327_ENSG00000237491`.|
|strand|The orientation of the sgRNA|
|short_seq|The sequence of the sgRNA|

Several columns are optional.

**Optional columns:**


| Header | Meaning |
|--------|---------|
|score|The efficiency score of the sgRNA, calculated by the [SSC](http://cistrome.dfci.harvard.edu/SSC) algorithm.|
|hg19_19_0|These are specificity measurements of the sgRNA. *hg19_19_0* means the number of hits of the sgRNA sequence mapped to the reference genome with 0 mismatch. |
|hg19_19_1| *hg19_19_0* means the number of hits of the sgRNA sequence mapped to the reference genome with 1 mismatch. |

The rest of the fields are not used in the current algorithm.


### gene information

This file is required by the *--gene-info* option, and contains the list of target lncRNAs. The file *data/scan_v8.transcriptinfo.txt* contains an example of how this file should look like:

```
ENSG00000230876 ENST00000414054_LINC00486       2       +       33050510        33151760
        33050510        33151760
ENSG00000230876 ENST00000435075_LINC00486       2       +       33107868        33151589
        33107868        33151589
```

The meanings of the columns are as follows. There is no header information, so you have to specify in the following order.

|Column|Meaning|
|------|-------|
|1     |lncRNA/gene ID|
|2     |Names of the lncRNA/gene|
|3     |The chromosome of the lncRNA/gene|
|4     |The orientation of the lncRNA/gene|
|5     |The start coordinate|
|6     |The end coordinate|
|7-8   |Currently not used in the program|



### gene annotation

Gene annotation file, required by the *--gene-annotation* option, provides coordinates of genes and lncRNAs. It follows the standard [GTF](http://useast.ensembl.org/info/website/upload/gff.html) format. The file *data/Homo_sapiens.CRCh37.75.chr1.gtf* contains an example of GTF annotation file (Homo Sapiens, GRCh37, chromosome 1).


### Output files

All output files by default are stored in the *results* directory. There are several output files generated:

| File name | Meanings |
|-----------|----------|
| design.bed|The pgRNA locations, in bed format|
| design.txt|All the information of pgRNAs (see below)|
|gene.txt| The summary of lncRNAs/genes, including the number of possible pairs and the number of pairs selected in the current design|
|goodsgrna.bed|All the sgRNAs that pass the quality check.|

An example of the *design.txt* is shown below.

```
ENSG00000228327_dist4588        4588    541.2   +GOODDIST;NDIST;+5UTR;+TSS;+TARGETEXON; chr1    713962  713981  TGCCCAGCTCCAGGCACCA     +       0.044485        chr1    718569  718588  CACGTGAACAGTTCTGGCC     +       0.107002        gene
```
 
Their meanings are specified as follows.

|Column|Content|Meaning|
|-----|-------|-------|
|1|ENSG00000228327_dist4588 |pgRNAs targeting the lncRNA|
|2|       4588 |The distance of the two sgRNAs in the pgRNA|
|3|   541.2 |The score of the pgRNA|
|4|  +GOODDIST;NDIST;+5UTR;+TSS;+TARGETEXON; |The code of the score of the pgRNA|
|5|chr1 |The chromosome of the first sgRNA|
|6|   713962 |The start coordinate of the first sgRNA | 
|7|713981 |The end coordinate |
|8| TGCCCAGCTCCAGGCACCA | The sequence of the first sgRNA |
|9|   + |The orientation of the first sgRNA |    
|10|  0.044485 |The efficiency of the first sgRNA |
|11|       chr1 | Columns 11-16 are the same as columns 5-10 for the second sgRNA|
|12|   718569 | |
|13| 718588 ||
|14| CACGTGAACAGTTCTGGCC  ||
|15|   + ||
|16 |     0.107002 ||
|17|       gene | The category of the pgRNA (gene,promoter, nontarget or aavs1)|



