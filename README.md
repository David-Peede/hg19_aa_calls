# Hg19 Ancestral Allele Calls

This repository contains the code to convert the chimp (panTro6) multiple sequence alignment and Enredo-Pecan-Ortheus (EPO) ancestral sequence alignment to an all-sites VCF file. If you don't want to run the code yourself, all of the processed VCF files and tables can be downloaded from my [Dropbox](https://www.dropbox.com/scl/fo/t4plqfno13fckii979hry/AGnubfwGr_81Fpn_USXRMkU?rlkey=n4hx7wa8ieelz7sfew4s213ke&st=a2nce2uw&dl=0). If you would like to run this code yourself, follow the instructions below.



## `data`

All of the data is publicly available for download.

```bash
# Download the latest Hg19 reference genome from UCSC.
wget -P ./data/hg19 https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz

# Download the panTro6 MAF file.
wget -P ./data/panTro https://hgdownload.soe.ucsc.edu/goldenPath/hg19/vsPanTro6/hg19.panTro6.synNet.maf.gz

# Download the Enredo-Pecan-Ortheus (EPO) ancestral sequences.
wget -P ./data https://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
# Extract and then delete the EPO tar file.
cd ./data
tar -xvjf homo_sapiens_ancestor_GRCh37_e71.tar.bz2
rm homo_sapiens_ancestor_GRCh37_e71.tar.bz2
```


## `aa_calls`

This directory containes two subdirectories `tables` and `vcfs`. It should be noted that the code to proccess the MAF files in this repo was inspired by Simon Martin's [`genomics_general` repo](https://github.com/simonhmartin/genomics_general) that I optimized for my own specific use.


### `tables` 

This directory contains gzipped CSV files with the Hg19 reference, EPO, and panTro6 alleles for comparison.

```bash
# Generate tables with the Hg19 reference sequence, EPO ancestral sequence, and panTro6 sequence.
for CHR in {1..22} X Y; do
    python ./aa_tools/hg19_epo_panTro6_table.py -c ${CHR}
done
```

### `vcfs`

This directory contains the all sites VCF files for the EPO and panTro6 alleles—note you will need [`tabix`](https://www.htslib.org/doc/tabix.html) to bgzip the VCF files.

```bash
# Generate VCF files for the panTro6 ancestral sequence.
for CHR in {1..22} X Y; do
    python ./aa_tools/hg19_panTro6_vcf.py -c ${CHR} | bgzip > ./aa_calls/vcfs/hg19_panTro6_chr${CHR}.vcf.gz
done

# Generate VCF files for the EPO ancestral sequence.
for CHR in {1..22} X Y; do
    python ./aa_tools/hg19_epo_vcf.py -c ${CHR} | bgzip > ./aa_calls/vcfs/hg19_epo_chr${CHR}.vcf.gz
done
```
