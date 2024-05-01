# Mhaemolytica_WGS
Pipeline for whole genome assembly and AMR analysis of Mannheimia haemolytica for Oxford nanopore reads. Requires input Path of the directory with the fastq files. Outputs a htmlfile with AMR,Genotyping,Serotyping and virulence factors found in the sample. 

Usage
```
nextflow run main.nf --input /data/fastqdir --out_dir Results_mannheimia_2 --trim_barcodes
```
```
Parameters:

--input		Path to directory containing sub-directory with fastq files
--out_dir	Path to Output directory

optional
--trim_barcodes barcode and adapter trimming using porechop
--medaka_model Add medaka model for polishing. The default model is r1041_e82_400bps_sup_g615.
``
