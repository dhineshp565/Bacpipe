# Mhaemolytica_WGS
Pipeline for whole genome assembly and AMR analysis of Mannheimia haemolytica for Oxford nanopore reads. Requires input file with SampleName and Path of the directory with the fastq files. Outputs a htmlfile with AMR,Genotyping,Serotyping and virulence factors found in the sample. 

Usage
```
nextflow run main.nf --input sampleslist.csv --out_dir Results_mannheimia_2 --trim_barcodes
```
```
Parameters:

--input		csv file with two columns with headers(SampleName,SamplePath).See samplelist.csv
--out_dir	Output directory

optional
--trim_barcodes barcode and adapter trimming using porechop
``
