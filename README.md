# Mhaemolytica_WGS
Pipeline for whole genome assembly and AMR analysis of Mannheimia haemolytica for Oxford nanopore reads. Requires input file with SampleName and Path of the directory with the fastq files. Outputs a htmlfile with AMR,Genotyping,Serotyping and virulence factors found in the sample. 

Usage
```
nextflow run main.nf --input sampleslist.csv --outdir Results_mannheimia_2 -profile docker --gsize 2.0M --trim_barcodes
```
```
Parameters:

--input		csv file with two columns with headers(SampleName,SamplePath).See samplelist.csv
--outdir	Output directory
--gsize 	genome size
--profile    	docker
optional
--trim_barcodes barcode and adapter trimming using porechop
``
