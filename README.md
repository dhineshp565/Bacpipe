# Mhaemolytica_WGS
Pipeline for whole genome assembly and AMR analysis of Mannheimia haemolytica fro Oxford nanopore reads. 


```
nextflow run main.nf \
		 --input /data/Dhinesh/Mannehimia/sampleslist.csv \
		 --outdir Results_mannheimia_2 \
		 --trim_barcodes \
		 -profile docker \
		 --gsize 2.0M \
		 -resume \
```
