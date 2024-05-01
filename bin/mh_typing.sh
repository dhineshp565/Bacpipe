#!/usr/bin/env bash

# This script performs insilico pcr using genotyping, serotyping and virulence factor primers on the given fasta file and outputs the typing results in csv format.
# $1 = SampleName
# $2 = input fasta file
# $3 = genotyping primer file
# $4 = serotyping primer file
# $5 = virulence factor primer file

# genotyping using AdhesinG primers
cat $2|seqkit amplicon -p $3 --bed > $1_seqkit_geno.csv
	

if grep -Fq "AdhesinG" "$1_seqkit_geno.csv"
then
	echo "$1	Genotype-2" > $1_genotyping_results.csv
else
	echo "$1	Genotype-1" > $1_genotyping_results.csv
fi
sed -i '1i SampleName	Genotype' "$1_genotyping_results.csv"

# Serotyping A1,A2,And A6 using multiplex primers

cat $2|seqkit amplicon -p $4 --bed > $1_seqkit_sero.csv
	

if grep -Fq "Serotype-1-Hypothetical_protein_Hyp" "$1_seqkit_sero.csv"
then
	echo "$1	Serotype-A1" > $1_serotyping_results.csv
fi
if grep -Fq "Serotype-2-Core-2_I-Branching_enzyme_Core2" $1_seqkit_sero.csv 
then
	echo "$1	Serotype-A2" > $1_serotyping_results.csv
fi
if grep -Fq "Serotype-6-TupA-like_ATPgrasp_TupA" $1_seqkit_sero.csv 
then
	echo "$1	Serotype-A6" > $1_serotyping_results.csv	

fi
sed -i '1i SampleName	Serotype' "$1_serotyping_results.csv"

# Virulence factors using VF primers

cat $2|seqkit amplicon -p $5 --bed > $1_seqkit_VF.csv
sed -i '1i SampleName	Start	End	Virulence genes	0	Strand	Sequence' "$1_seqkit_VF.csv"
	