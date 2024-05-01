#!/usr/bin/env bash

# $1= fasta file,$2=SampleName,$3=custom_description,$4=dbdir

# Perform abricate analysis using the card database and add a custom description
abricate --db card $1 > $2_CARD.csv

# Check if no AMR genes found
numblines=$(< "$2_CARD.csv" wc -l)

if [ ${numblines} -gt 1 ]
then
	# Run AMR_description.py script
	AMR_description.py "$3" "$2_CARD.csv" "$2"
else 
	echo "No AMR Genes Found" >> "$2_CARD.csv"
	cp "$2_CARD.csv" "$2_CARD.txt"
fi

# Perform abricate analysis using the custom Mhaemolytica_ICEBERG database
abricate --datadir $4 --db Mhaemolytica_ICEBERG -minid 80  -mincov 70 --quiet $1 1> $2_ICE.txt

if [ $(< "$2_ICE.txt" wc -l) -eq 1 ]
then
	echo "No ICE found" >> "$2_ICE.txt"
fi

# Perform abricate analysis using the custom Mhaemolytica_typedb database
abricate --datadir $4 --db Mhaemolytica_typedb -minid 80  -mincov 70 --quiet $1 1> $2_typing.csv

# Create header for $2_type.txt file
echo -e "SAMPLE\tGENOTYPE\tSEROTYPE\tDESCRIPTION" > $2_type.txt

# Check for specific genotypes and serotypes and write to Samplename_type.txt file
if grep -Fq "Genotype-2" "$2_typing.csv" && grep -Fq "Serotype-1" "$2_typing.csv";then
	echo -e "$2\tGenotype-2\tSerotype-1\tGenotype-2(Serotype-A1 & A6) strains are predominantly isolated from lower respiratory tract of cattle with  Bovine respiratory diseases complex" >> $2_type.txt
fi
if grep -Fq "Genotype-2" "$2_typing.csv" && grep -Fq "Serotype-6" "$2_typing.csv"
then
	echo -e "$2\tGenotype-2\tSerotype-6\tGenotype-2(Serotype-A1 & A6) strains are predominantly isolated from lower respiratory tract of cattle with  Bovine respiratory disease complex" >> $2_type.txt
fi
if grep -vFq "AdhesinG" "$2_typing.csv" && grep -Fq "Serotype-2" $2_typing.csv
then
	echo -e "$2\tGenotype-1\tSerotype-2\tGenotype-1(Serotype-A2) strains are frequently isolated fromfrom the upper respiratory tract of cattle without signs of Bovine respiratory disease complex" >> $2_type.txt
fi

# Perform abricate analysis using the custom Mhaemolytica_vfdb database
abricate --datadir $4 --db Mhaemolytica_vfdb -minid 80  -mincov 70 --quiet $1 1> $2_VF.txt

# Modify the filenames in the output files
sed -i 's/_flye.fasta//g' *.txt
sed -i 's/_/ /g' *.txt