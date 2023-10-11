#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.outdir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}")
	
	shell:
	"""
	count=\$(ls -1 $SamplePath/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]];
		then
			cat $SamplePath/*.fastq.gz > ${SampleName}.fastq.gz
		
		else
			cat $SamplePath/*.fastq > ${SampleName}.fastq
		fi
	"""
}
process porechop {
	label "medium"
	publishDir "${params.outdir}/trimmed",mode:"copy",overwrite: false
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
	"""
}
process dragonflye {
    label "high"
    publishDir "${params.outdir}/Assembly",mode:"copy"
    input:
    tuple val(SampleName),path(SamplePath)
    val(gsize)
    output:
    val(SampleName),emit:sample
	path("${SampleName}_flye.fasta"),emit:assembly
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --model r1041_e82_400bps_sup_g615 --gsize 2.4M --nanohq --medaka 1
    # rename fasta file with samplename
    mv "${SampleName}_assembly"/flye.fasta "${SampleName}"_flye.fasta
    # rename fasta header with samplename
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye.fasta"
     # rename flyeinfo file and contnents
    mv "${SampleName}_assembly"/flye-info.txt "${SampleName}"_flye-info.txt
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye-info.txt"
    """
}
process abricate {
	label "high"
	publishDir "${params.outdir}/abricate",mode:"copy"
	input:
	val(SampleName)
	path(fastafile)
	output:
	val(SampleName),emit:sample
	path("${SampleName}_CARD.csv"),emit:card
	script:
	"""
	abricate --db card ${fastafile} > ${SampleName}_CARD.csv
	"""
}
process abricate_summary {
	label "high"
	publishDir "${params.outdir}/abricate_summary",mode:"copy"
	input:
	path(card)
	output:
	path("All_CARD.csv"),emit:card
	script:
	"""
	abricate --summary ${card}> All_CARD.csv
	# Change header of the summary file
	sed -i 's/#FILE/SampleName/g' "All_CARD.csv"
	sed -i 's/NUM_FOUND/No\sof\sAMR\sgenes/g' "All_CARD.csv"
	"""
}

process seqkit_typing {
	label "high"
	publishDir "${params.outdir}/seqkit",mode:"copy"
	input:
	val(SampleName)
	file(fasta)
	path(geno_primerfile)
	path (sero_primerfile)
	path(vf_primerfile)
	output:
	val(SampleName),emit:sample
	path("${SampleName}_seqkit_geno.csv")
	path("${SampleName}_genotyping_results.csv"),emit:geno
	path("${SampleName}_seqkit_sero.csv")
	path("${SampleName}_serotyping_results.csv"),emit:sero
	path("${SampleName}_seqkit_VF.csv"),emit:vf

	script:
	"""
	# genotyping using AdhesinG primers
	cat ${fasta}|seqkit amplicon -p ${geno_primerfile} --bed > ${SampleName}_seqkit_geno.csv
	

	if grep -Fq AdhesinG "${SampleName}_seqkit_geno.csv"
	then
		echo "${SampleName}	Genotype-2" > ${SampleName}_genotyping_results.csv
	else
		echo "${SampleName}	Genotype-1" > ${SampleName}_genotyping_results.csv
	fi
	sed -i '1i SampleName	Type' "${SampleName}_genotyping_results.csv"

	# Serotyping A1,A2,And A6 using multiplex primers

	cat ${fasta}|seqkit amplicon -p ${sero_primerfile} --bed > ${SampleName}_seqkit_sero.csv
	

	if grep -Fq "Serotype-1-Hypothetical_protein_Hyp" "${SampleName}_seqkit_sero.csv"
	then
		echo "${SampleName}	Serotype-A1" > ${SampleName}_serotyping_results.csv
	fi
	if grep -Fq "Serotype-2-Core-2_I-Branching_enzyme_Core2" ${SampleName}_seqkit_sero.csv 
	then
		echo "${SampleName}	Serotype-A2" > ${SampleName}_serotyping_results.csv
	fi
	if grep -Fq "Serotype-6-TupA-like_ATPgrasp_TupA" ${SampleName}_seqkit_sero.csv 
	then
		echo "${SampleName}	Serotype-A6" > ${SampleName}_serotyping_results.csv	

	fi
	sed -i '1i SampleName	Type' "${SampleName}_serotyping_results.csv"

	# Virulence factors using VF primers

	cat ${fasta}|seqkit amplicon -p ${vf_primerfile} --bed > ${SampleName}_seqkit_VF.csv
	"""
}
process summarize_csv {
	label "low"
	publishDir "${params.outdir}/summary",mode:"copy"
	input:
	path(card)
	path(geno)
	path(sero)
	path(vf)
	output:
	path("card_all.csv")
	path("geno_all.csv")
	path("sero_all.csv")
	path("vf_all.csv")
	script:
	"""
	awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' ${card} > card_all.csv
	awk 'FNR==1 && NR!=1 { while (/^SampleName/) getline; } 1 {print}' ${geno} > geno_all.csv
	awk 'FNR==1 && NR!=1 { while (/^SampleName/) getline; } 1 {print}' ${sero} > sero_all.csv
	cat ${vf} > vf_all.csv
	sed -i '1i SampleName	Start	End	Virulence_genes	0	Strand	Sequence' "vf_all.csv"
	"""


}
process make_report {
	label "high"
	publishDir "${params.outdir}/reports",mode:"copy"
	input:
	path(rmdfile)
	path(amr_summary)
	path(amr)
	path(geno)
	path(sero)
	path(vf)

	output:
	path("MH_report.html")

	script:

	"""
	cp ${rmdfile} rmdfile_copy.Rmd
	cp ${amr_summary} amr_summary.csv
	cp ${amr} amr.csv
	cp ${geno} geno.csv
	cp ${sero} sero.csv
	cp ${vf} vf.csv

	Rscript -e 'rmarkdown::render(input="rmdfile_copy.Rmd",params=list(card ="amr.csv",card_summary="amr_summary.csv",geno = "geno.csv",sero= "sero.csv",vf="vf.csv"),output_file="MH_report.html")'

	"""

}


workflow {
    data=Channel
	.fromPath(params.input)
	.splitCsv(header:true)
    .map { row-> tuple(row.SampleName,row.SamplePath) }

     merge_fastq(data)
    if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out,params.gsize) 
	} else {
        dragonflye(merge_fastq.out,params.gsize)           
    }
	abricate(dragonflye.out.sample,dragonflye.out.assembly)
	card=abricate.out.card.collect()
	abricate_summary(card)
	geno_primerfile=file("${baseDir}/MH_genotyping_primers_geno.tsv")
	sero_primerfile=file("${baseDir}/Mannheimia_serotyping_primers.tsv")
	vf_primerfile=file("${baseDir}/MH_VF_primers.tsv")
	seqkit_typing(dragonflye.out.sample,dragonflye.out.assembly,geno_primerfile,sero_primerfile,vf_primerfile)
	rmdfile=file("${baseDir}/MH_report.Rmd")
	summarize_csv(abricate.out.card.collect(),seqkit_typing.out.geno.collect(),seqkit_typing.out.sero.collect(),seqkit_typing.out.vf.collect())
	make_report(rmdfile,abricate_summary.out.card,summarize_csv.out)
}
