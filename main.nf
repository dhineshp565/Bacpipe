#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// make csv file with headers from the given input

process make_csv {
	publishDir "${params.out_dir}"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	makecsv.sh ${fastq_input}

	"""

}


//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.out_dir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}")
	
	shell:

	"""
	count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]]
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
					
		else
			count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
			if [[ "\${count}" != "0" ]]
			then
				cat ${SamplePath}/*.fastq|gzip > ${SampleName}.fastq.gz
				
			fi
		fi
	"""
	
}
process porechop {
	label "high"
	publishDir "${params.out_dir}/trimmed",mode:"copy"
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
    publishDir "${params.out_dir}/Assembly",mode:"copy"
    input:
    tuple val(SampleName),path(SamplePath)
	val(medaka_model)
    output:
    val(SampleName),emit:sample
	tuple val(SampleName),path("${SampleName}_flye.fasta"),emit:assembly
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --model ${medaka_model} --gsize 2.4M --nanohq --medaka 1
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
	label "low"
	publishDir "${params.out_dir}/abricate",mode:"copy"
	input:
	tuple val(SampleName),path(fastafile)
	path (custom_description)
	path (dbdir)
	output:
	val(SampleName),emit:sample
	path("${SampleName}_CARD.csv"),emit:card
	path("${SampleName}_CARD.txt"),emit:cardtxt
	path("${SampleName}_VF.txt"),emit:vf
	path("${SampleName}_ICE.txt"),emit:ice
	path("${SampleName}_type.txt"),emit:typing
	script:
	"""
	MH_abricate.sh ${fastafile} ${SampleName} ${custom_description} ${dbdir}
	"""
}

process abricate_summary {
	label "low"
	publishDir "${params.out_dir}/abricate_summary",mode:"copy"
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
	label "low"
	publishDir "${params.out_dir}/seqkit",mode:"copy"
	input:
	tuple val(SampleName),path(fasta)
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
	mh_typing.sh ${SampleName} ${fasta} ${geno_primerfile} ${sero_primerfile} ${vf_primerfile}
	"""
}
process summarize_csv {
	label "low"
	publishDir "${params.out_dir}/summary",mode:"copy"
	input:
	path(card)
	path(geno)
	path(sero)
	path(vf)
	path(ice)
	path (type)
	output:
	path ("typing_results")
	script:
	"""
	mkdir typing_results
	cp ${card} typing_results/
	cp ${geno} typing_results/
	cp ${sero} typing_results/
	cp ${vf} typing_results/
	cp ${ice} typing_results/
	cp ${type} typing_results/
	"""


}
process make_report {
	label "low"
	publishDir "${params.out_dir}",mode:"copy"
	input:
	path(rmdfile)
	path(types)
	path(csv)

	output:
	path("MH_results_report.html")

	script:

	"""
	
	cp ${rmdfile} rmdfile_copy.Rmd
	cp ${types}/* ./
	cp ${csv} sample.csv

	Rscript -e 'rmarkdown::render(input="rmdfile_copy.Rmd",params=list(csv="sample.csv"),output_file="MH_results_report.html")'
	"""

}


workflow {
    data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	// Merge fastq files for each sample

	// based on the optional argument trim barcodes using porechop and assemble using dragonflye
    if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out,params.medaka_model) 
	} else {
        dragonflye(merge_fastq.out,params.medaka_model)           
    }
	abricate_db="${baseDir}/Mhaemolytica_db"
	AMR_description=file("${baseDir}/AMR_descriptions.txt")
	//AMR gene finding using abricate and CARD database
	abricate(dragonflye.out.assembly,AMR_description,abricate_db)
	card=abricate.out.card.collect()
	abricate_summary(card)
	//custom_description(abricate.out.sample,abricate.out.card,AMR_description)

	// genotyping serotyping and virulence typing using seqkit amplicon
	geno_primerfile=file("${baseDir}/MH_genotyping_primers_geno.tsv")
	sero_primerfile=file("${baseDir}/Mannheimia_serotyping_primers.tsv")
	vf_primerfile=file("${baseDir}/MH_VF_primers.tsv")
	seqkit_typing(dragonflye.out.assembly,geno_primerfile,sero_primerfile,vf_primerfile)
	rmdfile=file("${baseDir}/MH_tabbed_report.Rmd")
	// summarize all AMR abd typing data
	summarize_csv(abricate.out.cardtxt.collect(),seqkit_typing.out.geno.collect(),seqkit_typing.out.sero.collect(),abricate.out.vf.collect(),abricate.out.ice.collect(),abricate.out.typing.collect())
	// make rmarkdown report from the the summarised files
	make_report(rmdfile,summarize_csv.out,make_csv.out)
}
