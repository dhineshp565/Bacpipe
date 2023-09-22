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
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --model r1041_e82_400bps_sup_g615 --gsize ${gsize} --nanohq --medaka 1
    mv "${SampleName}_assembly"/flye.fasta "${SampleName}"_flye.fasta
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye.fasta"
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
	path("${SampleName}_resfinder.csv"),emit:resfinder
	path("${SampleName}_CARD.csv"),emit:card
	path("${SampleName}_megares.csv"),emit:megares
	path("${SampleName}_ncbi.csv"),emit:ncbi
	path("${SampleName}_argannot.csv"),emit:argannot
	path("${SampleName}_vfdb.csv"),emit:vfdb
	script:
	"""
	abricate --db resfinder ${fastafile} > ${SampleName}_resfinder.csv
	abricate --db card ${fastafile} > ${SampleName}_CARD.csv
	abricate --db megares ${fastafile} > ${SampleName}_megares.csv
	abricate --db ncbi ${fastafile} > ${SampleName}_ncbi.csv
	abricate --db argannot ${fastafile} > ${SampleName}_argannot.csv
	abricate --db vfdb ${fastafile} > ${SampleName}_vfdb.csv
	
	"""
}
process abricate_summary {
	label "high"
	publishDir "${params.outdir}/abricate_summary",mode:"copy"
	input:
	path(resfinder)
	path(card)
	path(megares)
	path(ncbi)
	path(argannot)
	path(vfdb)
	output:
	path("All_resfinder.csv")
	path("All_CARD.csv")
	path("All_megares.csv")
	path("All_ncbi.csv")
	path("All_argannot.csv")
	path("All_vfdb.csv")
	script:
	"""
	abricate --summary ${resfinder} > All_resfinder.csv
	abricate --summary ${card}> All_CARD.csv
	abricate --summary ${megares} > All_megares.csv
	abricate --summary ${ncbi} > All_ncbi.csv
	abricate --summary ${argannot} > All_argannot.csv
	abricate --summary ${vfdb} > All_vfdb.csv
	
	"""
}
process mlst {
	label "high"
	publishDir "${params.outdir}/mlst",mode:"copy"
	input:
	val(SampleName)
	path(fastafile)
	output:
	path("${SampleName}_MLST.csv")
	script:
	"""
	mlst ${fastafile} > ${SampleName}_MLST.csv
	"""
}
process rgi_amr {
	label "high"
	publishDir "${params.outdir}/rgi_amr",mode:"copy"
	input:
	val(SampleName)
	path(fastafile)
	output:
	path("${SampleName}_rgi.txt")
	script:
	"""
	rgi main --input_sequence ${fastafile} --output_file ${SampleName}_rgi -t contig --include_loose
	"""

}
process prokka {
	label "medium"
	publishDir "${params.outdir}/prokka",mode: "copy"
	input:
	val(SampleName)
	path(fastafile)
	output:
	path ("${SampleName}_prokka")
	path("${SampleName}_prokka.gff"),emit:gff
	script:
	"""
	prokka --outdir ${SampleName}_prokka --locustag ${SampleName} --prefix ${SampleName} ${fastafile}
	cat "${SampleName}_prokka"/*.gff > ${SampleName}_prokka.gff

	"""
}
process roary{
	label "high"
	publishDir "${params.outdir}/roary",mode:"copy"
	input:
	path(gff)
	output:
	path("Roary")
	script:
	"""
	roary -f Roary ${gff}
	"""

}

process genotyping_minimap {
	label "high"
	publishDir "${params.outdir}/minimap2",mode:"copy"
	input:
	val(SampleName)
	path(fastafile)
	file(adhesinG)
	output:
	val(SampleName)
	path("${SampleName}.sam")
	script:
	"""
	minimap2 -a ${adhesinG} ${fastafile} > ${SampleName}.sam
	"""
}
process genotyping_samtools {
	label "high"
	publishDir "${params.outdir}/samtools",mode:"copy"
	input:
	val(SampleName)
	path(samfile)
	output:
	val(SampleName)
	path("${SampleName}_mappedreads.csv")
	script:
	"""
	samtools view -b -F 256 "${samfile}" > ${SampleName}_mapped.bam
	samtools index ${SampleName}_mapped.bam
	samtools idxstats ${SampleName}_mapped.bam > ${SampleName}_idxstats.csv
	awk '{if (\$3!=0) print \$1,\$2,\$3}' "${SampleName}_idxstats.csv" > ${SampleName}_mappedreads.csv
	"""
}

process genotyping_seqkit{
	label "high"
	publishDir "${params.outdir}/seqkit",mode:"copy"
	input:
	val(SampleName)
	file(fasta)
	path(primerfile)
	output:
	val(SampleName)
	path("${SampleName}_seqkit.csv")
	path("${SampleName}_genotyping_results.csv")
	script:
	"""
	cat ${fasta}|seqkit amplicon -p ${primerfile} --bed > ${SampleName}_seqkit.csv
	

	if grep -Fq AdhesinG ${SampleName}_seqkit.csv 
	then
		echo "${SampleName}	Genotype-2" > ${SampleName}_genotyping_results.csv
	else
		echo "${SampleName}	Genotype-1" > ${SampleName}_genotyping_results.csv
	fi
	sed -i '1i SampleName	Type' "${SampleName}_genotyping_results.csv"
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
	resfind=abricate.out.resfinder.collect()
	car=abricate.out.card.collect()
	mega=abricate.out.megares.collect()
	ncb=abricate.out.ncbi.collect()
	arga=abricate.out.argannot.collect()
	vfd=abricate.out.vfdb.collect()
	abricate_summary(resfind,car,mega,ncb,arga,vfd)
	mlst(dragonflye.out.sample,dragonflye.out.assembly)
	rgi_amr(dragonflye.out.sample,dragonflye.out.assembly)
	prokka(dragonflye.out.sample,dragonflye.out.assembly)
	gff=prokka.out.gff.collect()
	roary(gff)
	adh=file("${baseDir}/AdhesinG_genotype2_Mannheimiah.fasta")
	genotyping_minimap(dragonflye.out.sample,dragonflye.out.assembly,adh)
	genotyping_samtools(genotyping_minimap.out)
	primerfile=file("${baseDir}/MH_genotypping_primers_geno.tsv")
	genotyping_seqkit(dragonflye.out.sample,dragonflye.out.assembly,primerfile)
}