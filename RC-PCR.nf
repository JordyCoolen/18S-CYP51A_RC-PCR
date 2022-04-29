#!/usr/bin/env nextflow
params.threads = 4
params.outDir = "./output"
params.reads = "$baseDir/test/test_OUT01_R{1,2}.fastq.gz"
params.UMILEN = 6
params.minreadlength = 100
params.database = "18S"
params.snp_database = "CYP51A"
params.abricate = true
params.visuals = true

NXF_OPTS="-Xms8g -Xmx8g"
_JAVA_OPTIONS="-Xms8g -Xmx8g"

// Parsing the input parameters
outDir           = "$params.outDir"
threads          = "$params.threads"
UMILEN           = "$params.UMILEN"
minreadlength    = "$params.minreadlength"
database         = "$params.database"
snp_database     = "$params.snp_database"
primerfile       = "${baseDir}/db/${database}/primers/${database}_primers.fasta"
ptrimfile        = "${baseDir}/db/${database}/primers/all_amplicon.txt"
bedfile          = "${baseDir}/db/CYP51A/primers/CYP51A_bedpe_primers.bed"
KMAdb            = "${baseDir}/db/${database}/KMA/${database}"
blastdbpath      = "${baseDir}/db/${database}/"
STARfasta        = "${baseDir}/db/${snp_database}/STAR/${snp_database}.fa"
def samplename   = file("$params.reads").simpleName[0].split('_')[0]
abricate         = "$params.abricate"
visuals          = "$params.visuals"

// Tools paths and command prefixes
reporter         = "$baseDir/scripts/final_report.py"
vcf2table        = "$baseDir/scripts/vcf2table.py"

// Parsing the input parameters
outDir           = "$params.outDir"+"/${samplename}"
threads          = "$params.threads"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .into  { reads_ch1; reads_ch2 }

log.info """

NEXTFLOW 18S/CYP51A RC-PCR V0.4
================================
samplename : $samplename
reads      : $params.reads
outDir     : $params.outDir
codeBase   : $baseDir
threads    : $params.threads
abricate   : $params.abricate
UMILEN     : $params.UMILEN
minreadlenth : $params.minreadlength
database   : $params.database

~~~~~~~~~~~Databases~~~~~~~~~~~
primerfile : $primerfile
bedfile    : $bedfile
blastdb    : $blastdbpath/blast
KMAdb      : $KMAdb
STARfasta  : $STARfasta

~~~~~~~~~~~Authors~~~~~~~~~~~~~~
        J.P.M. Coolen
            V. Karvink
            T. Baltussen
================================
"""

// Get UMI from cleaned reads
process '1A_UMI_reads' {
    tag '1A'
    conda 'bioconda::fastp=0.23.2 bioconda::pyfastx=0.6.12 conda-forge::simplejson=3.17.0'
    publishDir outDir + '/report/UMI', mode: 'copy'
    input:
        set pairID, file(reads) from reads_ch1
    output:
        set file("${samplename}_R1_fastp.fastq.gz"), file("${samplename}_R2_fastp.fastq.gz") into fastp_2A
        file "${samplename}.fastp.json"
        file "${samplename}.fastp.html"
        file ".command.*"
    script:
        """
        fastp -i ${reads[0]} -o ${samplename}_R1_fastp.fastq.gz \
                                -I ${reads[1]} -O ${samplename}_R2_fastp.fastq.gz \
                                --umi_len=${UMILEN} --umi --umi_loc=per_read --umi_prefix=UMI \
                                --html ${samplename}.fastp.html --json ${samplename}.fastp.json \
                                --length_required ${minreadlength} --trim_poly_g --trim_poly_x
        """
}

// Clean reads (adapter and read length filter)
process '1B_clean_reads' {
    tag '1B'
    conda 'bioconda::fastp=0.23.2 bioconda::pyfastx=0.6.12 conda-forge::simplejson=3.17.0'
    publishDir outDir + '/fastp', mode: 'copy'
    input:
        set pairID, file(reads) from reads_ch2
    output:
        set file("${samplename}_R1_fastp.fastq.gz"), file("${samplename}_R2_fastp.fastq.gz") into fastp_1C
        file "${samplename}.fastp.json"
        file "${samplename}.fastp.html"
        file ".command.*"
    script:
        """
        fastp -i ${reads[0]} -o ${samplename}_R1_fastp.fastq.gz \
                                -I ${reads[1]} -O ${samplename}_R2_fastp.fastq.gz \
                                --html ${samplename}.fastp.html --json ${samplename}.fastp.json \
                                --length_required ${minreadlength} --trim_poly_g --trim_poly_x
        """
}

// Clean reads (adapter and read length filter)
process '1C_remove_primers' {
    tag '1C'
    conda 'bioconda::ptrimmer=1.3.3'
    publishDir outDir + '/fastp', mode: 'copy'
    input:
        file reads from fastp_1C
    output:
        set file("${samplename}_R1_fastp.trim.fastq.gz"), file("${samplename}_R2_fastp.trim.fastq.gz") into fastp_3A, fastp_3C, fastp_6
        file "*"
    script:
        """
        ptrimmer -t pair --ampfile ${ptrimfile} --read1 ${reads[0]} --trim1 ${samplename}_R1_fastp.trim.fastq \
        --read2 ${reads[1]} --trim2 ${samplename}_R2_fastp.trim.fastq --summary ${samplename}.summary.ampcount.txt

        gzip ${samplename}_R1_fastp.trim.fastq
        gzip ${samplename}_R2_fastp.trim.fastq
        """
}

// Measure the UMI to get number of amplicons
process '2A_measure_amplicons' {
    tag '2A'
    conda 'conda-forge::pandas=1.2.4 bioconda::pysam=0.19.0 bioconda::htslib=1.15.1 anaconda::openpyxl'
    publishDir outDir + '/report/UMI', mode: 'copy'
    input:
        set file(forward_read), file(reverse_read) from fastp_2A
    output:
        file "*.csv" into UMI_2B
        file ".command.*"
    script:
        """
        python $baseDir/scripts/main.py --input $forward_read --primers ${primerfile}
        """
}

// Clean reads (adapter and read length filter)
process '2B_plot_amplicons' {
    tag '2B'
    conda "${baseDir}/conda/env-nodejs"
    publishDir outDir + '/report/UMI', mode: 'copy'
    input:
        file UMI from UMI_2B
    output:
        file ".command.*"
    script:
        """
        python ${baseDir}/conda/env-nodejs/sunburst/data-processing/sunburst-parsing.py -i $UMI \
        -o ${baseDir}/conda/env-nodejs/sunburst/Data.csv

        cd ${baseDir}/conda/env-nodejs/sunburst/

        #excute nodejs code
        npm run build

        # copy files manual to report folder
        mkdir -p ${outDir}/report/UMI
        cp index.html ${outDir}/report/UMI
        cp -r dist ${outDir}/report/UMI/
        cp style.css ${outDir}/report/UMI
        cp ${baseDir}/conda/env-nodejs/sunburst/Data.csv ${outDir}/report/UMI

        """
}


// create KMA tool to detect 16S
// Process 3A: KMA
process '3A_KMA' {
  tag '3A'
  time "30m"
  conda 'bioconda::kma=1.3.28'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file reads from fastp_3A
  output:
    file "${samplename}*"
    file "${samplename}.fsa" into consensus_4C
    file "${samplename}.sam" into kma_3B
    file "${samplename}.res" into classify_7B
    file ".command.*"
  script:
    """
    kma -t_db ${KMAdb} -ipe ${reads[0]} ${reads[1]} -t ${threads} -1t1 -a -ex_mode -ef -and -apm f -o ${samplename} -sam 4 > ${samplename}.sam 2>/dev/null || exit 0
    #kma -t_db ${KMAdb} -ipe ${reads[0]} ${reads[1]} -t ${threads} -a -ex_mode -ef -dense -1t1 -ref_fsa -mem_mode -and -apm f -o kma 2>/dev/null || exit 0
    """
}

// create KMA tool to match primers
// Process 3B: KMA
process '3B_process_KMA' {
  tag '3B'
  conda 'bioconda::samtools=1.12'
  errorStrategy 'ignore'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file sam from kma_3B
  output:
    file "${samplename}.sorted.bam"
    file "${samplename}.sorted.bam.bai"
    file ".command.*"
  script:
    """
    # sam --> bam
    samtools view -b ${sam} > ${samplename}.bam
    # sort bam
    samtools sort ${samplename}.bam > ${samplename}.sorted.bam
    # index bam
    samtools index ${samplename}.sorted.bam

    # additionally filter reference on only hits
    # this would make it possible to quickly evaluate the results

    """
}

// create KMA tool to detect 16S
// Process 3A: KMA
process '3C_STAR' {
  tag '3C'
  time "30m"
  conda 'bioconda::star=2.7.9a bioconda::samtools=1.12'
  //conda 'bioconda::star=2.7.10a'
  publishDir outDir + '/star', mode: 'copy'
  input:
  file reads from fastp_3C
  output:
    file "*.out"
    file "${samplename}.final.bam" into bam_4A, bam_4B
    file "${samplename}.final.bam.bai" into bamindex_4A, bamindex_4B
    file "${samplename}.unique.sorted.bam" into bam_5B
    file "${samplename}.unique.sorted.bam.bai" into bamindex_5B
    file ".command.*"
  script:
    """
    ## create database
    #STAR --runMode genomeGenerate --genomeDir "${baseDir}/db/CYP51A/STAR/" --genomeFastaFiles $STARfasta --genomeSAindexNbases 4
    STAR --genomeDir "${baseDir}/db/${snp_database}/STAR/" --readFilesCommand zcat --readFilesIn ${reads[0]} ${reads[1]} \
    --scoreDelOpen 0 --scoreDelBase 0 --scoreInsOpen 0 --scoreInsBase 0 \
    --alignIntronMax 100 \
    --runThreadN ${threads} --outFileNamePrefix ${samplename} --limitBAMsortRAM 1001609349 --outSAMtype BAM SortedByCoordinate

    #samtools rmdup ${samplename}Aligned.sortedByCoord.out.bam ${samplename}.unique.bam
    # sort bam file unique
    #samtools sort ${samplename}.unique.bam > ${samplename}.unique.sorted.bam

    # rename bam file
    mv ${samplename}Aligned.sortedByCoord.out.bam ${samplename}.final.bam
    cp ${samplename}.final.bam ${samplename}.unique.sorted.bam

    # index bam files
    samtools index ${samplename}.unique.sorted.bam
    samtools index ${samplename}.final.bam
    """
}

// Process 4A: primerdepth
process '4A_primerdepth' {
    tag '4A'
    conda 'bioconda::mosdepth=0.3.1'
    publishDir outDir + '/mosdepth', mode: 'copy'
    input:
        file bam from bam_4A
        file bamindex from bamindex_4A
    output:
        file "*" into data_6
        file ".command.*"
    script:
        """
        mosdepth --fast-mode --no-per-base --threads $threads ${samplename} ${bam}
        """
}

// Process 4B: freebayes
process '4B_freebayes' {
  tag '4B'
  conda 'bioconda::freebayes=1.3.6'
  publishDir outDir + '/freebayes', mode: 'copy'
  input:
    file bam from bam_4B
    file bamindex from bamindex_4B
  output:
    file "${samplename}.vcf" into vcf_5A
    file ".command.*"
    file "*"
  script:
    """
    freebayes -f ${STARfasta} --ploidy 1 ${bam} > ${samplename}.vcf
    """
}

// create abricate to detect 16S
// Process 4C: abricate
process '4C_abricate' {
  tag '4C'
  conda 'bioconda::abricate=1.0.1'
  publishDir outDir + '/abricate', mode: 'copy'
  input:
    file consensus from consensus_4C
  output:
    file "${samplename}_blast.csv" into blast_7B
    file ".command.*"
  script:
    if(abricate==true)
    """
    abricate --datadir ${blastdbpath} --db blast ${consensus} --mincov 30 --minid 60 --threads ${threads}  > ${samplename}_blast.csv
    """
    else if(abricate==false)
    """
    echo 'none' > ${samplename}_blast.txt
    """
}

// 5A: annotation of the genome/consensus fasta
process '5A_annotation' {
    tag '5A'
    conda 'bioconda::snpeff=5.0 bioconda::bcftools=1.12'
    publishDir outDir + '/annotation', mode: 'copy'
    input:
        file vcf from vcf_5A
    output:
        file "${samplename}.final.vcf"
        file "${samplename}_annot_table.txt" into annotation_7
        file ".command.*"
  script:
        """
        bcftools view -f . ${vcf} > ${samplename}.pass.vcf
        bcftools reheader -f "${STARfasta}.fai" -o ${samplename}.pass.correct.vcf ${samplename}.pass.vcf
        snpEff CYP51A ${samplename}.pass.correct.vcf -hgvs1LetterAa > ${samplename}.final.vcf

        ${baseDir}/conda/env-variantcalling/bin/python $vcf2table ${samplename}.final.vcf --sample ${samplename} \
        -ad -e -o ${samplename}_annot_table.txt
        """
}

// 5B: create visualization
process '5B_visualization' {
    tag '5B'
    conda 'bioconda::jvarkit-bam2svg'
    publishDir outDir + '/report', mode: 'copy'
    input:
        file bam from bam_5B
        file bamindex from bamindex_5B
    output:
        file "${samplename}_CYP51A.svg" into SVG_7B
        file "file.list"
        file ".command.*"
  script:
        if(visuals==true)
        """
        echo "${bam}" > file.list
        bam2svg.sh -i "CYP51A:0-2000" --groupby readgroup -w 800 -R ${STARfasta} -o ${samplename}_CYP51A.svg file.list
        """
        else if(visuals!=true)
        """
        echo "none" > "${samplename}_CYP51A.svg"
        echo "none" > file.list
        """
}

// Process 6: multiQC
process '6_multiQC' {
  tag '6'
  conda 'bioconda::multiqc'
  publishDir outDir + '/QC', mode: 'copy'
  input:
  file reads from fastp_6
  file data from data_6
  output:
    file "${samplename}.mosdepth.summary.txt" into mosdepth_7B
    file "*.html"
    file ".command.*"
  script:
    """
    multiqc ${outDir}/fastp/ ${outDir}/mosdepth/
    """
}

// Process 7A: obtain run parameters
process '7A_parameters' {
    tag '7A'
    publishDir outDir + '/report', mode: 'copy'
    input:
    output:
        file "parameters.txt" into params_7B
    script:
        """
        touch parameters.txt
        echo "Parameter\tValue" >> parameters.txt
        echo "Reads\t$params.reads" >> parameters.txt
        echo "Database\t$params.database" >> parameters.txt
        echo "abricate:\t$params.abricate" >> parameters.txt
        echo "minreadlength:\t$params.minreadlength" >> parameters.txt
        echo "UMI length:\t$params.UMILEN" >> parameters.txt
        """
}

// Process 7B: generate a report for interpretation by the clinician (or for research purposes)
process '7B_report' {
    tag '7B'
    conda "${baseDir}/conda/env-025066a104bf8ce5621e328d8009733a"
    publishDir outDir + '/report', mode: 'copy'
    input:
        file annotation from annotation_7
        file params from params_7B
        file mosdepth from mosdepth_7B
        file classify from classify_7B
        file blast from blast_7B
        file CYP51SVG from SVG_7B
    output:
        file "${samplename}.html"
        file "${samplename}.pdf"
        file ".command.*"
    script:
        """
        $reporter --sampleName ${samplename} \
        --annotation ${annotation} --params ${params} --mosdepth ${mosdepth} \
        --classification ${classify} --blast ${blast} --visualization ${CYP51SVG}
        """
}