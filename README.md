#RC-PCR pipeline V0.2 (BETA)

# install
1. obtain docker image
```bash
docker pull jonovox/nextflowcentos:latest
```

2. download this project

3. extract this project

4. navigate to project

5. Download prebuild conda.tar.gz
https://surfdrive.surf.nl/files/index.php/s/5q2feFVult4v81k
   
6. Extract conda.tar.gz in folder

## single sample docker
```bash
sh docker/run.sh RC jonovox/nextflowcentos:latest
```

## batch run docker
```bash
# USAGE
cd project
# bash run_batch_docker.sh <inputpath> <file_extension> <database> <threads> <image>
# bash run_batch_docker.sh ${1}             ${2}          ${3}        ${4}     ${5}
# Example:
bash run_batch_docker.sh /file/location/ _001.fastq.gz CYP51A 8 jonovox/nextflowcentos:latest
# <file_extension> most common _001.fastq.gz
```

## FLOW-DIAGRAM
![Alt text](flowchart.png?raw=true "Flowdiagram")

##conda environments

  * 1A_clean_reads\
  env-f07c78eef9e8319c7eb087d931e36003
  * 2A_measure_amplicons\
    env-cd4ea0676bf53b5d7e5c6c6c523f0013
  * 3A_KMA (version 1.3.28)\
    env-84e06c5335c0a958ed012db619fdfceb
  * 3B_process_KMA\
    env-4ca5b26b8a059c60e73996439311c22f
  * 4A_abricate\
    env-b415f051979c22cdef40a3cbee1f0aa3 
  * 5A_annotation\
    env-9f6b61e20675ae28786fdb538092d4db
  * 6_multiQC\
    env-3abca7a24ea4d6c708bf4c6cea6413d2

#output

```bash
.
├── QC
│   └── multiqc_report.html
├── abricate
│   └── test_blast.txt        #abricate/blast result
├── annotation
│   └── test.final.vcf        #annotated vcf file
├── fastp
│   └── test.fastp.json
├── kma
│   ├── test.aln
│   ├── test.frag.gz
│   ├── test.frag_raw.gz
│   ├── test.fsa
│   ├── test.mapstat
│   ├── test.res              #KMA result file
│   ├── test.sam
│   ├── test.sorted.bam       #bam for genomebrowser
│   ├── test.sorted.bam.bai
│   └── test.vcf.gz
└── test_UMI_counttable.xlsx  #primer count table
```

#database structure

```bash
db
├── databasename
│   ├── KMA
│   │   ├── databasename.comp.b
│   │   ├── databasename.length.b
│   │   ├── databasename.name
│   │   └── databasename.seq.b
│   ├── blast
│   │   ├── sequences.fasta
│   │   ├── sequences.fasta.fai
│   │   ├── sequences.nhr
│   │   ├── sequences.nin
│   │   └── sequences.nsq
│   └── primers
│       └── databasename_primers.fasta
```

# Extra

# blastdb
18S
```bash
makeblastdb -in /workflow/db/blast_db/18S/sequences.fasta -title 18S -dbtype nucl -out /workflow/db/blast_db/18S/sequences
```

CYP51A (Afu4g06890)
```bash
makeblastdb -in /workflow/db/blast_db/CYP51A/sequences.fasta -title CYP51A -dbtype nucl -out /workflow/db/blast_db/CYP51A/sequences
```

# CYP51A
18S
```bash
kma_index -i /workflow/db/KMA/18S.fa -o /workflow/db/KMA/18S
```
CYP51A (Afu4g06890)
```bash
kma_index -i /workflow/db/KMA/CYP51A.fa -o /workflow/db/KMA/CYP51A
```

# snpEff
manual CYP51A (Afu4g06890)
```bash
snpEff build -gff3 CYP51A
```

# SILVA database
SILVA_138.1_SSURef_NR99_tax_silva_trunc
convert rRNA to DNA
```bash
perl -pe 'tr/tU/uT/ unless(/>/)' < db/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta > SILVA_138.1_SSURef_NR99_tax_silva_trunc_DNA.fasta
```

# aggregation example of UMI counttables
copy UMI counttable files to new folder
```bash
cp <pathto>/*/report/UMI/*counttable.csv <pathto>/summary
```
create summary on files
```bash
python scripts/aggregation/aggregate_UMI.py -i <pathto>/summary -o <pathto>/summary
```

# NOTES
# not implemented
splitting of samples on UMI using seqkit
seqkit grep -irp UMI samplename.fastq.gz > output.fastq