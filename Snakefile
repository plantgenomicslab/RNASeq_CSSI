import sys
import os
import pandas as pd
SAMPLE_FILE = pd.read_table('samples.tsv', sep="\s+", dtype=str).set_index("sample", drop=False)  # enforce str in index
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()
THREADS = 16
configfile: "./config.json"

os.makedirs("output/cluster/logs/", exist_ok=True)
for path in SAMPLE_LIST:
	os.makedirs("output/" + path + "/stat/", exist_ok=True)
	os.makedirs("output/" + path + "/RAW/", exist_ok=True)




###################
##CHECK reference file, if there's no reference, please check fasta file and gff file. If the fasta and GFF is there, need to create with below commands
#chrLength_test = os.path.exists('chrLength.txt')
#chrNameLength_test = os.path.exists('chrNameLength.txt')
#chrName_test = os.path.exists('chrName.txt')
#exonGTI_test = os.path.exists('exonGeTrInfo.tab')
#exonInfo_test = os.path.exists('exonInfo.tab')
REF_FILES = ['chrNameLength.txt', 'chrName.txt', 'exonGeTrInfo.tab', 'exonInfo.tab', 'geneInfo.tab', 
					'Genome', 'genomeParameters.txt', 'SA', 'SAindex', 'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab', 
					'sjdbList.out.tab', 'transcriptInfo.tab']
missing_files = 0
os.makedirs('ref', exist_ok=True)
#for ref in REF_FILES:
#	if not os.path.exists('ref/' + ref):
#		missing_files = 1
#if missing_files:
#	os.system('curl -L https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff > ref/TAIR10_GFF3_genes.gff; curl -L https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > ref/TAIR10_chr_all.fas')

#https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
#https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

##REF = 'reference/'
##We will use this at mapping step

#os.system("gffread ref/TAIR10_GFF3_genes.gff -T -o ref/TAIR10.gtf")
#FASTA_LINES = open('ref/TAIR10_chr_all.fas', 'r').readlines()
#FASTA_OUTPUT = open('ref/TAIR10_chr_all.fas.fixed', 'w')
#chr_count = 0
#for row in FASTA_LINES:
#	if row[0] == '>':
#		FASTA_OUTPUT.write('>Chr' + row[1].upper() + '\n')
#	else:
#		FASTA_OUTPUT.write(row)
#FASTA_OUTPUT.close()

##need to create refrence folder
#os.system("STAR --runThreadN 6 --genomeDir ./reference/ --genomeFastaFiles ref/TAIR10_chr_all.fas.fixed --sjdbGTFfile ref/TAIR10.gtf --runMode genomeGenerate --sjdbOverhang 99 --genomeSAindexNbases 10")


###################


#os.makedirs("output/cluster/logs/", exist_ok=True)
#for path in SAMPLE_LIST:
#        os.makedirs("output/" + path + "/stat/", exist_ok=True)
#        os.makedirs("output/" + path + "/RAW/", exist_ok=True)


rule all:
	input:
		"output/Final_analysis/diffExpr.P1e-3_C1.matrix",
		"output/featureCount.cnt_for_tpm.tpm.tab",
		expand("output/{rundeg_input}", rundeg_input = config["rundeg-output"]),
		#"output/featureCount.cnt.fixed",
		expand("output/{sample}/{sample}Aligned.sortedByCoord.out.bam", sample = SAMPLE_LIST),
		#expand("output/{sample}/{rundeg_input}", sample = SAMPLE_LIST, rundeg_input = config["rundeg-output"])
		expand("output/{sample}/TRIM/{sample}_pass_2_val_2.fq.gz", sample = SAMPLE_LIST),
		expand("output/{sample}/TRIM/{sample}_pass_1_val_1.fq.gz", sample = SAMPLE_LIST)
################
### please implement
rule fastqdump:
	message: "-~- Downloading fastq files... -~-"
	output:
		fwd="output/{sample}/RAW/{sample}_pass_1.fastq.gz",
		rev="output/{sample}/RAW/{sample}_pass_2.fastq.gz"
	shell:
		"fastq-dump --outdir output/{wildcards.sample}/RAW/  --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.sample}"



################

rule trimming:
	message: "-~- Trimming fastq files... -~-"
	input:
		fwd="output/{sample}/RAW/{sample}_pass_1.fastq.gz",
		rev="output/{sample}/RAW/{sample}_pass_2.fastq.gz"
	output:
		fwd="output/{sample}/TRIM/{sample}_pass_1_val_1.fq.gz",
		rev="output/{sample}/TRIM/{sample}_pass_2_val_2.fq.gz"
	threads: THREADS
	run:
		shell('mkdir -p output/{wildcards.sample}/TRIM')
		shell('trim_galore \
				--paired \
				--three_prime_clip_R1 10 \
				--three_prime_clip_R2 10 \
				--cores {threads} \
				--max_n 50 \
				--gzip \
				-o output/{wildcards.sample}/TRIM/ \
				{input.fwd} \
				{input.rev}')
				#module = ['LSC', 'IRA', 'IRB', 'SSC']

rule mapping:
	threads: THREADS
	input:
		fwd="output/{sample}/TRIM/{sample}_pass_1_val_1.fq.gz",
		rev="output/{sample}/TRIM/{sample}_pass_2_val_2.fq.gz"
	output:
		file = "output/{sample}/{sample}Aligned.sortedByCoord.out.bam"
	message: "STAR mapping"
	run:
		#shell('mkdir -p output/{wildcards.sample}/')
		shell('STAR --runMode alignReads --runThreadN 10 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 30 --alignIntronMax 55000 --genomeDir ./reference  --readFilesIn {input.fwd},{input.rev}  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.sample}/{wildcards.sample}')

##################
rule readscount:
	input:
		bam = expand("output/{sample}/{sample}Aligned.sortedByCoord.out.bam", sample = SAMPLE_LIST)
	output:
		ftcount = "output/featureCount.cnt.fixed"
	shell:
		"featureCounts -p -Q 10 -M --fraction -s 0 -T 16 -o output/featureCount.cnt  -a ./ref/TAIR10.gtf {input.bam}; python src/fix_featCnt_header.py output/featureCount.cnt"

 ### EXAMPLE
 #Cont_1Aligned.sortedByCoord.out.bam Cont_2Aligned.sortedByCoord.out.bam Cont_3Aligned.sortedByCoord.out.bam Hours2_1Aligned.sortedByCoord.out.bam Hours2_2Aligned.sortedByCoord.out.bam Hours2_3Aligned.sortedByCoord.out.bam Hours4_1Aligned.sortedByCoord.out.bam Hours4_2Aligned.sortedByCoord.out.bam Hours4_3Aligned.sortedByCoord.out.bam Hours6_1Aligned.sortedByCoord.out.bam Hours6_2Aligned.sortedByCoord.out.bam Hours6_3Aligned.sortedByCoord.out.bam Hours8_1Aligned.sortedByCoord.out.bam Hours8_3Aligned.sortedByCoord.out.bam

###############


###############
rule rundeg:
	input:
		ftcount = "output/featureCount.cnt.fixed"
	output:
		rd_out = directory(expand("output/{rundeg_input}", rundeg_input = config["rundeg-output"]))
	shell:
		"run_DE_analysis.pl -m output/featureCount.cnt.fixed --method DESeq2 --samples_file ./samplefile --output {output.rd_out} --contrasts ./contrasts"
##############


###############
rule TPM:
	input:
		ftcount = "output/featureCount.cnt.fixed"
	output:
		ft_tab = "output/featureCount.cnt_for_tpm.tpm.tab"
	shell:
		"cut -f 1,6- {input.ftcount} | egrep -v '#' >> output/featureCount.cnt_for_tpm; python tpm_raw_exp_calculator.py -count output/featureCount.cnt_for_tpm"


##############

##############
rule diffexpr:
## Need to go to output folder
	input:
		tmp_tab = "output/featureCount.cnt_for_tpm.tpm.tab"
	output:
		analyze = "output/Final_analysis/diffExpr.P1e-3_C1.matrix"
	params:
		rundeg_file = expand("output/{rundeg_input}", rundeg_input = config["rundeg-output"])
	shell:
		"cd {params.rundeg_file}; analyze_diff_expr.pl --matrix ../featureCount.cnt_for_tpm.tpm.tab -P 1e-3 -C 1 --max_genes_clust 40000; mkdir -p ../Final_analysis/; mv diffExpr.P1e-3_C1* ../Final_analysis/; mv *.subset ../Final_analysis; mv *.samples ../Final_analysis; mv *.matrix ../Final_analysis; cd ../.."

#############
