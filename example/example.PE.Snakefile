configfile: "config.yaml"

rule all:
	input:
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/raw/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R1_paired_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R2_paired_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('bam/{sample}.bam', sample=config['samples']),
		expand('bam/{sample}.bam.bai', sample=config['samples']),
		expand('bam/{sample}.bamqc', sample=config['samples']),
		expand('bam/{sample}.infer_strand_spec', sample=config['samples']),
		expand('stat/bamqc_stat.tsv'),
#		expand('track/{sample}.tdf', sample=config['samples']),
		expand('count/{sample}_cnt.tsv', sample=config['samples']),
#		expand('count/{sample}_junction_counts', sample=config['samples']),
		'count/all_sample_cnt.tsv',
#		expand('AS/{sample}', sample=config['samples']),
#		expand('AS/{sample}/SE.txt', sample=config['samples']),
#		expand('AS/{sample}/RI.txt', sample=config['samples']),
#		expand('AS/{sample}/MXE.txt', sample=config['samples']),
#		expand('AS/{sample}/AFE.txt', sample=config['samples']),
#		expand('AS/{sample}/ALE.txt', sample=config['samples']),
#		expand('AS/{sample}/A3SS.txt', sample=config['samples']),
#		expand('AS/{sample}/A5SS.txt', sample=config['samples']),
#		expand('AS/all_sample_{AStype}_fisher_test.tsv', AStype=config['AStype']),
		expand('table/CPM_table_FDR0.05_FC{fc}_all.tsv', fc=config['fold_change']),
		expand('table/CPM_table_FDR0.05_FC{fc}_DEG.tsv', fc=config['fold_change']),
		expand('table/RPKM_table_FDR0.05_FC{fc}_all.tsv', fc=config['fold_change']),
		expand('table/RPKM_table_FDR0.05_FC{fc}_DEG.tsv', fc=config['fold_change']),
		expand('figure/DEG_barplot_FDR0.05_FC{fc}.pdf', fc=config['fold_change']),
		expand('RData/edgeR_output_FDR0.05_FC{fc}.RData', fc=config['fold_change']),
		'figure/PCA.pdf',
		'figure/DEG_volcano_and_pie.pdf',
		'figure/DEG_venn.pdf',
		'figure/DEG_pheatmap.pdf',
#		if no replicate is available:
#		'table/expr_table_cpm_no_replicate.tsv',

rule fastqc_raw_PE:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/raw/{sample}_R1_fastqc.html',
		'fastqc/raw/{sample}_R2_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -t 2 -o fastqc/raw {input}'

rule trimmomatic_PE:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		conda = config['conda_path'],
		adapter = config['adapter']
	shell:
		'{params.conda}/trimmomatic PE -threads 3 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule fastqc_clean_PE:
	input:
		'clean/{sample}_R1_paired.fastq.gz',
		'clean/{sample}_R2_paired.fastq.gz'
	output:
		'fastqc/clean/{sample}_R1_paired_fastqc.html',
		'fastqc/clean/{sample}_R2_paired_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -t 2 -o fastqc/clean {input}'

rule fastqc_stat_PE:
	input:
		['fastqc/raw/{sample}_R1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/raw/{sample}_R2_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R1_paired_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R2_paired_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule hisat2_PE:
	input:
		r1 = 'clean/{sample}_R1_paired.fastq.gz',
		r2 = 'clean/{sample}_R2_paired.fastq.gz'
	output:
		bam = 'bam/{sample}.bam'
	params:
		conda = config['conda_path'],
		prefix = 'bam/{sample}',
		cpu = config['cpu'],
		index = config['index'],
		strandness_hisat2 = config['strandness_hisat2']
	shell:
		"{params.conda}/hisat2 --rna-strandness {params.strandness_hisat2} -p {params.cpu} --dta -x {params.index} -1 {input.r1} -2 {input.r2} |{params.conda}/samtools view -Shub|{params.conda}/samtools sort - -T {params.prefix} -o {output.bam}"

rule bam_idx:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bai = 'bam/{sample}.bam.bai'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/samtools index {input.bam} {output.bai}'

rule bam_qc:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bamqc_dir = 'bam/{sample}.bamqc',
		bamqc_html = 'bam/{sample}.bamqc/qualimapReport.html'
	params:
		conda = config['conda_path'],
		cpu = config['cpu']
	shell:
		"{params.conda}/qualimap bamqc --java-mem-size=10G -nt {params.cpu} -bam {input.bam} -outdir {output.bamqc_dir}"

rule bam_qc_stat:
	input:
		['bam/{sample}.bamqc/qualimapReport.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/bamqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']		
	shell:
		"{params.Rscript} script/mapping_stat_by_bamqc.R"

rule infer_strand_spec:
	input:
		'bam/{sample}.bam'
	output:
		'bam/{sample}.infer_strand_spec'
	params:
		bed = config['bed'],
		rseqc = config['rseqc_path']
	shell:
		'{params.rseqc}/infer_experiment.py -r {params.bed} -i {input} > {output}'

rule sort_bam_by_name:
	input:
		'bam/{sample}.bam'
	output:
		'count/{sample}_sort_name.bam'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/samtools sort -n {input} -o {output}'

rule htseq:
	input:
		bam = 'count/{sample}_sort_name.bam'
	output:
		cnt = 'count/{sample}_cnt.tsv'
	params:
		conda = config['conda_path'],
		gtf = config['gtf'],
		strandness_htseq = config['strandness_htseq']
	shell:
		"{params.conda}/htseq-count --format=bam --order=name --stranded={params.strandness_htseq} {input.bam} {params.gtf} > {output.cnt}"

rule bamtobed12:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bed12 = 'bam/{sample}.bed12'
	params:
		conda = config['conda_path']
	shell:
		"{params.conda}/bedtools bamtobed -tag NH -split -bed12 -i {input} > {output}"

rule bed12tobed6:
	input:
		bed12 = 'bam/{sample}.bed12'
	output:
		block_2_3_bed12 = 'bam/{sample}_block_2_3.bed12',
		block_2_3_bed6 = 'bam/{sample}_block_2_3.bed6'
	params:
		conda = config['conda_path']
	shell:
		"cat {input} |awk '{{if($5==1 && ($10==2 || $10==3))print}}' > {output.block_2_3_bed12}; {params.conda}/bedtools bed12tobed6 -i {output.block_2_3_bed12} -n > {output.block_2_3_bed6}"

rule cal_junc_counts:
	input:
		block_2_3_bed6 = 'bam/{sample}_block_2_3.bed6'
	output:
		junction_counts = 'count/{sample}_junction_counts'
	params:
		conda = config['conda_path']
	shell:
		"python script/cal_junc_counts.py {input} > {output}"

rule ASFinder:
	input:
		bam = 'bam/{sample}.bam',
		junction_counts = 'count/{sample}_junction_counts'
	output:
		outdir = 'AS/{sample}',
		SE = 'AS/{sample}/SE.txt',
		RI = 'AS/{sample}/RI.txt',
		MXE = 'AS/{sample}/MXE.txt',
		AFE = 'AS/{sample}/AFE.txt',
		ALE = 'AS/{sample}/ALE.txt',
		A3SS = 'AS/{sample}/A3SS.txt',
		A5SS = 'AS/{sample}/A5SS.txt'
	params:
		conda = config['conda_path'],
		bed = config['AS'],
		bed_RI = config['AS_RI']
	shell:
		"python script/ASFinder.py {input.bam} {input.junction_counts} {params.bed} {params.bed_RI} {output.outdir}"

rule bam2count:
	input:
		bam = 'bam/{sample}.bam'
	output:
		cnt = 'bam/{sample}.cnt'
	params:
		conda = config['conda_path']
	shell:
		"{params.conda}/samtools view -c -F 4 {input.bam} > {output.cnt}"
		
rule bam2bedgraph:
	input:
		bam = 'bam/{sample}.bam',
		cnt = 'bam/{sample}.cnt'
	output:
		bg = 'track/{sample}.bedgraph'
	params:
		conda = config['conda_path']
	shell:
		"X=`awk '{{print 1/$1*1000000}}' {input.cnt}`; "
		"{params.conda}/bedtools genomecov -ibam {input.bam} -bga -scale $X -split > {output.bg}"

rule bedgraph2tdf:
	input:
		bg = 'track/{sample}.bedgraph'
	output:
		tdf = 'track/{sample}.tdf'
	params:
		conda = config['conda_path'],
		IGV = config['IGV']
	shell:
		"{params.conda}/igvtools toTDF {input.bg} {output.tdf} {params.IGV}"

rule all_sample_cnt:
	input:
		['count/{sample}_cnt.tsv'.format(sample=x) for x in config['samples']]
	output:
		count_all = 'count/all_sample_cnt.tsv',
	params:
		conda = config['conda_path'],
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/all_sample_cnt.R {params.config_file} count {output}'

rule merge_AS_tables:
	input:
		['AS/{sample}/SE.txt'.format(sample=x) for x in config['samples']],
		['AS/{sample}/RI.txt'.format(sample=x) for x in config['samples']],
		['AS/{sample}/MXE.txt'.format(sample=x) for x in config['samples']],
		['AS/{sample}/AFE.txt'.format(sample=x) for x in config['samples']],
		['AS/{sample}/ALE.txt'.format(sample=x) for x in config['samples']],
		['AS/{sample}/A3SS.txt'.format(sample=x) for x in config['samples']],
		['AS/{sample}/A5SS.txt'.format(sample=x) for x in config['samples']]
	output:
		'AS/all_sample_{AStype}.tsv'
	params:
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/merge_AS_tables.R {params.config_file}'

rule fisher_test:
	input:
		rules.merge_AS_tables.output
	output:
		'AS/all_sample_{AStype}_fisher_test.tsv'
	params:
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/AS_fisher_test.R {params.config_file} {input} {output}'

ruleorder: fisher_test > merge_AS_tables 

rule edgeR:
	input:
		count_all = 'count/all_sample_cnt.tsv',
	output:
		['table/CPM_table_FDR0.05_FC{fc}_all.tsv'.format(fc=x) for x in config['fold_change']],
		['table/CPM_table_FDR0.05_FC{fc}_DEG.tsv'.format(fc=x) for x in config['fold_change']],
		['table/RPKM_table_FDR0.05_FC{fc}_all.tsv'.format(fc=x) for x in config['fold_change']],
		['table/RPKM_table_FDR0.05_FC{fc}_DEG.tsv'.format(fc=x) for x in config['fold_change']],
		['figure/DEG_barplot_FDR0.05_FC{fc}.pdf'.format(fc=x) for x in config['fold_change']],
		['RData/edgeR_output_FDR0.05_FC{fc}.RData'.format(fc=x) for x in config['fold_change']],
	params:
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/edgeR.R {params.config_file} {input}'

rule PCA:
	input:
		cpm_all = 'table/CPM_table_FDR0.05_FC1.5_all.tsv'
	output:
		pca_output = 'figure/PCA.pdf'
	params:
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/PCA.R {params.config_file} {input} {output}'

rule volcano:
	input:
		'table/RPKM_table_FDR0.05_FC2_all.tsv'
	output:
		'figure/DEG_volcano_and_pie.pdf'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/DEG_volcano_and_pie.R {input} {output}'

rule venn:
	input:
		'table/RPKM_table_FDR0.05_FC2_all.tsv'
	output:
		'figure/DEG_venn.pdf'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/DEG_venn.R {input} {output}'

rule pheatmap:
	input:
		'table/RPKM_table_FDR0.05_FC2_all.tsv'
	output:
		'figure/DEG_pheatmap.pdf'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/DEG_pheatmap.R {input} {output}'

rule edgeR_no_replicate:
	input:
		count_all = 'count/all_sample_cnt.tsv',
	output:
		cpm_all = 'table/expr_table_cpm_no_replicate.tsv',
	params:
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/edgeR_no_replicate.R {params.config_file} {input}'

