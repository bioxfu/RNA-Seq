configfile: "config.yaml"

rule all:
	input:
		expand('clean/{sample}_clean.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_clean_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('bam/{sample}.bam', sample=config['samples']),
		expand('bam/{sample}.bam.bai', sample=config['samples']),
		expand('bam/{sample}.bamqc', sample=config['samples']),
		expand('bam/{sample}.infer_strand_spec', sample=config['samples']),
		expand('stat/bamqc_stat.tsv'),
		expand('track/{sample}.tdf', sample=config['samples']),
		expand('count/{sample}_cnt.tsv', sample=config['samples']),
		'count/all_sample_cnt.tsv',
		expand('table/CPM_table_FDR0.05_FC{fc}_all.tsv', fc=config['fold_change']),
		expand('table/CPM_table_FDR0.05_FC{fc}_DEG.tsv', fc=config['fold_change']),
		expand('table/RPKM_table_FDR0.05_FC{fc}_all.tsv', fc=config['fold_change']),
		expand('table/RPKM_table_FDR0.05_FC{fc}_DEG.tsv', fc=config['fold_change']),
		expand('figure/DEG_barplot_FDR0.05_FC{fc}.pdf', fc=config['fold_change']),
		expand('RData/edgeR_output_FDR0.05_FC{fc}.RData', fc=config['fold_change']),
		'figure/PCA.pdf',
		#'figure/DEG_volcano_and_pie.pdf',
		#'figure/DEG_venn.pdf',
		#'figure/DEG_pheatmap.pdf',
#		if no replicate is available:
#		'table/expr_table_cpm_no_replicate.tsv',

rule fastqc_raw_SE:
	input:
		config['path']+'/{sample}.fastq.bz2'
	output:
		'fastqc/raw/{sample}_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -o fastqc/raw {input}'

rule trimmomatic_SE:
	input:
		config['path']+'/{sample}.fastq.bz2'
	output:
		'clean/{sample}_clean.fastq.gz'
	params:
		conda = config['conda_path'],
		adapter = config['adapter']
	shell:
		'{params.conda}/trimmomatic SE -threads 3 -phred33 {input} {output} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule fastqc_clean_SE:
	input:
		'clean/{sample}_clean.fastq.gz'
	output:
		'fastqc/clean/{sample}_clean_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -o fastqc/clean {input}'

rule fastqc_stat_SE:
	input:
		['fastqc/raw/{sample}_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_clean_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		conda = config['conda_path'],
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule hisat2_SE:
	input:
		'clean/{sample}_clean.fastq.gz'
	output:
		bam = 'bam/{sample}.bam'
	params:
		conda = config['conda_path'],
		prefix = 'bam/{sample}',
		cpu = config['cpu'],
		index = config['index'],
		strandness_hisat2 = config['strandness_hisat2']
	shell:
		"{params.conda}/hisat2 --rna-strandness {params.strandness_hisat2} -p {params.cpu} --dta -x {params.index} -U {input} |samtools view -Shub|samtools sort - -T {params.prefix} -o {output.bam}"

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
		cpu = config['cpu'],
		conda = config['conda_path']
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
		bed = config['bed']
	shell:
		'/cluster/home/xfu/Python/rseqc/bin/infer_experiment.py -r {params} -i {input} > {output}'

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
		gtf = config['gtf'],
		strandness_htseq = config['strandness_htseq'],
		conda = config['conda_path']		
	shell:
		"{params.conda}/htseq-count --format=bam --order=name --stranded={params.strandness_htseq} {input.bam} {params.gtf} > {output.cnt}"

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
		IGV = config['IGV'],
		conda = config['conda_path']		
	shell:
		"{params.conda}/igvtools toTDF {input.bg} {output.tdf} {params.IGV}"

rule all_sample_cnt:
	input:
		['count/{sample}_cnt.tsv'.format(sample=x) for x in config['samples']]
	output:
		count_all = 'count/all_sample_cnt.tsv',
	params:
		config_file = 'config.yaml',
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/all_sample_cnt.R {params.config_file} count {output}'

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
