source activate gmatic
#conda env export > doc/environment.yml

if [ ! -d fastqc ]; then
	mkdir -p fastq fastqc/raw fastqc/clean clean bam count track table stat figure RData AS
fi
