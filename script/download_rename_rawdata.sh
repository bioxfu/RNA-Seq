###cd /var/data/raw_data/ZHM/CSL/
scp -r lpeng@10.41.25.251:$1 ./
ls */*/* |./rush 'mv {} {/%}_R{%@.+_(\d).fq.gz}.fastq.gz' --dry-run
ls */*/* |./rush 'mv {} {/%}_R{%@.+_(\d).fq.gz}.fastq.gz'
mkdir fastq
mv *.fastq.gz fastq

