DIR=$1
mkdir fastq
scp -r xfu@10.41.25.100:$DIR ./fastq
ls ./fastq/*/*/* |./script/rush -k 'mv {} fastq/{/%}_{%@.+_(R\d).fq.gz}.fastq.gz' --dry-run
ls ./fastq/*/*/* |./script/rush -k 'mv {} fastq/{/%}_{%@.+_(R\d).fq.gz}.fastq.gz'
