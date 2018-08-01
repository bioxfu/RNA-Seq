DIR=$1
mkdir fastq

rsync -var $DIR ./fastq/

ls ./fastq/*/*/* |./script/rush -k 'mv {} fastq/s{/%}_{%@.+[_.](R\d).fa*s*t*q.gz}.fastq.gz' --dry-run

ls ./fastq/*/*/* |./script/rush -k 'mv {} fastq/s{/%}_{%@.+[_.](R\d).fa*s*t*q.gz}.fastq.gz' 
