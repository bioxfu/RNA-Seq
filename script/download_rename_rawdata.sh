DIR=$1
mkdir fastq

host=xfu@10.41.25.100
rsync -var $host:$DIR ./fastq/

ls ./fastq/*/*/* |./script/rush -k 'mv {} fastq/s{/%}_{%@.+[_.](R\d).fa*s*t*q.gz}.fastq.gz' --dry-run

ls ./fastq/*/*/* |./script/rush -k 'mv {} fastq/s{/%}_{%@.+[_.](R\d).fa*s*t*q.gz}.fastq.gz' 
