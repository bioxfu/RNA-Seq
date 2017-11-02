## edit the following three lines
SRRs=(SRR001 SRR002 SRR003 SRR004 SRR005 SRR006)
names=(WT_1 WT_2 WT_3 KO_1 KO_2 KO_3)
LAYOUT=SE

cd fastq

## download .sra files using sratools
for SRR in ${SRRs[@]}; do
	prefetch -v $SRR
done

## convert the sra files to fastq files
if [ $LAYOUT == SE ]; then
	find ~/ncbi/public/sra/*.sra|parallel --gnu "fastq-dump -B -I --split-spot --skip-technical --gzip {}"
elif [ $LAYOUT == PE ]; then
	find ~/ncbi/public/sra/*.sra|parallel --gnu "fastq-dump -B -I --split-spot --split-files --skip-technical --gzip {}"
fi

## make md5sum for converted fastq files
md5sum * > checksum

## rename the fastq files
if [ $LAYOUT == SE ]; then
	for idx in ${!SRRs[@]}; do
		mv ${SRRs[$idx]}.fastq.gz ${names[$idx]}.fastq.gz
	done
elif [ $LAYOUT == PE ]; then
	for idx in ${!SRRs[@]}; do
		mv ${SRRs[$idx]}_R1.fastq.gz ${names[$idx]}_R1.fastq.gz
		mv ${SRRs[$idx]}_R2.fastq.gz ${names[$idx]}_R2.fastq.gz
	done
fi
