DRY=$1

if [[ "$DRY" == "--dry-run" ]] || [[ "$DRY" == "" ]]; then
	ls ./fastq/*/*/*.gz |./script/rush -k 'mv {} fastq/{/%}_R{%@.+[_.]R*([12]).fa*s*t*q.gz}.fastq.gz' $DRY
else
	echo 'the parameter is --dry-run or empty'
fi
