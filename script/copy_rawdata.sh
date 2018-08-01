DIR=$1
DOWNLOAD=$2

if [ "$DIR" == "" ]; then
	echo "usage: $0 data_path [--download]"
	exit
fi

mkdir fastq

host=xfu@10.41.25.100

if [ "$DOWNLOAD" == "" ]; then
	rsync -var $DIR ./fastq/
elif [ "$DOWNLOAD" == "--download" ]; then
	rsync -var $host:$DIR ./fastq/
else
	echo 'the 2nd parameter is --download or empty'
fi
