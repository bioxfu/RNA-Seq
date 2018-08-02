## run on HPC ##
nohup snakemake --cluster qsub -j 32 -rp --latency-wait 3600 >> nohup.log 2>&1 &
