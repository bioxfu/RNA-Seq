## run on HPC ##
nohup snakemake --cluster qsub -j 32 -rp >> nohup.log 2>&1 &
