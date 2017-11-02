library(yaml)

# set parameters
argv <- commandArgs(T)
config_file <- argv[1]
count_dir <- argv[2]
count_all <- argv[3]

# load data
config <- yaml.load_file(config_file)
files <- paste0(count_dir, '/', config$samples, '_cnt.tsv')
count <- do.call("cbind", lapply(files, read.table, sep='\t', row.names=1))
colnames(count) <- config$samples
count <- count[-grep('^__', rownames(count)), ]
write.table(count, count_all, col.names=NA, sep='\t', quote=F)
