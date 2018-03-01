library(yaml)
# set parameters
argv <- commandArgs(T)
config_file <- argv[1]
count_file <- argv[2]
output <- argv[3]
#config_file <- '../config.yaml'
#count_file <- '../table/AS_summary.tsv'
print(count_file)
# load data
config <- yaml.load_file(config_file)
count <- read.table(count_file, sep='\t')

anno <- read.table(config$gene_anno, sep='\t', header = T, quote = '', row.names = 1)
count_table <- merge(count, anno, by.x = 5, by.y = 0, all.x = T)
colnames(count_table) <- c('Gene', 'Position', 'Substitution', 'Frequency', 'P-value', 'Sample', 'GeneModeType', 'Description')

write.table(count_table, output, row.names=F, sep='\t', quote=F)
