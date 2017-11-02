library(yaml)

argv <- commandArgs(T)
config_file <- argv[1]
input_file <- argv[2]
output_file <- argv[3]

config <- yaml.load_file(config_file)

gene_anno <- function(x) {
  dfm <- read.table(x, header = T, skip=9, comment.char = '')
  anno <- read.table(config$gene_anno, sep='\t', header = T)
  dfm_anno <- merge(dfm, anno, by.x = 1, by.y = 1, all.x = T)
  return(dfm_anno)
}

write.table(gene_anno(input_file), output_file, quote = F, row.names = F, sep = '\t')
