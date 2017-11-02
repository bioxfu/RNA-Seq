library(yaml)

argv <- commandArgs(T)
config_file <- argv[1]
cpm_all <- argv[2]
cpm_DEG <- argv[3]
cpm_all_anno <- argv[4]
cpm_DEG_anno <- argv[5]

config <- yaml.load_file(config_file)

gene_anno <- function(x) {
  dfm <- read.table(x, header = T)
  anno <- read.table(config$gene_anno, sep='\t', header = T)
  dfm_anno <- merge(dfm, anno, by.x = 1, by.y = 1, all.x = T)
  return(dfm_anno)
}

write.table(gene_anno(cpm_all), cpm_all_anno, quote = F, row.names = F, sep = '\t')
write.table(gene_anno(cpm_DEG), cpm_DEG_anno, quote = F, row.names = F, sep = '\t')
