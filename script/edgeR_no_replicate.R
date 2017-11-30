library(edgeR)
library(RColorBrewer)
library(yaml)
library(plotrix)

# set parameters
argv <- commandArgs(T)
config_file <- argv[1]
count_file <- argv[2]
#config_file <- '../config.yaml'
#count_file <- '../count/all_sample_cnt.tsv'
  
# load data
config <- yaml.load_file(config_file)
count <- read.table(count_file, header = 1, row.names = 1)

group <- factor(config$groups)

## reading the data into the DGEList object ##
print(group)
y <- DGEList(counts = count, group = group)

## filtering lowly expressed genes ##
minGroupSize <- min(table(group))
keep <- rowSums(cpm(y) > 1) >= minGroupSize
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

## Normalizing ##
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=T)
#y <- estimateTagwiseDisp(y)

expr_cols <- round(cpm(y), 2)
expr_table <- as.data.frame(expr_cols)

anno <- read.table(config$gene_anno, sep='\t', header = T)
expr_table <- merge(expr_table, anno, by.x = 0, by.y = 1, all.x = T)
colnames(expr_table)[1] <- 'Gene'

write.table(expr_table, 'table/expr_table_cpm_no_replicate.tsv', row.names=F, sep='\t', quote=F)

