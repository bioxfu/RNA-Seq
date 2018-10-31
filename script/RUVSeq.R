library(RUVSeq)
library(RColorBrewer)
library(yaml)
colors <- brewer.pal(3, 'Set2')

config_file <- 'config.yaml'
count_file <- 'count/all_sample_cnt.tsv'

# load data
config <- yaml.load_file(config_file)
count <- read.table(count_file, header = T, row.names = 1)

if(config$design_table == 'none') {
  VS <- as.data.frame(t(combn(unique(config$groups), 2)), stringsAsFactors = F)
  colnames(VS) <- c('Ctrl', 'Expt')
}else {
  VS <- read.table(config$design_table, sep='\t', header = T, stringsAsFactors = F)
}
geneLength <- read.table(config$gene_length, row.names = 1, header = T)
anno <- read.table(config$gene_anno, sep='\t', header = T, quote = '', row.names = 1, comment.char = '')

for (i in 1:nrow(VS)) {
  
  idx <- c(which(config$groups %in% VS[i, 'Ctrl']), which(config$groups %in% VS[i, 'Expt']))
  dat <- count[idx]
  minGroupSize <- min(table(config$groups))
  filter <- apply(dat, 1, function(x){length(x[x > 5]) >= minGroupSize})
  filtered <- dat[filter, ]
  
  x <- factor(config$groups[idx], levels = c(VS[i, 'Ctrl'], VS[i, 'Expt']))
  set <- newSeqExpressionSet(as.matrix(filtered), 
                             phenoData = data.frame(x, row.names = colnames(filtered)))
  
  ## Empirical control genes
  design <- model.matrix(~x, data = pData(set))
  y <- DGEList(counts=counts(set), group=x)
  y <- calcNormFactors(y, method = 'upperquartile')
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n=nrow(set))$table
  empirical <- rownames(set)[which(!rownames(set) %in% rownames(top)[1:5000])]
  set2 <- RUVg(set, empirical, k=1)
  
  pdf(paste0('figure/remove_batch_effect_', VS[i, 'Expt'], '_vs_', VS[i, 'Ctrl'], '.pdf'), wid=12)
  par(mfrow=c(1, 2))
  plotPCA(set, col=colors[x], main='before removing batch effect', cex=1.5)
  plotPCA(set2, col=colors[x], main='after removing batch effect', cex=1.5)
  dev.off()
  
  ## Differential expression analysis
  design <- model.matrix(~x + W_1, data = pData(set2))
  y <- DGEList(counts=counts(set2), group=x, genes = geneLength[rownames(counts(set2)), , drop=F])
  y <- calcNormFactors(y, method = 'upperquartile')
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n=nrow(set))$table[c('logFC', 'FDR')]
  
  rpkm <- round(rpkm(y), 2)
  
  result <- merge(merge(rpkm, top, by.x = 0, by.y = 0), anno, by.x = 1, by.y = 0, all.x = T)
  colnames(result)[1] <- 'geneID'
  DEG <- result[result$FDR < 0.05, ]
  
  ## output
  write.table(result, paste0('table/remove_batch_effect_', VS[i, 'Expt'], '_vs_', VS[i, 'Ctrl'], '_all.tsv'), row.names=F, sep='\t', quote=F)
  write.table(DEG, paste0('table/remove_batch_effect_', VS[i, 'Expt'], '_vs_', VS[i, 'Ctrl'], '_FDR0.05.tsv'), row.names=F, sep='\t', quote=F)

}
