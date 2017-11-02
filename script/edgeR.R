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

# Differentially Expressed Genes (DEGs)
# fc: fold-change of expression level
fc <- 2
p_cutoff <- 0.05
fc_cutoff <- log2(fc)
group <- factor(config$groups)
VS <- as.data.frame(t(combn(unique(config$groups), 2)), stringsAsFactors = F)
colnames(VS) <- c('Ctrl', 'Expt')

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
y <- estimateTagwiseDisp(y)

## exact test ##
expr_cols <- round(cpm(y), 2)
regulate_cols <- NULL
logfc_cols <- NULL
fdr_cols <- NULL

for(i in 1:(nrow(VS))){
  ctr <- VS$Ctrl[i]
  exp <- VS$Expt[i]
  cat(paste('testing Exp:', exp, '- Ctr:', ctr, '\t'))
  et <- exactTest(y, pair=c(ctr, exp)) 
  # Note that the first group listed in the pair is the baseline for the comparison—
  # so if the pair is c("A","B") then the comparison is B - A, so genes with 
  # positive log-fold change are up-regulated in group B compared with group A 
  # (and vice versa for genes with negative log-fold change).
  de <- decideTestsDGE(et, p.value=p_cutoff, lfc=fc_cutoff, adjust.method='fdr')
  p <- sprintf('%.3e', p.adjust(et$table$PValue, method='fdr'))
  regulate_cols <- cbind(regulate_cols, de)
  logfc_cols <- cbind(logfc_cols, round(et$table$logFC, 3))
  fdr_cols <- cbind(fdr_cols, p)
  colnames(regulate_cols)[i] <- paste0(exp, '_vs_', ctr)      
  colnames(logfc_cols)[i] <- paste0(exp, '_vs_', ctr, '_logFC')      
  colnames(fdr_cols)[i] <- paste0(exp, '_vs_', ctr, '_FDR')      
  cat(paste('UP-regulated in', exp, sum(de > 0), 'genes\t'))
  cat(paste('DOWN-regulated in', exp, sum(de < 0), 'genes\n'))
}

expr_table <- as.data.frame(cbind(expr_cols, regulate_cols, logfc_cols, fdr_cols), stringsAsFactors = F)

anno <- read.table(config$gene_anno, sep='\t', header = T)
expr_table <- merge(expr_table, anno, by.x = 0, by.y = 1, all.x = T)
colnames(expr_table)[1] <- 'Gene'

DEG_table <- expr_table[rowSums(regulate_cols != 0) >= 1,]

write.table(expr_table, 'table/expr_table_cpm_all.tsv', row.names=F, sep='\t', quote=F)
write.table(DEG_table, 'table/expr_table_cpm_DEG.tsv', row.names=F, sep='\t', quote=F)

save(list = c('expr_table', 'DEG_table'), file = 'RData/edgeR_output.RData')

## Number of DEGs (barplot)
regulate.stat <- apply(as.matrix(regulate_cols), 2, function(x){table(x)[c('1','-1')]})
rownames(regulate.stat) <- c('up','down')
cols <- brewer.pal(3,'Set1')
upmax <- 1.2*max(regulate.stat['up',], na.rm=T)
dnmax <- 1.2*max(regulate.stat['down',], na.rm=T)

pdf('figure/DEG_barplot.pdf')
par(mar=c(8,4,4,2))
barplot(regulate.stat[1,], ylim=c(-dnmax,upmax), col=cols[2], border=cols[2], yaxt='n', las=2)
bp <- barplot(-regulate.stat[2,], add=T, names=NA, col=cols[3], border=cols[3], ylab='The number of DEGs', yaxt='n')
axis(2,at=pretty(0:upmax), label=pretty(0:upmax))
axis(2,at=-pretty(0:dnmax), label=pretty(0:dnmax))
text(bp,regulate.stat[1,], regulate.stat[1,], pos=3)
text(bp,-regulate.stat[2,], regulate.stat[2,], pos=1)
legend('topright', c('up','down'), fill=cols[2:3], bty='n',border=F)
dev.off()

## Number of DEGs (matrix)
col_reds <- brewer.pal(9, 'Reds')
col_blues <- brewer.pal(9, 'Blues')
mat_up <- NULL
mat_dn <- NULL
for (i in 1:nrow(VS)) {
  mat_up <- c(mat_up, rep(as.character(VS[i, ]), regulate.stat['up', i]))
  mat_dn <- c(mat_dn, rep(as.character(VS[i, ]), regulate.stat['down', i]))
}
mat_up <- matrix(mat_up, ncol=2, byrow=T)
mat_dn <- matrix(mat_dn, ncol=2, byrow=T)
diff_table_up <- table(as.data.frame(mat_up))
diff_table_dn <- table(as.data.frame(mat_dn))
gl <- levels(group)
n <- length(gl)
diff_table_all <- matrix(0, nrow=n, ncol=length(gl), dimnames = list(gl, gl))
diff_table_all[1:(n-1), 2:n] <- diff_table_all[1:(n-1), 2:n] + diff_table_up
diff_table_all[2:n, 1:(n-1)] <- diff_table_all[2:n, 1:(n-1)] + -t(diff_table_dn)
x <- diff_table_all
cellcol<-matrix(rep("#FFFFFF", n^2), nrow=n)
cellcol[x > 0] <- color.scale(x[x > 0], extremes=c(col_reds[1], col_reds[7]))
cellcol[x < 0] <- color.scale(x[x < 0], extremes=c(col_blues[7], col_blues[1]))

pdf('figure/DEG_matrix.pdf')
par(xpd=T)
par(mar=c(5,4,4,9))
color2D.matplot(x, cellcolors=cellcol, border='gray', axes = F, xlab='', ylab='', main='Number of DEGs')
mtext(side = 1, at=1:n-0.5, gl, line=1)
mtext(side = 2, at=1:n-0.5, rev(gl), line=1, las=2)
for (i in 1:n) {
  for (j in 1:n) {
    text(i-0.5, j-0.5, abs(x[n+1-j, i]))
  }
}
legend(n,n, c('up-regulated', 'down-regulated'), fill=c(col_reds[7], col_blues[7]), border = F, bty='n')
dev.off()
