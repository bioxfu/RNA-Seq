library(pheatmap)
library(RColorBrewer)
library(scales)
library(yaml)
options(stringsAsFactors = F)

# set parameters
argv <- commandArgs(T)
config_file <- argv[1]
cpm_DEG <- argv[2]
profiles_output <- argv[3]
kmeans_all_output <- argv[4]
kmeans_TF_output <- argv[5]

# set color
col_set <- c(brewer.pal(8, 'Dark2'), brewer.pal(7, 'BuPu')[7])

# load data
config <- yaml.load_file(config_file)
tf <- read.table(config$TF, row.names=1)
colnames(tf) <- 'Family'
cpm <- read.table(cpm_DEG, head=T, row.names = 1)
cpm <- cpm[-grep('_vs_', colnames(cpm))]

groups <- factor(config$groups)
gl <- levels(groups)

cpm_average <- matrix(nrow=nrow(cpm), ncol=length(gl), dimnames = list(rownames(cpm), gl))
for (i in 1:length(gl)) {
  cpm_average[, i] <- rowMeans(cpm[groups == gl[i]])
}

cpm_scale <- t(scale(t(log2(cpm_average + 1))))

k <- config$kmean
set.seed(123)
kc <- pheatmap(cpm, scale='row', cluster_cols=F, show_rownames = F, border_color = NA, kmeans_k = k, silent = TRUE)
col_kc <- kc$kmeans$cluster
for (i in 1:k) {
  col_kc[col_kc == i] <- col_set[i]
}

pdf(profiles_output, wid=3, hei=8)
par(mfrow=c(6,1))
par(mar=c(2,5,0,4))
for (i in 1:k) {
  m <- cpm_scale[kc$kmeans$cluster==i,]
  m_median <- apply(m, 2, median)
  y_range <- c(-2, 2)
  plot(m_median, type='b', col=col_set[i], ylim=y_range, lwd=4, yaxt='n' ,xaxt='n', bty='l', ylab='Z-score', cex.lab=1.2, cex.axis=1.2)
  axis(2, at=c(-2,0,2), cex.axis=1.2)
  axis(1, at=1:length(gl), lab=gl, cex.axis=1.1)
  m_tf <- m[intersect(rownames(m), rownames(tf)), ]
  for (j in 1:nrow(m)) {
    points(m[j, ], type='l', col=alpha(col_set[i], 0.1), lwd=0.5)
  }
  points(m_median, type='b', col='gray', ylim=y_range, lwd=2, yaxt='n' ,xaxt='n', bty='l', ylab='Z-score', cex.lab=1.2, cex.axis=1.2)
  mtext(text=paste(nrow(m), 'genes'), side=4, line = 1, cex=0.8)
  mtext(text=paste(nrow(m_tf), 'TFs'), side=4, line = 2.5, cex=0.8)
}
dev.off()

## write to table
gene2cluster <- data.frame(cluster=kc$kmeans$cluster)
cpm_tf <- cpm[intersect(rownames(cpm), rownames(tf)),]
cpm_tf_cluster <- merge(cpm_tf, gene2cluster, by.x=0, by.y=0)
cpm_tf_cluster_fam <- merge(cpm_tf_cluster, tf, by.x=1, by.y=0)
colnames(cpm_tf_cluster_fam)[1] <- 'gene_id'

write.table(gene2cluster, kmeans_all_output, sep='\t', col.names = NA, quote = F)
write.table(cpm_tf_cluster_fam, kmeans_TF_output, sep='\t', row.names = F, quote = F)
