library(RColorBrewer)
library(pheatmap)
options(stringsAsFactors=F)

# set parameters
argv <- commandArgs(T)
cpm_DEG <- argv[1]
hclust_output <- argv[2]
kmean_withinss <- argv[3]
k4 <- argv[4]
k6 <- argv[5]
k8 <- argv[6]

# set colors
set1_cols <- brewer.pal(9, 'Set1')
heat_cols <- colorRampPalette(c(set1_cols[2],"white",set1_cols[1]))(50)

# load data
cpm <- read.table(cpm_DEG, head=T, row.names = 1)
cpm <- cpm[-grep('_vs_', colnames(cpm))]
cpm <- log2(cpm+1)

png(hclust_output, hei=800, wid=400)
pheatmap(cpm, scale='row', cluster_cols=F, color = heat_cols, show_rownames = F, border_color = NA)
dev.off()

# k-means
wss <- NULL
for (i in 2:12) {
  set.seed(123)
  wss[i] <- kmeans(cpm,centers = i)$tot.withinss
}
pdf(kmean_withinss)
plot(wss, type='b', xlab='Number of clusters', ylab='Total within-cluster sum of squares')
dev.off()

png(k4, hei=400, wid=800)
set.seed(123)
pheatmap(cpm, scale='row', cluster_cols=F, color = heat_cols, show_rownames = T, border_color = NA, kmeans_k=4)
dev.off()

png(k6, hei=400, wid=800)
set.seed(123)
pheatmap(cpm, scale='row', cluster_cols=F, color = heat_cols, show_rownames = T, border_color = NA, kmeans_k=6)
dev.off()

png(k8, hei=400, wid=800)
set.seed(123)
pheatmap(cpm, scale='row', cluster_cols=F, color = heat_cols, show_rownames = T, border_color = NA, kmeans_k=8)
dev.off()
