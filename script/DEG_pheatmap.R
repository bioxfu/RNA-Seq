library(pheatmap)
library(RColorBrewer)

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

dfm <- read.table(input, header = T, sep = '\t', quote = '', comment.char = '', row.names = 1)
rpkm <- dfm[1:(grep('_vs_', colnames(dfm))[1]-1)]

# set color
col_set <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

pdf(output, wid=5, hei=8)
pheatmap(rpkm, scale='row', cluster_cols=T, show_rownames = F, border_color = NA, color = col_set)
dev.off()
