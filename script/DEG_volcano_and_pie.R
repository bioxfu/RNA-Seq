library(RColorBrewer)
col_set <- brewer.pal(n = 3, 'Set1')

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

dfm <- read.table(input, header = T, sep = '\t', quote = '', comment.char = '')

pdf(output)
vs <- grep('_vs_', colnames(dfm), value = T)
for (i in 1:(length(vs)/3)) {
  print(vs[i])
  idx <- dfm[, vs[i]]
  logfc <- dfm[, paste0(vs[i], '_logFC')]
  fdr <- dfm[, paste0(vs[i], '_FDR')]
  gene_num <- c(dn=sum(idx == -1), up=sum(idx == 1))

  pt_col <- idx
  pt_col[pt_col == 1] <- col_set[1]  # up: red
  pt_col[pt_col == -1] <- col_set[2] # down: blue
  pt_col[pt_col == 0] <- 'gray'

  layout(matrix(c(1,1,2,3),nrow=2), wid=c(5, 3))
  par(mar=c(5,5,4,0))
  plot(logfc, -log10(fdr), col = pt_col, pch=16, ylab=expression('-log'[10]*'(FDR)'), xlab=expression('log'[2]*'(fold change)'), cex.lab=1.5, cex.axis=1.5, main=vs[i])
  abline(h=-log10(0.05), lty=2, lwd=1.5)
  par(mar=c(5,0,4,0))
  plot.new()
  legend('left', c('Up-regulated', 'Down-regulated'), fill=col_set[1:2], bty='n', border = NA, cex=1.5)
  pie(gene_num, col = col_set[2:1], border = 'white', labels = gene_num, cex=1.2, init.angle = 90)
}
dev.off()
