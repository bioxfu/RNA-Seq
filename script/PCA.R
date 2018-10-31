library(RColorBrewer)
library(yaml)

# set parameters
argv <- commandArgs(T)
config_file <- argv[1]
cpm_all <- argv[2]
#config_file <- '../config.yaml'
#cpm_all <- '../table/expr_table_cpm_all.tsv'
pca_output <- argv[3]

# load data
config <- yaml.load_file(config_file)
cpm <- read.table(cpm_all, head=T, row.names = 1, sep='\t', quote = '', comment.char = '')
cpm <- cpm[1:(grep('_vs_', colnames(cpm))[1] - 1)]
cpm <- log2(cpm+1)
pca <- prcomp(t(cpm), center = T, scale. = T)

x <- pca$x
prop_var <- round(summary(pca)$importance[2,1:2]*100,0)

set2_cols <- brewer.pal(8, 'Set1')
cols <- rep(set2_cols, each=config$seq_info$replicate)

pdf(pca_output, hei=7, wid=7)
layout(matrix(c(1,2),nrow=1), wid=c(5, 2))
par(mar=c(5,4,4,0))
par(xpd=TRUE)
plot(x[,1], x[,2], xlim=range(x[,1])*1.1, ylim=range(x[,2])*1.1, col=cols, cex=1, type='n',
     xlab=paste0('PC1 (',prop_var[1],'% of Variance)'),
     ylab=paste0('PC2 (',prop_var[2],'% of Variance)')
)
text(x[,1], x[,2], colnames(cpm), col=cols)
par(mar=c(5,0,4,0))
plot.new()
legend('left', unique(config$groups), pch = 1, col=set2_cols, bty='n')
dev.off()

