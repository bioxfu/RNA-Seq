library(gplots)
library(eulerr)
library(RColorBrewer)
library(yaml)

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

config_file <- 'config.yaml'
config <- yaml.load_file(config_file)
VS <- read.table(config$venn_table, sep='\t', header = T, stringsAsFactors = F)

dfm <- read.table(input, header = T, sep = '\t', quote = '', comment.char = '')

VS_list <- split(VS, VS$Venn)

pdf(output, wid=7, hei=7)
for (i in 1:length(VS_list)) {
  vs <- VS_list[[i]]
  glist <- list()
  for (n in 1:nrow(vs)) {
    glist[[n]] <- dfm$Gene[dfm[, vs[n, 'VS']]==vs[n, 'Index']]
  }
  names(glist) <- vs$Group
  vn <- venn(glist, show.plot = F)
  x <- sapply(attr(vn, 'intersections'), length)
  names(x) <- sub(':', '&', names(x))
  
  col_set=brewer.pal(4, 'Set3')
  print(plot(euler(x), quantities=T, legend=T, fills=list(fill=col_set, alpha=1)))
}
dev.off()

