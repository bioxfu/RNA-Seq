library(rvest)

parse_html <- function(html) {
  html <- read_html(html)
  x <- html %>% html_nodes('.table-summary') %>% html_text()
  strsplit(x[3], '\n')[[1]][c(6,8)]
}

dfm <- data.frame(sample_name = sub('.bamqc', '', dir('bam', pattern = '*.bamqc')),
                  file = paste0(dir('bam', pattern = '*.bamqc', full=TRUE), '/qualimapReport.html'))


stat <- apply(dfm, 1, function(x){parse_html(x[2])}) %>% t()
rownames(stat) <- dfm$sample_name
colnames(stat) <- c('Number of reads', 'Mapped reads')

write.table(stat, 'stat/bamqc_stat.tsv', quote=F, col.names = NA, sep='\t')
