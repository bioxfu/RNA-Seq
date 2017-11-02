library(fastqcr)
library(magrittr)

qc_raw <- qc_aggregate('fastqc/raw') %>% qc_stats()
qc_clean <- qc_aggregate('fastqc/clean') %>% qc_stats()
qc <- cbind(qc_raw, qc_clean)
qc <- qc[, c(1,4,9)]
colnames(qc) <- c('sample', 'tot.seq', 'tot.clean.seq')
qc$percentage <- round(as.numeric(qc$tot.clean.seq) / as.numeric(qc$tot.seq) * 100, 1)
  
write.table(qc, 'stat/fastqc_stat.tsv', sep='\t', row.names = F, quote=F)
