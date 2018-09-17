# https://www.biostars.org/p/83901/
# First, import the GTF-file that you have also used as input for htseq-count
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/home/xfu/Gmatic7/gene/tair10/tair10.gtf",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb, by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sapply(exons.list.per.gene, function(x){sum(width(reduce(x)))})

dfm <- data.frame(Length = exonic.gene.sizes)

write.table(dfm, 'tair10_merged_exons_length.tsv', col.names = NA, quote = F, sep = '\t')
