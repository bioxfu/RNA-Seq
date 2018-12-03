library(pheatmap)
library(RColorBrewer)
source('script/topGO.R')
library(clusterProfiler)
library(pathview)

dfm <- read.table('table/RPKM_table_FDR0.05_FC2_DEG.tsv', header = T, quote = '', comment.char = '', sep = '\t', stringsAsFactors = F)

# KEGG
loc_up <- dfm$AthID[dfm$loc_TYLCV_vs_loc_EV == 1]
loc_dn <- dfm$AthID[dfm$loc_TYLCV_vs_loc_EV == -1]
sys_up <- dfm$AthID[dfm$sys_TYLCV_vs_sys_EV == 1]
sys_dn <- dfm$AthID[dfm$sys_TYLCV_vs_sys_EV == -1]
loc_up <- loc_up[!is.na(loc_up)]
loc_dn <- loc_dn[!is.na(loc_dn)]
sys_up <- sys_up[!is.na(sys_up)]
sys_dn <- sys_dn[!is.na(sys_dn)]
olp_up <- intersect(loc_up, sys_up)
olp_dn <- intersect(loc_dn, sys_dn)

loc_up_kegg <- as.data.frame(enrichKEGG(gene=loc_up, organism = 'ath', pvalueCutoff = 0.05))[c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
loc_dn_kegg <- as.data.frame(enrichKEGG(gene=loc_dn, organism = 'ath', pvalueCutoff = 0.05))[c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
sys_up_kegg <- as.data.frame(enrichKEGG(gene=sys_up, organism = 'ath', pvalueCutoff = 0.05))[c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
sys_dn_kegg <- as.data.frame(enrichKEGG(gene=sys_dn, organism = 'ath', pvalueCutoff = 0.05))[c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
olp_up_kegg <- as.data.frame(enrichKEGG(gene=olp_up, organism = 'ath', pvalueCutoff = 0.05))[c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
olp_dn_kegg <- as.data.frame(enrichKEGG(gene=olp_dn, organism = 'ath', pvalueCutoff = 0.05))[c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust')]

kegg_lst <- list(loc_up=loc_up_kegg, loc_dn=loc_dn_kegg,
                 sys_up=sys_up_kegg, sys_dn=sys_dn_kegg,
                 olp_up=olp_up_kegg, olp_dn=olp_dn_kegg)

write2xlsx(kegg_lst, 'table/KEGG_enrichment.xlsx')

all_mat <- combine_kegg_list(kegg_lst = kegg_lst, pvalue = 0.05)
pdf('figure/KEGG_enrich_heatmap.pdf', height = 7)
pheatmap(all_mat, cluster_cols=F, cluster_rows=F, border_color = 'gray', col=rev(colorRampPalette(brewer.pal(n=7, name="RdYlBu"))(100)), legend_breaks=c(2.1,3,4,5,6), legend_labels = c('2', '3', '4','5','>=6'))
dev.off()

## KEGG pathway view
kegg <- c('ath04016', 'ath04075', 'ath04626')
loc_fc <- unlist(tapply(dfm$loc_TYLCV_vs_loc_EV_logFC, dfm$AthID, mean, simplify = F))
sys_fc <- unlist(tapply(dfm$sys_TYLCV_vs_sys_EV_logFC, dfm$AthID, mean, simplify = F))

for (k in kegg) {
  pathview(gene.data = loc_fc, gene.idtype = 'kegg', pathway.id = k, species = 'ath', out.suffix='loc_pathway', kegg.native = T)
  #pathview(gene.data = sys_fc, gene.idtype = 'kegg', pathway.id = k, species = 'ath', out.suffix='sys_pathway', kegg.native = T)
}

