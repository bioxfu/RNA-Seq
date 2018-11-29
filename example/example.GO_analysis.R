library(pheatmap)
library(RColorBrewer)
source('script/topGO.R')

dfm <- read.table('table/RPKM_table_FDR0.05_FC2_DEG.tsv', header = T, quote = '', comment.char = '', sep = '\t', stringsAsFactors = F)

loc_up <- dfm$Gene[dfm$loc_TYLCV_vs_loc_EV == 1]
loc_dn <- dfm$Gene[dfm$loc_TYLCV_vs_loc_EV == -1]
sys_up <- dfm$Gene[dfm$sys_TYLCV_vs_sys_EV == 1]
sys_dn <- dfm$Gene[dfm$sys_TYLCV_vs_sys_EV == -1]
olp_up <- intersect(loc_up, sys_up)
olp_dn <- intersect(loc_dn, sys_dn)

loc_up_go <- topGO(loc_up, species = 'niben101')
loc_dn_go <- topGO(loc_dn, species = 'niben101')
sys_up_go <- topGO(sys_up, species = 'niben101')
sys_dn_go <- topGO(sys_dn, species = 'niben101')
olp_up_go <- topGO(olp_up, species = 'niben101')
olp_dn_go <- topGO(olp_dn, species = 'niben101')

go_lst <- list(loc_up=loc_up_go, loc_dn=loc_dn_go,
               sys_up=sys_up_go, sys_dn=sys_dn_go,
               olp_up=olp_up_go, olp_dn=olp_dn_go)

write2xlsx(go_lst, 'table/GO_enrichment.xlsx')

#########
all_mat <- combine_go_list(go_lst = go_lst, pvalue = 0.01)
all_mat[all_mat >=6 ] <- 6

pdf('figure/GO_enrich_heatmap.pdf', hei=13)
pheatmap(all_mat, cluster_cols=F, cluster_rows=F, border_color = 'gray', col=rev(colorRampPalette(brewer.pal(n=7, name="RdYlBu"))(100)), legend_breaks=c(2.1,3,4,5,6), legend_labels = c('2', '3', '4','5','>=6'))
dev.off()

# KEGG
library(clusterProfiler)
loc_up <- dfm$AthID[dfm$loc_TYLCV_vs_loc_EV == 1]
loc_dn <- dfm$AthID[dfm$loc_TYLCV_vs_loc_EV == -1]
sys_up <- dfm$AthID[dfm$sys_TYLCV_vs_sys_EV == 1]
sys_dn <- dfm$AthID[dfm$sys_TYLCV_vs_sys_EV == -1]
olp_up <- intersect(loc_up, sys_up)
olp_dn <- intersect(loc_dn, sys_dn)

loc_up_kegg <- as.data.frame(enrichKEGG(gene=loc_up, organism = 'ath', pvalueCutoff = 0.05))[c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
loc_dn_kegg <- as.data.frame(enrichKEGG(gene=loc_dn, organism = 'ath', pvalueCutoff = 0.05))[c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
sys_up_kegg <- as.data.frame(enrichKEGG(gene=sys_up, organism = 'ath', pvalueCutoff = 0.05))[c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
sys_dn_kegg <- as.data.frame(enrichKEGG(gene=sys_dn, organism = 'ath', pvalueCutoff = 0.05))[c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
olp_up_kegg <- as.data.frame(enrichKEGG(gene=olp_up, organism = 'ath', pvalueCutoff = 0.05))[c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]
olp_dn_kegg <- as.data.frame(enrichKEGG(gene=olp_dn, organism = 'ath', pvalueCutoff = 0.05))[c('Description', 'GeneRatio', 'BgRatio', 'p.adjust')]

kegg_lst <- list(loc_up=loc_up_kegg, loc_dn=loc_dn_kegg,
                 sys_up=sys_up_kegg, sys_dn=sys_dn_kegg,
                 olp_up=olp_up_kegg, olp_dn=olp_dn_kegg)

write2xlsx(kegg_lst, 'table/KEGG_enrichment.xlsx')

all_mat <- combine_kegg_list(kegg_lst = kegg_lst, pvalue = 0.05)
pdf('figure/KEGG_enrich_heatmap.pdf', height = 7)
pheatmap(all_mat, cluster_cols=F, cluster_rows=F, border_color = 'gray', col=rev(colorRampPalette(brewer.pal(n=7, name="RdYlBu"))(100)), legend_breaks=c(2.1,3,4,5,6), legend_labels = c('2', '3', '4','5','>=6'))
dev.off()
