library(GO.db)
library(topGO)

gomap <- c(tair10='script/tair10_gene2go.map', niben101='script/Niben101_gene2go.map')
geneid <- c(tair10='script/tair10_gene2go.geneid', niben101='script/Niben101_gene2go.geneid')

### GO annotation
topGO <- function(myGenes, category='BP', p_cutoff=0.05, species='tair10'){
  geneID2GO <- readMappings(file = gomap[species])
  geneNames <- read.table(geneid[species], stringsAsFactors=F)$V1
  geneList <- factor(as.integer(geneNames %in% myGenes))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata", ontology=category, allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
  allRes <- GenTable(GOdata,pvalue=resultFisher,topNodes=100)
  allRes$pvalue[grep('<',allRes$pvalue)] <- "1e-30"
  allRes$pvalue <- as.numeric(allRes$pvalue)
  allRes <- allRes[order(allRes$pvalue,decreasing=F),]
  allRes$catagory <- category
  allRes <- allRes[allRes$pvalue < p_cutoff, ]
  allRes$Term <- apply(allRes, 1, function(x){Term(GOTERM[[x[1]]])})
  return(allRes)
}

combine_go_list <- function(go_lst, pvalue=0.01) {
  go_lst_filt <- lapply(go_lst, function(x){x[x$pvalue<pvalue,]})
  all_term <- unique(do.call(rbind, go_lst_filt)$Term)
  all_mat <- matrix(nrow=length(all_term), ncol=length(go_lst_filt), dimnames=list(all_term, names(go_lst)))
  for (i in 1:length(go_lst_filt)) {
    all_mat[go_lst_filt[[i]]$Term, i] <- -log10(go_lst_filt[[i]]$pvalue)
  }
  return(all_mat)
}
