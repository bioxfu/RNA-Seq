library(topGO)
library(GO.db)
library(xlsx)
library(yaml)

topGO <- function(myGenes, category='BP', p_cutoff=0.05, gomap, geneid){
  geneID2GO <- readMappings(file = gomap)
  geneNames <- read.table(geneid, stringsAsFactors=F)$V1
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
  return(allRes)
}

# set parameters
argv <- commandArgs(T)
config_file <- argv[1]
kmeans_all <- argv[2]
output_rdata <- argv[3]
output_xlsx <- argv[4]

# load data
config <- yaml.load_file(config_file)
k_cluster <- read.table(kmeans_all, head=T, stringsAsFactors=F)
gene_id <- rownames(k_cluster)

go <- list()
for (i in 1:config$kmean) {
  go[[i]] <- topGO(gene_id[k_cluster$cluster == i], gomap=config$gomap, geneid=config$geneid)
}

# add full term
for (i in 1:length(go)) {
  go[[i]]$Term <- apply(go[[i]], 1, function(x){Term(GOTERM[[x[1]]])})
}

save(go, file=output_rdata)

names(go) <- paste0('cluster', 1:config$kmean)
wb <- createWorkbook()
for (i in 1:length(go)) {
  dfm <- go[[i]]
  sheet <- createSheet(wb, sheetName=names(go)[i])
  addDataFrame(dfm, sheet, row.names = FALSE)
}
saveWorkbook(wb, output_xlsx)

