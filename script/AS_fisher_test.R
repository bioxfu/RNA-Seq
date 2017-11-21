library(yaml)

AS_fisher_test <- function(dfm_ctrl, dfm_expt, AS, pref) {
  print(pref)
  if (AS == 'SE') {
    input <- cbind((dfm_ctrl[1] + dfm_ctrl[2] + dfm_ctrl[3]), dfm_ctrl[4], (dfm_expt[1] + dfm_expt[2] + dfm_expt[3]), dfm_expt[4])
  }
  if(AS == 'RI') {
    input <- cbind(dfm_ctrl, dfm_expt)
  }
  if(AS == 'A5SS' || AS == 'A3SS') {
    input <- cbind((dfm_ctrl[1] + dfm_ctrl[2]), dfm_ctrl[3], (dfm_expt[1] + dfm_expt[2]), dfm_expt[3])
  }
  if(AS == 'AFE' || AS == 'ALE') {
    input <- cbind((dfm_ctrl[1] + dfm_ctrl[2]), (dfm_ctrl[3] + dfm_ctrl[4]), (dfm_expt[1] + dfm_expt[2]), (dfm_expt[3] + dfm_expt[4]))
  }
  if(AS == 'MXE') {
    input <- cbind((dfm_ctrl[1] + dfm_ctrl[2] + dfm_ctrl[3]), (dfm_ctrl[4] + dfm_ctrl[5] + dfm_ctrl[6]), (dfm_expt[1] + dfm_expt[2] + dfm_expt[3]), (dfm_expt[4] + dfm_expt[5] + dfm_expt[6]))
  }
  input <- round(input)
  pvalue <- apply(input, 1, function(x) {
    if ( (x[1] + x[2]) > 0 & (x[3] + x[4]) > 0 ) {
      return(fisher.test(matrix(x, nrow=2, ncol=2))$p.value)
    }
    else {
      return(1)
    }})
  fdr <- p.adjust(pvalue, method='fdr')
  result <- data.frame(pvalue = round(pvalue, 4), fdr = round(fdr, 4))
  colnames(result) <- paste0(pref, '_', colnames(result))
  return(result)
}

#=====================
argv <- commandArgs(trailingOnly = T)
config_file <- argv[1]
input_file <- argv[2]
output_file <- argv[3]

# config_file <- '../config.yaml'
# input_file <- '../AS/all_sample_SE.tsv'

config <- yaml.load_file(config_file)
AS_type <- sub('.+all_sample_', '', input_file)
AS_type <- sub('.tsv', '', AS_type)

dfm <- read.table(input_file, head=T)
grp <- config$groups

test_result <- list()

for (vs in config$group_vs) {
  x <- strsplit(vs, ' ')[[1]]
  grp_ctrl <- grp[as.numeric(x[1])]
  grp_expt <- grp[as.numeric(x[2])]
  dfm_no_ratio <- dfm[, -grep('inclusion_ratio', colnames(dfm))]
  dfm_ctrl <- dfm_no_ratio[, grep(paste0(grp_ctrl, '.*_'), colnames(dfm_no_ratio))]
  dfm_expt <- dfm_no_ratio[, grep(paste0(grp_expt, '.*_'), colnames(dfm_no_ratio))]

  r <- config$seq_info$replicate
  
  if (r > 1) {
    n <- ncol(dfm_ctrl)
    l <- n/r
    s <- seq(1, n, by=l)
    dfm_ctrl2 <- dfm_ctrl[, 1:l]
    for (n in 2:length(s)) {
      dfm_ctrl2 <- dfm_ctrl2 + dfm_ctrl[, s[n]:(s[n]+l-1)]
    }
    dfm_ctrl <- dfm_ctrl2

    n <- ncol(dfm_expt)
    l <- n/r
    s <- seq(1, n, by=l)
    dfm_expt2 <- dfm_expt[, 1:l]
    for (n in 2:length(s)) {
      dfm_expt2 <- dfm_expt2 + dfm_expt[, s[n]:(s[n]+l-1)]
    }
    dfm_expt <- dfm_expt2
  }
  test_result <- c(test_result, AS_fisher_test(dfm_ctrl, dfm_expt, AS_type, paste0(grp_ctrl, '_', grp_expt)))
} 
  
dfm2 <- cbind(dfm, do.call(cbind, test_result))
write.table(dfm2, output_file, sep='\t', quote = F, col.names = NA)

