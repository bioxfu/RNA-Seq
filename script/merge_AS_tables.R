library(yaml)

# set parameters
argv <- commandArgs(T)
config_file <- argv[1]

config <- yaml.load_file(config_file)
files <- config$samples
AS_type <- c('SE', 'RI', 'MXE', 'AFE', 'ALE', 'A3SS', 'A5SS')

calculate_inclusion_ratio <- function(x, AS) {
    
  if(AS == 'SE') {
    y = c((x[1] + x[2] + x[3]), x[4])
  }
  if(AS == 'RI') {
    y = x
  }
  if(AS == 'A5SS' || AS == 'A3SS') {
    y = c((x[1] + x[2]), x[3])
  }
  if(AS == 'AFE' || AS == 'ALE') {
    y = c((x[1] + x[2]), (x[3] + x[4]))
  }
  if(AS == 'MXE') {
    y = c((x[1] + x[2] + x[3]), (x[4] + x[5] + x[6]))
  }
  
  ratio <- round(y[1] / sum(y), 3)
  return(ratio)
}

for (AS in AS_type) {
  
  file <- files[1]
  tab1 <- read.table(paste0('AS/', file, '/', AS, '.txt'), head=T, sep='\t', row.names = 1)
  tab1$inclusion_ratio <- apply(tab1, 1, calculate_inclusion_ratio, AS)
  colnames(tab1) <- paste0(file,'_', colnames(tab1))
  
  for (file in files[-1]) {
    print(paste('merging', file, AS, '....'))
    tab <- read.table(paste0('AS/', file, '/', AS, '.txt'), head=T, sep='\t', row.names = 1)
    tab$inclusion_ratio <- apply(tab, 1, calculate_inclusion_ratio, AS)
    colnames(tab) <- paste0(file,'_', colnames(tab))
    
    tab1 <- merge(tab1, tab, by.x=0, by.y=0, all=T)
    rownames(tab1) <- tab1[,1]
    tab1 <- tab1[,-1]
  }

  counts <- rowSums(tab1[-grep('inclusion_ratio', colnames(tab1))])
  write.table(tab1[counts > 0,], paste0('AS/all_sample_',AS,'.tsv'), col.names=NA, sep='\t', quote=F)

}
