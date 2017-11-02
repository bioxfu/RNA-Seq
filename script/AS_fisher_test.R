#! /usr/bin/Rscript

args = commandArgs(trailingOnly = T)
AS_type = c('SE', 'A5SS', 'A3SS', 'AFE', 'ALE', 'MXE', 'RI')

#args = c('../../AS/W1', '../../AS/W2', 'W1', 'W2')
#AS='RI'

for (AS in AS_type) {
    print(AS)
    input1 = read.table(paste(args[1],'/',AS,'.txt',sep=''), header = T, row.names = 1)
    input2 = read.table(paste(args[2],'/',AS,'.txt',sep=''), header = T, row.names = 1)
    colnames(input1) = paste(args[3], '_', colnames(input1), sep='')
    colnames(input2) = paste(args[4], '_', colnames(input2), sep='')
    input3 = round(cbind(input1, input2))
    
    if(AS == 'SE') {
        input = cbind((input1[1] + input1[2] + input1[3]), input1[4], (input2[1] + input2[2] + input2[3]), input2[4])
    }
    if(AS == 'RI') {
        input = input3
    }
    if(AS == 'A5SS' || AS == 'A3SS') {
        input = cbind((input1[1] + input1[2]), input1[3], (input2[1] + input2[2]), input2[3])
    }
    if(AS == 'AFE' || AS == 'ALE') {
        input = cbind((input1[1] + input1[2]), (input1[3] + input1[4]), (input2[1] + input2[2]), (input2[3] + input2[4]))
    }
    if(AS == 'MXE') {
        input = cbind((input1[1] + input1[2] + input1[3]), (input1[4] + input1[5] + input1[6]), (input2[1] + input2[2] + input2[3]), (input2[4] + input2[5] + input2[6]))
    }
    
    filter = (input[1]+input[2]) > 0 & (input[3]+input[4]) > 0
    input = input[filter, ]
    input3 = input3[filter, ]
    input3$diff_ratio = round(input[,3]/(input[,3]+input[,4]) - input[,1]/(input[,1]+input[,2]),3)
    pvalue = apply(input, 1, function(x) {fisher.test(matrix(x,nrow=2,ncol=2))$p.value})
    fdr = p.adjust(pvalue,method='fdr')
    input3$pvalue = pvalue
    input3$fdr = fdr
    output = input3[order(input3$pvalue), ]
    output_up = output[output$diff_ratio > 0, ]
    output_dn = output[output$diff_ratio < 0, ]
    output_up$change = 'up'
    output_dn$change = 'dn'
    output_all = rbind(output_up, output_dn)
    output_sig = output_all[output_all$pvalue < 0.05,]
    output_nonsig = output_all[output_all$pvalue > 0.95,]
    output_sig$pvalue = sprintf('%.1e', output_sig$pvalue)
    output_sig$fdr = sprintf('%.1e', output_sig$fdr)
    output_nonsig$pvalue = sprintf('%.1e', output_nonsig$pvalue)
    output_nonsig$fdr = sprintf('%.1e', output_nonsig$fdr)
    write.table(output_sig, paste(args[5],'_',AS,'_sig.txt',sep=''), quote = F, col.names = NA, sep = '\t')
    write.table(output_all, paste(args[5],'_',AS,'_all.txt',sep=''), quote = F, col.names = NA, sep = '\t')
    #write.table(output_nonsig, paste(args[5],'_',AS,'_nonsig.txt',sep=''), quote = F, col.names = NA, sep = '\t')
    #output_all = input3[abs(input3$diff_ratio) > 0.1, ]
    #output_all = output_all[order(output_all$diff_ratio), ]
    #write.table(output_all, paste(args[4], '_diff_ratio_0.1.txt', sep=''), quote = F, col.names = NA, sep = '\t')
}


