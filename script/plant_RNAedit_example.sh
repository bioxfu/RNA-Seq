
source ~/Python/REDItools/bin/activate
REDItoolDenovo.py -i bam/svr75_1.nodup.bam -f /home/xfu/Gmatic5/genome/tair10/tair10.fa -o redit/svr75_1
REDItoolDenovo.py -i bam/svr75_2.nodup.bam -f /home/xfu/Gmatic5/genome/tair10/tair10.fa -o redit/svr75_2
REDItoolDenovo.py -i bam/svr75_3.nodup.bam -f /home/xfu/Gmatic5/genome/tair10/tair10.fa -o redit/svr75_3
REDItoolDenovo.py -i bam/var2-9_1.nodup.bam -f /home/xfu/Gmatic5/genome/tair10/tair10.fa -o redit/var2-9_1
REDItoolDenovo.py -i bam/var2-9_2.nodup.bam -f /home/xfu/Gmatic5/genome/tair10/tair10.fa -o redit/var2-9_2
REDItoolDenovo.py -i bam/var2-9_3.nodup.bam -f /home/xfu/Gmatic5/genome/tair10/tair10.fa -o redit/var2-9_3


cat redit/svr75_1/denovo*/outTable* |grep -v 'Region'|awk -F "\t" '{if($9>0 && $10<0.05)print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > redit/RNA_edit_svr75_1 
cat redit/svr75_2/denovo*/outTable* |grep -v 'Region'|awk -F "\t" '{if($9>0 && $10<0.05)print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > redit/RNA_edit_svr75_2
cat redit/svr75_3/denovo*/outTable* |grep -v 'Region'|awk -F "\t" '{if($9>0 && $10<0.05)print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > redit/RNA_edit_svr75_3 
cat redit/var2-9_1/denovo*/outTable*|grep -v 'Region'|awk -F "\t" '{if($9>0 && $10<0.05)print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > redit/RNA_edit_var2-9_1 
cat redit/var2-9_2/denovo*/outTable*|grep -v 'Region'|awk -F "\t" '{if($9>0 && $10<0.05)print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > redit/RNA_edit_var2-9_2 
cat redit/var2-9_3/denovo*/outTable*|grep -v 'Region'|awk -F "\t" '{if($9>0 && $10<0.05)print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > redit/RNA_edit_var2-9_3 


bedtools intersect -a redit/RNA_edit_svr75_1 -b /home/xfu/Gmatic5/genome/tair10/tair10_gene.bed -wa -wb|cut -f1,3,9-11,15|sed 's/\t/|/' > redit/RNA_edit_svr75_1_gene
bedtools intersect -a redit/RNA_edit_svr75_2 -b /home/xfu/Gmatic5/genome/tair10/tair10_gene.bed -wa -wb|cut -f1,3,9-11,15|sed 's/\t/|/' > redit/RNA_edit_svr75_2_gene
bedtools intersect -a redit/RNA_edit_svr75_3 -b /home/xfu/Gmatic5/genome/tair10/tair10_gene.bed -wa -wb|cut -f1,3,9-11,15|sed 's/\t/|/' > redit/RNA_edit_svr75_3_gene
bedtools intersect -a redit/RNA_edit_var2-9_1 -b /home/xfu/Gmatic5/genome/tair10/tair10_gene.bed -wa -wb|cut -f1,3,9-11,15|sed 's/\t/|/' > redit/RNA_edit_var2-9_1_gene
bedtools intersect -a redit/RNA_edit_var2-9_2 -b /home/xfu/Gmatic5/genome/tair10/tair10_gene.bed -wa -wb|cut -f1,3,9-11,15|sed 's/\t/|/' > redit/RNA_edit_var2-9_2_gene
bedtools intersect -a redit/RNA_edit_var2-9_3 -b /home/xfu/Gmatic5/genome/tair10/tair10_gene.bed -wa -wb|cut -f1,3,9-11,15|sed 's/\t/|/' > redit/RNA_edit_var2-9_3_gene

cat redit/RNA_edit_svr75_*_gene|sort -k1,1|groupBy -g 1 -c 2,3,4,5 -o distinct,max,min,first|awk '{print $0"\tsvr75"}' > table/RNA_edit_svr75_gene
cat redit/RNA_edit_var2-9_*_gene|sort -k1,1|groupBy -g 1 -c 2,3,4,5 -o distinct,max,min,first|awk '{print $0"\tvar2-9"}' > table/RNA_edit_var2-9_gene

cat table/RNA_edit_var2-9_gene table/RNA_edit_svr75_gene|sort -k1,1|groupBy -g 1 -c 2,3,4,5,6 -o collapse,collapse,collapse,collapse,collapse|grep -v ',' > table/RNA_edit_var2-9_vs_svr75

Rscript script/plant_RNAedit_anno.R config.yaml table/RNA_edit_var2-9_vs_svr75 table/RNA_edit_var2-9_vs_svr75_anno.tsv

