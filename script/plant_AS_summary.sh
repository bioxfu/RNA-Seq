#VS=var2-9_vs_svr75
#VS=WT_vs_var2-T-DNA
VS=WT_vs_var2-9

echo -e "AS_event\tgene\tpvalue\tAS_type" > table/AS_summary_${VS}.tsv
grep -v 'pvalue' AS/all_sample_A3SS_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tA3SS"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv
grep -v 'pvalue' AS/all_sample_A5SS_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tA5SS"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv 
grep -v 'pvalue' AS/all_sample_AFE_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tAFE"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv
grep -v 'pvalue' AS/all_sample_ALE_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tALE"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv
grep -v 'pvalue' AS/all_sample_MXE_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tMXE"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv
grep -v 'pvalue' AS/all_sample_RI_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tRI"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv
grep -v 'pvalue' AS/all_sample_SE_fisher_test.tsv|awk '{if($2<0.05)print $1"\t"$2"\tSE"}'|sed -r 's/\.[0-9]\t/\t/'|sed 's/|AT/\tAT/' >> table/AS_summary_${VS}.tsv

Rscript script/plant_AS_anno.R config.yaml table/AS_summary_${VS}.tsv table/AS_summary_anno_${VS}.tsv

