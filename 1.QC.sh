# 过滤缺失率大于0.1的样本和SNP位点（删掉12个）
plink  --bfile all_data --missing
plink --bfile  all_data --geno 0.1 --make-bed --out all_data_1
plink --bfile  all_data_1 --mind 0.1 --make-bed --out all_data_2
# 过滤性别不一致的样本（删掉4个）
plink --bfile all_data_2 --check-sex 
grep "PROBLEM"  plink.sexcheck | awk '{print$1,$2}'> sex_discrepancy.txt
plink --bfile  all_data_2 --remove sex_discrepancy.txt --make-bed --out all_data_3
# 过滤MAF<0.01位点，并提取1~22号染色体位点
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' all_data_3.bim > snp_1_22.txt
plink --bfile  all_data_3 --extract snp_1_22.txt --make-bed --out all_data_4
plink --bfile  all_data_4 --freq --out MAF_check
Rscript --no-save MAF_check.R
plink --bfile  all_data_4 --maf 0.01 --make-bed --out all_data_5
# 删除HWE p-value<0.00001的位点
plink --bfile  all_data_5 --hardy
awk '{ if ($9 <0.00001) print $0 }' plink.hwe >plinkzoomhwe.hwe
Rscript --no-save hwe.R
plink --bfile  all_data_5 --hwe 1e-5 --make-bed --out all_data_6
# 根据LD删除杂合率偏离均值超过3倍标准差的个体（去掉8个样本）
plink --bfile  all_data_6  --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile  all_data_6   --extract indepSNP.prune.in --het --out R_check
Rscript --no-save check_heterozygosity_rate.R
Rscript --no-save heterozygosity_outliers_list.R
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
plink --bfile  all_data_6 --remove het_fail_ind.txt --make-bed --out all_data_7
# IBD计算样本之间亲缘关系，过滤IBD>0.1875样本（去掉4个样本）
plink --bfile  all_data_7  --extract indepSNP.prune.in --genome --min 0.1875 --out pihat_min0.1875
awk '{ if ($8 >0.9) print $0 }' pihat_min0.1875.genome >zoom_pihat.genome
plink --bfile  all_data_7 --filter-founders --make-bed --out all_data_8
plink --bfile  all_data_8  --extract indepSNP.prune.in --genome --min 0.1875 --out pihat_min0.1875_in_founders
plink --bfile  all_data_8 --missing
plink --bfile  all_data_8 --remove 0.1875_low_call_rate_pihat.txt --make-bed --out all_data_9
####
finally ,494929 variants and 726 people pass filters and QC，includ 275 cases and 451 controls.
