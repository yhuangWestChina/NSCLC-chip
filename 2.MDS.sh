###分染色体下载1000g数据库vcf文件,合并
#for id in `seq 1 21` X; do axel -n 40 http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$id.1kg.phase3.v5a.vcf.gz; done 
#/opt/bio/vcftools_0.1.13/bin/vcf-concat chr1.1kg.phase3.v5a.vcf.gz chr2.1kg.phase3.v5a.vcf.gz chr3.1kg.phase3.v5a.vcf.gz chr4.1kg.phase3.v5a.vcf.gz chr5.1kg.phase3.v5a.vcf.gz chr6.1kg.phase3.v5a.vcf.gz chr7.1kg.phase3.v5a.vcf.gz chr8.1kg.phase3.v5a.vcf.gz chr9.1kg.phase3.v5a.vcf.gz chr10.1kg.phase3.v5a.vcf.gz chr11.1kg.phase3.v5a.vcf.gz chr12.1kg.phase3.v5a.vcf.gz chr13.1kg.phase3.v5a.vcf.gz chr14.1kg.phase3.v5a.vcf.gz chr15.1kg.phase3.v5a.vcf.gz chr16.1kg.phase3.v5a.vcf.gz chr17.1kg.phase3.v5a.vcf.gz chr18.1kg.phase3.v5a.vcf.gz chr19.1kg.phase3.v5a.vcf.gz chr20.1kg.phase3.v5a.vcf.gz chr21.1kg.phase3.v5a.vcf.gz chr22.1kg.phase3.v5a.vcf.gz > 1000g.vcf
###去除重复和indel位点
perl filter_indel.pl 1000g.vcf >1000g2.vcf
perl filterID.pl 1000g2.vcf >1000g3.vcf
##转plink格式
plink --vcf 1000g3.vcf --make-bed --out ALL.2of4intersection.20100804.genotypes
### 处理没rs编号位点
plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs
### 去掉缺失>0.2位点
plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS
### 样本和位点进行QC
plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2
plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3
plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4
plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5
### 提取本项目rs位点
awk '{print$2}' HapMap_3_r3_12.bim > HapMap_SNPs.txt
plink  --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6
#### 提取千人rs位点
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
plink --bfile HapMap_3_r3_12 --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out HapMap_MDS
#####  重构map文件
awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
plink  --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
#### 整合千人与本项目位点集合
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
plink  --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp
sort 1kGMDS7_tmp HapMap-adj_tmp |uniq -u > all_differences.txt
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS2
plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8
plink --bfile HapMap_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2
####MDS分析
plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2
plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2
awk '{print$1,$1,$3}' integrated_call_samples_v3.20130502.ALL.panel.txt  > race_1kG.txt
awk '{print$1,$2,"OWN"}' HapMap_MDS.fam >racefile_own.txt
cat race_1kG.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt
Rscript MDS_merged.R
####剔除与东亚分群样本
awk '{ if ($4 <-0.04 && $5 >0.02) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2 
plink --bfile HapMap_3_r3_12 --keep EUR_MDS_merge2 --make-bed --out HapMap_3_r3_13
### 创建GWAS分析协变量文件
plink --bfile HapMap_3_r3_13 --extract indepSNP.prune.in --genome --out HapMap_3_r3_13
plink --bfile HapMap_3_r3_13 --read-genome HapMap_3_r3_13.genome --cluster --mds-plot 10 --out HapMap_3_r3_13_mds
awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' HapMap_3_r3_13_mds.mds > covar_mds.txt