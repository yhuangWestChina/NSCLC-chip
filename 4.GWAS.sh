perl 1.pl ../all_data.fam ../covar_mds.txt >c1
perl 2.pl age.txt c1 >c2
perl 3.pl smoking.txt c2 >cov.txt
/its1/GB_BT2/huangliang/bio/plink/plink  --bfile ../all_data --ci 0.95  -covar cov.txt --logistic  --out logistic_all_result
/its1/GB_BT2/huangliang/bio/plink/plink  --bfile ../all_data --ci 0.95  -covar cov.txt --logistic --hide-covar --out logistic_result
perl d2-add.pl logistic_result.assoc.logistic >qq.data
perl /its1/GB_BT1/zhuangzhenhua/work/2021/chip/an/GWAS/sort_p.pl qq.data >sort
cp /its1/GB_BT1/zhuangzhenhua/project/2020/CLP/GWAS/figure/qq.r .
cp /its1/GB_BT1/zhuangzhenhua/project/2020/CLP/GWAS/figure/qq.sh .
perl g_m.pl logistic_result.assoc.logistic >m.data
cp /its1/GB_BT1/zhuangzhenhua/work/2021/chip/an/GWAS/test/manhattan.r .
/opt/bio/R/bin/Rscript manhattan.r
perl deal_result.pl logistic_result.assoc.logistic >r1
perl deal_result2.pl /its1/GB_BT1/zhuangzhenhua/work/2021/chip/db/ASA-CHIA-ANNOVA-GWAS.txt  r1 >GWAS.xls