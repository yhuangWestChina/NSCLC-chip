plink --bfile  /its1/GB_BT1/zhuangzhenhua/work/2021/chip/an/all_data_9 --pca 20 --out pca
twstats  -t twtable -i pca.eigenval -o eigenvaltw.out