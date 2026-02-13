#/bin/bash
awk \
-F, \
"BEGIN{OFS="\t"} NR == 1 || $5=="22A.." || $6=="22A.." {print $1,$3,$4,$7,$8,$9}" \
/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA/data/UKB/raw_data/Weight_gp_clinical.csv >\
 /ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA/data/UKB/raw_data/OnlyWeights_gp.tsv