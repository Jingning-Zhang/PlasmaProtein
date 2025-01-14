
rm(list=ls())

library(readr)
library(stringr)

r=2

for (i in 1:22){
  tmp <- try(read.table(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/conditional/",r,"/chr",i,"/conditional.txt")), silent=TRUE)
  if ('try-error' %in% class(tmp)) {
    print(paste0("chr",i,": no significant"))
  }else{
    tmp1 <- tmp[tmp$V19==1,]
    if(i==1){
      res <- tmp1
    }else{
      res <- rbind(res, tmp1)
    }
  }

  print(i)
}
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/Random_10times/conditional/",r,"/allsig.txt"), col_names = F)

length(unique(res$V1))
# 1486



#/dcl01/chatterj/data/jzhang2/TOOLS/plink/plink2 \
#--bfile /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/window2M_pre/byseq/SeqId_6919_3 \
#--pheno /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/invrankpheno/120/SeqId_6919_3.pheno \
#--glm allow-no-covars \
#--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results_GRCh38/White/pQTL/White_matchNblack/all_sample_peers/summary_stat/SeqId_6919_3
