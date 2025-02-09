
rm(list = ls())

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(plink2R))

args <- commandArgs(T)

for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

annota.t <- read_tsv(paste0('/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR/GTEXv8.',tissue,'.pos'))
annota.t <- annota.t[!is.na(annota.t$CHR),]

geneENSG <- substr(annota.t$WGT, start = 2*nchar(annota.t$PANEL[1])+3, stop = nchar(annota.t$WGT)-9)

allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F

  snp[["flip"]] = a1 != ref1

  return(snp)
}

allgene <- length(geneENSG)

pred_matrix = matrix(nrow=498,ncol=allgene)

for (i in 1:allgene) {

  WGTi <- paste0(annota.t$PANEL[1], "/", annota.t$PANEL[1],".", geneENSG[i],".wgt.RDat")
  WGTi <- gsub("GTExv8","GTEXv8",WGTi)
  chr <- annota.t$CHR[annota.t$WGT==WGTi]

  load(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.EUR/", WGTi))
  wgt.matrix[is.na(wgt.matrix)] = 0
  m <- wgt.matrix[ , "enet" ]!=0

  if(sum(m)==0){
    pred_matrix[,i] = 0
  }else{
    wgt.matrix = wgt.matrix[m,,drop=F]
    snps = snps[m,,drop=F]
    writeLines(snps$V2,
               paste0("/fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i,".txt"))

    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2",
                  " --bfile /fastscratch/myscratch/jzhang2/1000G/GRCh38/EUR/chr",chr,
                  " --extract /fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i,".txt",
                  " --make-bed",
                  " --out /fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i),
           ignore.stdout = TRUE)

    # load in genotype files by chromosome, restrict to matching SNPs and combine
    tmp <- tryCatch({genos = read_plink(paste0("/fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i),impute="avg") }, error=function(e) e, warning=function(w) w)
    if(is(tmp,"warning")){ pred_matrix[,i] = 0; next }
    sdtmp <- apply(genos$bed, MARGIN = 2, FUN=sd)
    a <- which(sdtmp==0)
    if(length(a)>0){
      genos$bed = scale(genos$bed[,-a])
      genos$bim = genos$bim[-a,]
      if(nrow(genos$bim)==0){
        pred_matrix[,i] = 0
        system(paste0("rm -rf /fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i,".*"))
        next
      }
    }else{
      genos$bed = scale(genos$bed)
    }

    if(sum(m)==1){

      # Match up the SNPs and weights
      m = match( snps[,2] , genos$bim[,2] )
      m.keep = !is.na(m)

      if(sum(m.keep)==0){
        pred_matrix[,i] = 0
        system(paste0("rm -rf /fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i,".*"))
        next
      }

      cur.genos = genos$bed[,m[m.keep]]
      cur.bim = genos$bim[m[m.keep],]
      # Flip WEIGHTS for mismatching alleles
      qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
      if(qc$flip){
        wgt.matrix = -1 * wgt.matrix
      }

      # Predict into reference
      pred_matrix[,i] = cur.genos * wgt.matrix[ "enet" ]

    }else{

      # Match up the SNPs and weights
      m = match( snps[,2] , genos$bim[,2] )
      m.keep = !is.na(m)

      if(sum(m.keep)==0){
        pred_matrix[,i] = 0
        system(paste0("rm -rf /fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i,".*"))
        next
      }

      snps = snps[m.keep,]
      wgt.matrix = wgt.matrix[m.keep,]

      cur.genos = genos$bed[,m[m.keep]]
      cur.bim = genos$bim[m[m.keep],]
      # Flip WEIGHTS for mismatching alleles
      qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )

      if(sum(m.keep)==1){
        if(qc$flip){
          wgt.matrix = -1 * wgt.matrix
        }

        # Predict into reference
        pred_matrix[,i] = cur.genos * wgt.matrix[ "enet" ]

      }else{
        wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

        # Predict into reference
        pred_matrix[,i] = cur.genos %*% wgt.matrix[ , "enet" ]
      }

      rm(list = c("genos","sdtmp","a","cur.genos","cur.bim"))
    }
  }
  system(paste0("rm -rf /fastscratch/myscratch/jzhang2/PWAS/1000G_imputed_FUSION_v8.EUR/tmp/",tissue,"-",i,".*"))

  print(paste0(i,"/",allgene))
}

colnames(pred_matrix) <- geneENSG

write_tsv(as.data.frame(pred_matrix), paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_FUSION_v8.EUR/",tissue,".txt"))

