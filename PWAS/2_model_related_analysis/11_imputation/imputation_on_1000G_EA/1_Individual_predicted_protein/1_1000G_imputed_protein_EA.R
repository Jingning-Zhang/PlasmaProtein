
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(plink2R))

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

## Read protein .pos file in the provided package
prot <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/PWAS_tutorial/Plasma_Protein_EA_hg38.pos")

pred_matrix = matrix(nrow=498,ncol=nrow(prot)) # number of samples by number of proteins

## for each protein, impute the levels and save in the ith column in the pred_matrix
for (i in 1:nrow(prot)) {
    chr <- prot$CHR[i]

    ## Load the weights in the provided package
    load(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/coefs/",prot$WGT[i])) 

    wgt.matrix[is.na(wgt.matrix)] = 0
    m <- wgt.matrix[ , "enet" ]!=0
    wgt.matrix = wgt.matrix[m,,drop=F]
    snps = snps[m,,drop=F]
    writeLines(snps$V2,
               paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/tmp/",prot$WGT[i],"-",i,".txt"))

    system(
      paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2",
             " --bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh38/EUR/chr",chr,
             " --extract /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/tmp/",prot$WGT[i],"-",i,".txt",
             " --make-bed --threads 1",
             " --out /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/tmp/",i,"-",prot$WGT[i]), ignore.stdout=T)

    # load in genotype files by chromosome, restrict to matching SNPs and combine
    tmp <- tryCatch({ genos = read_plink(paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/tmp/",i,"-",prot$WGT[i]),impute="avg") }, error=function(e) e, warning=function(w) w)
    if(is(tmp,"warning")){ pred_matrix[,i] = 0; next }

    ## Please note it is very IMPORTANT to scale genotype matrix here, because the weights are for scaled genotype matrix!!
    sdtmp <- apply(genos$bed, MARGIN = 2, FUN=sd)
    a <- which(sdtmp==0)
    if(length(a)>0){
        genos$bed = scale(genos$bed[,-a]) ## scale genotype matrix
        genos$bim = genos$bim[-a,]
    }else{
        genos$bed = scale(genos$bed) ## scale genotype matrix
        genos$bim = genos$bim
    }

    # Match up the SNPs and weights
    m = match( snps[,2] , genos$bim[,2] )
    m.keep = !is.na(m)
    snps = snps[m.keep,,drop=F]

    wgt.matrix = wgt.matrix[m.keep,,drop=F]

    cur.genos = genos$bed[,m[m.keep],drop=F]
    cur.bim = genos$bim[m[m.keep],,drop=F]
    # Flip WEIGHTS for mismatching alleles
    qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
    wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

    # Predict into reference
    pred_matrix[,i] = cur.genos %*% wgt.matrix[ , "enet" ]

    rm(list = c("genos","sdtmp","a","cur.genos","cur.bim"))

    print(i)
}
colnames(pred_matrix) <- unlist(lapply(prot$WGT,FUN=function (x){substr(x,start=1, stop=nchar(x)-9)}))

write_tsv(as.data.frame(pred_matrix), paste0("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/1000G_imputed_prot.txt"))

system("rm -rf /dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/1000G_imputed/tmp/")
