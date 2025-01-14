library(readr)
annota <- read_tsv("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/prot.anno_autosomal.txt")

summary.AA <- readRDS("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/Black/PWAS/summary.black.rds")
summary.EA <- readRDS("/dcs04/nilanjan/data/jzhang2/pwas/pipeline/Results_GRCh38/White/PWAS/summary.white.rds")

tmp <- summary.AA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.AA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]
tmp <- summary.EA$RDat
tmp <- gsub(".wgt.RDat","",tmp)
summary.EA$gene <- annota$entrezgenesymbol[match(tmp,annota$seqid_in_sample)]

annota.Liver <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.Liver.geneid.pos"))
annota.Whole_Blood <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/WEIGHTS_v8.ALL/GTEXv8.Whole_Blood.geneid.pos"))

summary.Liver <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v8.ALL/summary.Liver.rds"))
summary.Whole_Blood <- readRDS("/dcs04/nilanjan/data/jzhang2/TWAS/fusion_twas-master/RESULTS_v8.ALL/summary.Whole_Blood.rds")

summary.Liver$gene <- annota.Liver$ID[match(paste0("Liver/",summary.Liver$RDat), annota.Liver$WGT)]
summary.Whole_Blood$gene <- annota.Whole_Blood$ID[match(paste0("Whole_Blood/",summary.Whole_Blood$RDat), annota.Whole_Blood$WGT)]

geneEALiver <- unique(intersect(summary.EA$gene, summary.Liver$gene))
geneEAWhole_Blood <- unique(intersect(summary.EA$gene, summary.Whole_Blood$gene))

geneAALiver <- unique(intersect(summary.AA$gene, summary.Liver$gene))
geneAAWhole_Blood <- unique(intersect(summary.AA$gene, summary.Whole_Blood$gene))

## AA Liver
df.hsqAALiver <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneAALiver]), 
                            group=c(rep("P in AA", sum(summary.AA$gene %in% geneAALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneAALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneAALiver)) ),
                            stringsAsFactors = F)

df.hsqAALiver$group <- factor(df.hsqAALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p1 <- ggplot(data = df.hsqAALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneAALiver)," overlapping cis-heritable\ngenes for expression in liver\nand protein in plasma in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  scale_fill_manual(values=c("#4a1486","#cb181d"))

ggsave(filename="AA_Liver.png", 
       plot=p1, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/cish2/", 
       width=3.7, height=5, units="in", dpi=500)


## AA Whole_Blood
df.hsqAAWhole_Blood <- data.frame(hsq = c(summary.AA$hsq.all[summary.AA$gene %in% geneAAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneAAWhole_Blood]), 
                                  group=c(rep("P in AA", sum(summary.AA$gene %in% geneAAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.AA$gene %in% geneAAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneAAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqAAWhole_Blood$group <- factor(df.hsqAAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p1 <- ggplot(data = df.hsqAAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneAAWhole_Blood)," overlapping cis-heritable\ngenes for expression in whole blood\nand protein in plasma in AA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  scale_fill_manual(values=c("#4a1486","#cb181d"))

ggsave(filename="AA_Whole_Blood.png", 
       plot=p1, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/cish2/", 
       width=3.7, height=5, units="in", dpi=500)




## EA Liver
df.hsqEALiver <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEALiver], 
                                    summary.Liver$hsq.all[summary.Liver$gene %in% geneEALiver]), 
                            group=c(rep("P in EA", sum(summary.EA$gene %in% geneEALiver)),
                                    rep("T in Liver", sum(summary.Liver$gene %in% geneEALiver)) ),
                            kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEALiver)), 
                                   rep("Gene expr-\nessions (T)", sum(summary.Liver$gene %in% geneEALiver)) ),
                            stringsAsFactors = F)

df.hsqEALiver$group <- factor(df.hsqEALiver$group, 
                              levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p1 <- ggplot(data = df.hsqEALiver, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneEALiver)," overlapping cis-heritable\ngenes for expression in liver\nand protein in plasma in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  scale_fill_manual(values=c("#4a1486","#cb181d"))

ggsave(filename="EA_Liver.png", 
       plot=p1, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/cish2/", 
       width=3.7, height=5, units="in", dpi=500)


## EA Whole_Blood
df.hsqEAWhole_Blood <- data.frame(hsq = c(summary.EA$hsq.all[summary.EA$gene %in% geneEAWhole_Blood], 
                                          summary.Whole_Blood$hsq.all[summary.Whole_Blood$gene %in% geneEAWhole_Blood]), 
                                  group=c(rep("P in EA", sum(summary.EA$gene %in% geneEAWhole_Blood)),
                                          rep("T in Whole_Blood", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  kind=c(rep("Plasma protein\nSOMAmers (P)", sum(summary.EA$gene %in% geneEAWhole_Blood)), 
                                         rep("Gene expr-\nessions (T)", sum(summary.Whole_Blood$gene %in% geneEAWhole_Blood)) ),
                                  stringsAsFactors = F)

df.hsqEAWhole_Blood$group <- factor(df.hsqEAWhole_Blood$group, 
                                    levels = c("T in Liver", "T in Whole_Blood", "P in AA", "P in EA"))
p1 <- ggplot(data = df.hsqEAWhole_Blood, aes(x = group)) + 
  geom_boxplot(alpha=0.6, notch = TRUE, notchwidth = 0.5, aes(y=hsq, fill=kind)) +
  coord_cartesian(ylim = c(0,0.5)) +  
  labs(y = expression(paste("cis-",h^2)), x=NULL, 
       title=paste0(length(geneEAWhole_Blood)," overlapping cis-heritable\ngenes for expression in whole blood\nand protein in plasma in EA")) +
  theme(legend.position="none",
        legend.title=element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(color = c("#4a1486","#cb181d")))+
  scale_fill_manual(values=c("#4a1486","#cb181d"))

ggsave(filename="EA_Whole_Blood.png", 
       plot=p1, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/sp/cish2/", 
       width=3.7, height=5, units="in", dpi=500)

