library(dplyr)
library(tidyr)
library(tidyverse)
library("ggplot2")
library(ggpubr)
library(enrichplot)
library(clusterProfiler)
library(rstatix)


RNA_dir<-"ED_Transcriptome_evol/output_group/"
DNA_dir<-"ED_resequencing/"
DNA_RNA_dir<-"ED_evol_SNP_SV_RNA/"
setwd(DNA_dir)

### Read in GO annotation and format
ED_gff_GO<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_GO.genename.csv",
                    header=TRUE)
colnames(ED_gff_GO)<-c("gene","gotermId","goName", "keytypes","goAcc")
ED_gff_GO$gene<-as.character(gsub(".*Exode1.","",ED_gff_GO$gene))

ED_gff_GO$keytypes[ED_gff_GO$keytypes=="biological_process"]<-"BP"
ED_gff_GO$keytypes[ED_gff_GO$keytypes=="cellular_component"]<-"CC"
ED_gff_GO$keytypes[ED_gff_GO$keytypes=="molecular_function"]<-"MF"

ED_gff_GO$gene[ED_gff_GO$gene=="2864"]<-"rad52" 
ED_gff_GO$gene[ED_gff_GO$gene=="3244"]<-"pks1"  
ED_gff_GO$gene[ED_gff_GO$gene=="8372"]<-"8115"

### Read in Interproscan annotation
ED_gff_IPR<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_IPR.csv",
                     header=TRUE)
colnames(ED_gff_IPR)[1]<-"gene"
ED_gff_IPR$gene<-as.character(ED_gff_IPR$gene)

ED_gff_IPR$gene[ED_gff_IPR$gene=="2864"]<-"rad52" 
ED_gff_IPR$gene[ED_gff_IPR$gene=="3244"]<-"pks1"  
ED_gff_IPR$gene[ED_gff_IPR$gene=="8372"]<-"8115"
### Read in KEGG annotation
ED_gff_KEGG<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_KEGG.csv",
                      header=TRUE)
colnames(ED_gff_KEGG)[1]<-"gene"
ED_gff_KEGG$gene<-as.character(ED_gff_KEGG$gene)

ED_gff_KEGG$gene[ED_gff_KEGG$gene=="2864"]<-"rad52" 
ED_gff_KEGG$gene[ED_gff_KEGG$gene=="3244"]<-"pks1"  
ED_gff_KEGG$gene[ED_gff_KEGG$gene=="8372"]<-"8115"
### Read in KOG annoation
ED_gff_KOG<-read.csv("../ED_Transcriptome/Exode1_GeneCatalog_proteins_20160902_KOG.csv",
                     header=TRUE)
colnames(ED_gff_KOG)[1]<-"gene"
ED_gff_KOG$gene<-as.character(ED_gff_KOG$gene)

ED_gff_KOG$gene[ED_gff_KOG$gene=="2864"]<-"rad52" 
ED_gff_KOG$gene[ED_gff_KOG$gene=="3244"]<-"pks1"  
ED_gff_KOG$gene[ED_gff_KOG$gene=="8372"]<-"8115"

term2gene_go<-ED_gff_GO[,c("goName","gene", "keytypes")]
term2gene_kegg<-ED_gff_KEGG[,c("gene", "pathway", "pathway_class")]
term2gene_kog<-ED_gff_KOG[,c("gene", "kogdefline", "kogClass")]

## Merge the AnnotationHub dataframe with the results 
ed_gff_go<-mutate(ED_gff_GO[c("gene","goName")])%>%
  group_by(gene)%>% summarise(goName = paste0(goName, collapse = ";")) %>%
  ungroup()
ed_gff_kegg<-mutate(ED_gff_KEGG[c("gene","pathway")])%>%
  group_by(gene)%>% summarise(pathway = paste0(pathway, collapse = ";")) %>%
  ungroup()
ed_gff_kog<-mutate(ED_gff_KOG[c("gene","kogdefline")])%>%
  group_by(gene)%>% summarise(kogdefline = paste0(kogdefline, collapse = ";")) %>%
  ungroup()
ed_gff_ipr<-mutate(ED_gff_IPR[c("gene","iprDesc")],)%>%
  group_by(gene)%>% summarise(iprDesc = paste0(iprDesc, collapse = ";")) %>%
  ungroup()

annot<-Reduce(function(...) merge(..., by='gene', all=TRUE), 
              list(ed_gff_go,ed_gff_kegg,ed_gff_kog,ed_gff_ipr))
genenames<-read.table("../ED_Transcriptome_evol/transcript_genenames.txt",
                      header=FALSE,col.names=c("gene", "genename"))
annot<-merge(genenames,annot, by="gene", all=TRUE)

#############################################
#### Read in SNP data #######################
#############################################
### read in SNP dataframe
SNPs_df0 <- as.data.frame(read.table("ED.multisample.WT_vs_PKS.variant.GT.txt", header = TRUE))
### read in SNPeff dataframe
SNPeff_vcf<-as.data.frame(read.table("ED.multisample.WT_vs_PKS.variant.snpeff.txt", header = TRUE)) #change header: gene	Type	EFF 
SNPeff <- as.data.frame(read.table("ED.multisample.WT_vs_PKS.variant.vcf.summary.genes.txt", header = TRUE))
colnames(SNPeff)[2]<-"gene"

SNPeff<-merge(SNPeff_vcf,SNPeff, by='gene', all=TRUE)
SNPeff$gene<-gsub("gene_","",SNPeff$gene)
SNPeff$gene[SNPeff$gene=="2864"]<-"rad52" 
SNPeff$gene[SNPeff$gene=="3244"]<-"pks1"  
SNPeff$gene[SNPeff$gene=="8372"]<-"8115"
SNPeff<-merge(SNPeff,annot,by="gene")
colnames(SNPeff)[84:97]<-paste0("orig-",colnames(SNPeff)[84:97])
colnames(SNPs_df0)[81:94]<-paste0("orig-",colnames(SNPs_df0)[81:94])

##################################################
#### Read in SV data #############################
##################################################
svcallers<-read.table("sample_merged.GT.vcf",header=TRUE)

colnames(svcallers)<-gsub(".*W","W",colnames(svcallers))
colnames(svcallers)<-gsub(".*P","P",colnames(svcallers))
colnames(svcallers)<-gsub(".*CHROM","CHROM",colnames(svcallers))
colnames(svcallers)<-gsub(".GT","",colnames(svcallers))
colnames(svcallers)[85:96]<-paste0("orig-",colnames(svcallers)[85:96])

###add genotype information
svcallers[svcallers=="./."] <- 0
svcallers[svcallers=="1/1"] <- 1
svcallers[6:96]<-mutate_all(svcallers[6:96], function(x) as.numeric(as.character(x)))

##################################################
################# SNP counts #####################
##################################################
my_comparisons <- list( c("WTC", "WT15"), c("PKSC", "PKS15"), 
                        c("WTC", "PKSC"), c("WT15", "PKS15"))
###Calculate frequencies of genotypes per Group
SNPs_sample <- data.frame(lapply(SNPs_df0, as.character), stringsAsFactors=FALSE)
SNPs_sample <- lapply(SNPs_sample[grep("WT|PKS",colnames(SNPs_sample))],table)
SNPs_sample <-plyr::ldply(SNPs_sample, rbind)
SNPs_sample[4][is.na(SNPs_sample[4])]<-0 ##no genotype 2
SNPs_sample$genotyped<-rowSums(SNPs_sample[3:4],na.rm=TRUE)
colnames(SNPs_sample)<-c("Sample","mis","count_0","count_1","genotyped")
SNPs_sample$alt<-unlist(SNPs_sample[4])
SNPs_sample$freq_0<-unlist(SNPs_sample[3]/SNPs_sample$genotyped)
SNPs_sample$freq_1<-unlist(SNPs_sample[4]/SNPs_sample$genotyped)
SNPs_sample$freq_alt<-unlist((SNPs_sample[4])/SNPs_sample$genotyped)

SNPs_sample$condition<-gsub("orig.*C.*","orig-control",SNPs_sample$Sample)
SNPs_sample$condition<-gsub("orig.*15.*","orig-evolved",SNPs_sample$condition)
SNPs_sample$condition<-gsub(".*C.*","control",SNPs_sample$condition)
SNPs_sample$condition<-gsub(".*15.*","evolved",SNPs_sample$condition)

SNPs_sample$Group<-gsub("^WT15.*","WT15",SNPs_sample$Sample)
SNPs_sample$Group<-gsub("^PKS15.*","PKS15",SNPs_sample$Group)
SNPs_sample$Group<-gsub("^WTC.*","WTC",SNPs_sample$Group)
SNPs_sample$Group<-gsub("^PKSC.*","PKSC",SNPs_sample$Group)
SNPs_sample$Group<-gsub("^orig.WTC.*","orig-WTC",SNPs_sample$Group)
SNPs_sample$Group<-gsub("^orig.PKSC.*","orig-PKSC",SNPs_sample$Group)
SNPs_sample$Group<-gsub("^orig.WT15.*","orig-WT15",SNPs_sample$Group)
SNPs_sample$Group<-gsub("^orig.PKS15.*","orig-PKS15",SNPs_sample$Group)
SNPs_sample$Group<-factor(SNPs_sample$Group, levels = c("WTC", "WT15","PKSC","PKS15",
                                                        "orig-WTC", "orig-WT15","orig-PKSC","orig-PKS15"))
SNPs_sample$lineage<-SNPs_sample$Sample
SNPs_sample$lineage[SNPs_sample$Group=="WT15"|SNPs_sample$Group=="PKS15"]<-gsub('\\..$',"",SNPs_sample$lineage[SNPs_sample$Group=="WT15"|SNPs_sample$Group=="PKS15"])

SNPs_sample_2<-aggregate(SNPs_sample$alt, by=list(SNPs_sample$lineage), FUN=mean)
SNPs_sample_2$Group<-gsub("^WT15.*","WT15",SNPs_sample_2$Group.1)
SNPs_sample_2$Group<-gsub("^PKS15.*","PKS15",SNPs_sample_2$Group)
SNPs_sample_2$Group<-gsub("^WTC.*","WTC",SNPs_sample_2$Group)
SNPs_sample_2$Group<-gsub("^PKSC.*","PKSC",SNPs_sample_2$Group)
SNPs_sample_2$Group<-gsub("^orig.WTC.*","orig-WTC",SNPs_sample_2$Group)
SNPs_sample_2$Group<-gsub("^orig.PKSC.*","orig-PKSC",SNPs_sample_2$Group)
SNPs_sample_2$Group<-gsub("^orig.WT15.*","orig-WT15",SNPs_sample_2$Group)
SNPs_sample_2$Group<-gsub("^orig.PKS15.*","orig-PKS15",SNPs_sample_2$Group)
SNPs_sample_2$condition<-gsub("orig.*C.*","orig-control",SNPs_sample_2$Group.1)
SNPs_sample_2$condition<-gsub("orig.*15.*","orig-evolved",SNPs_sample_2$condition)
SNPs_sample_2$condition<-gsub(".*C.*","control",SNPs_sample_2$condition)
SNPs_sample_2$condition<-gsub(".*15.*","evolved",SNPs_sample_2$condition)

colnames(SNPs_sample_2)<-c("lineage","alt","Group","condition")

########## How many SNPs missing per isolate ###########
summary(SNPs_sample$mis)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#1.00   15.75   56.50  163.99  166.25 1232.00      16  
summary(SNPs_sample$genotyped)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1510    2612    2712    2607    2735    2742 

### Plot per sample SNP count
snp_counts_sample<-ggplot(SNPs_sample, aes(x=Sample, y=alt, fill=condition)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(breaks = c("control", "evolved",
                               "orig-control","orig-evolved"), 
                    values=c("blue", "red","dark blue","dark red"))+
  coord_flip() +ylab("SNP counts") #+ scale_x_discrete(limits=rev)


##################################################
#### SV counts ###################################
##################################################
### Calculate stats for barplots
svcallers.df<-data.frame()
for (svtype in unique(svcallers$PE)) {
  m<-svcallers[svcallers$PE==svtype,]
  svcallers_stats<-data.frame(colSums(m[6:96]))
  svcallers_stats$sample<-rownames(svcallers_stats)
  colnames(svcallers_stats)<-c("counts","sample")
  svcallers_stats$isolate<-gsub("15.*","15",svcallers_stats$sample)
  svcallers_stats$isolate<-gsub("C.*","C",svcallers_stats$isolate)
  svcallers_stats$svtype<-svtype
  svcallers.df<-rbind(svcallers.df,data.frame(svcallers_stats))
}
svcallers.df$condition<-svcallers.df$isolate
svcallers.df$condition[svcallers.df$condition=="PKSC" | svcallers.df$condition=="WTC"]<-"control"
svcallers.df$condition[svcallers.df$condition=="PKS15" | svcallers.df$condition=="WT15"]<-"evolved"
svcallers.df$condition[svcallers.df$condition=="orig-PKSC" | svcallers.df$condition=="orig-WTC"]<-"orig-control"
svcallers.df$condition[svcallers.df$condition=="orig-PKS15" | svcallers.df$condition=="orig-WT15"]<-"orig-evolved"
svcallers.df$condition <- factor(svcallers.df$condition, levels = c("control","evolved","orig-control","orig-evolved"))

my_comparisons <- list( c("WTC", "WT15"), c("PKSC", "PKS15"), 
                        c("WTC", "PKSC"), c("WT15", "PKS15"))

mean(svcallers.df$counts[svcallers.df$svtype=="TRA"])
mean(svcallers.df$counts[svcallers.df$svtype=="INV"])
mean(svcallers.df$counts[svcallers.df$svtype=="DUP"])
mean(svcallers.df$counts[svcallers.df$svtype=="DEL"])

svcallers.df_orig<-svcallers.df[svcallers.df$condition=="orig-evolved"|svcallers.df$condition=="orig-control",]
ggplot(svcallers.df_orig, 
       aes(x=svtype, y=counts)) + 
  geom_boxplot()+
  coord_flip()+
  stat_compare_means(comparisons = list( c("TRA", "INV"), c("TRA", "DUP"),c("TRA", "DEL"), 
                                         c("INV", "DUP"),c("INV", "DEL"),c("DUP","DEL")), 
                     method="wilcox.test")
mean(svcallers.df_orig$counts[svcallers.df_orig$svtype=="TRA"])
mean(svcallers.df_orig$counts[svcallers.df_orig$svtype=="INV"])
mean(svcallers.df_orig$counts[svcallers.df_orig$svtype=="DUP"])
mean(svcallers.df_orig$counts[svcallers.df_orig$svtype=="DEL"])

svcallers.df_new<-svcallers.df[svcallers.df$condition=="evolved"|svcallers.df$condition=="control",]
ggplot(svcallers.df_new, 
       aes(x=svtype, y=counts)) + 
  geom_boxplot()+
  coord_flip()+
  stat_compare_means(comparisons = list( c("TRA", "INV"), c("TRA", "DUP"),c("TRA", "DEL"), 
                                         c("INV", "DUP"),c("INV", "DEL"),c("DUP","DEL")), 
                     method="wilcox.test")
mean(svcallers.df_new$counts[svcallers.df_new$svtype=="TRA"])
mean(svcallers.df_new$counts[svcallers.df_new$svtype=="INV"])
mean(svcallers.df_new$counts[svcallers.df_new$svtype=="DUP"])
mean(svcallers.df_new$counts[svcallers.df_new$svtype=="DEL"])

### Counts of SVs per Sample
sv_count_sample<-ggplot(svcallers.df, aes(x=sample, y=counts, fill=condition)) + 
  geom_bar(stat="identity")+
  facet_grid(~svtype, scales="free_x")+
  coord_flip()+ #scale_x_discrete(limits=rev)+
  scale_fill_manual(breaks = c("control", "evolved", "orig-control", "orig-evolved"), 
                    values=c("blue", "red","dark blue","dark red"))


ggplot(svcallers.df_new, aes(x=isolate, y=counts, fill=condition)) + 
  geom_boxplot()+
  facet_wrap(~svtype, scales="free_y")+
  stat_compare_means(comparisons = my_comparisons, 
                     method="wilcox.test")+
  scale_fill_manual(breaks = c("control", "evolved", "orig-control", "orig-evolved"), 
                    values=c("blue", "red","dark blue","dark red"))

####################################################################################
### SNP accumulation (##also, check if descendants have same mutation as ancestor##)
####################################################################################
lineages<-c("WTC","PKSC",
            "WT15.1.2","WT15.2.1","WT15.2.2","WT15.4.1","WT15.4.2",
            "PKS15.1.1","PKS15.2.1","PKS15.2.2","PKS15.3.2","PKS15.4.2")
### Find new SNPs after 8000Gy 
AncDesc<-SNPs_df0
for (i in lineages) {
  ancestor <-grep(paste0("orig-",i), names(SNPs_df0))
  descendants <- grep(paste0("^",i), names(SNPs_df0))
  ##df<-SNPs_df0[apply(SNPs_df0[ancestor] == 1, 1, any), ]
  ##df2<-df[apply(df[descendants] == 0, 1, any), ]
  ##print(paste(i,length(df2$POS)))
  df<-SNPs_df0
  df[,3:94]<-lapply(df[,3:94], function(x) {as.numeric(as.character(x))})
  df2<-df
  df2<-apply(df[ancestor],2,
             function(x) df[descendants] - x )
  AncDesc[,descendants]<-as.data.frame(df2)
}

SNPaccum_df<-as.data.frame(colSums(AncDesc[grep("^WT|^PKS",colnames(AncDesc))],na.rm = TRUE))
SNPaccum_df$Sample<-row.names(SNPaccum_df)
colnames(SNPaccum_df)[1]<-"accumulated_SNPs"
SNPaccum_df$condition<-gsub(".*C.*","control",SNPaccum_df$Sample)
SNPaccum_df$condition<-gsub(".*15.*","evolved",SNPaccum_df$condition)

SNPaccum_df$Group<-gsub("^WT15.*","WT15",SNPaccum_df$Sample)
SNPaccum_df$Group<-gsub("^PKS15.*","PKS15",SNPaccum_df$Group)
SNPaccum_df$Group<-gsub("^WTC.*","WTC",SNPaccum_df$Group)
SNPaccum_df$Group<-gsub("^PKSC.*","PKSC",SNPaccum_df$Group)
SNPaccum_df$Group<-factor(SNPaccum_df$Group, levels = c("WTC", "WT15","PKSC","PKS15"))
SNPaccum_df$lineage<-SNPaccum_df$Sample
SNPaccum_df$lineage[SNPaccum_df$Group=="WT15"|SNPaccum_df$Group=="PKS15"]<-gsub('\\..$',"",SNPaccum_df$lineage[SNPaccum_df$Group=="WT15"|SNPaccum_df$Group=="PKS15"])

SNPaccum_df_lineage<-aggregate(SNPaccum_df$accumulated_SNPs, by=list(SNPaccum_df$lineage), FUN=mean)
SNPaccum_df_lineage$Group<-gsub("^WT15.*","WT15",SNPaccum_df_lineage$Group.1)
SNPaccum_df_lineage$Group<-gsub("^PKS15.*","PKS15",SNPaccum_df_lineage$Group)
SNPaccum_df_lineage$Group<-gsub("^WTC.*","WTC",SNPaccum_df_lineage$Group)
SNPaccum_df_lineage$Group<-gsub("^PKSC.*","PKSC",SNPaccum_df_lineage$Group)
SNPaccum_df_lineage$condition<-gsub(".*C.*","control",SNPaccum_df_lineage$Group)
SNPaccum_df_lineage$condition<-gsub(".*15.*","evolved",SNPaccum_df_lineage$condition)
SNPaccum_df_lineage$type<-"SNP"
colnames(SNPaccum_df_lineage)[1:2]<-c("lineage","accumulated_SNPs")

kw.SNPaccum_populations<-SNPaccum_df_lineage%>% kruskal_test(accumulated_SNPs ~ Group)
#.y.                  n statistic    df     p method        
#accumulated_SNPs    24      1.39     3 0.709 Kruskal-Wallis

SNPaccum_df$lineage<-gsub("C.*","C",SNPaccum_df$lineage)
SNPaccum_df$type<-"SNP"

SNPaccum_8000Gy_lineage<-ggplot(SNPaccum_df, aes(x=lineage, y=accumulated_SNPs, fill=condition)) +
  geom_boxplot()+ xlab("Lineage")+ylab("No. of SNPs after 8000 Gy")+
  facet_wrap(~type)+
  scale_fill_manual(breaks = c("control", "evolved"), 
                    values=c("blue", "red"))

kw.SNPaccum_lineage<-SNPaccum_df%>% kruskal_test(accumulated_SNPs ~ lineage)
#y.                  n statistic    df     p method        
#accumulated_SNPs    78      13.1     8 0.109 Kruskal-Wallis

####################################################################################
############################ SVs accumulation ######################################
####################################################################################

AncDesc_SV<-svcallers
for (i in lineages) {
  ancestor <-grep(paste0("orig-",i), names(svcallers))
  descendants <- grep(paste0("^",i), names(svcallers))
  ##df<-svcallers[apply(svcallers[ancestor] == 1, 1, any), ]
  ##df2<-df[apply(df[descendants] == 0, 1, any), ]
  ##print(paste(i,length(df2$POS)))
  df<-svcallers
  df[,3:94]<-lapply(df[,3:94], function(x) {as.numeric(as.character(x))})
  df2<-df
  df2<-apply(df[ancestor],2,
             function(x) df[descendants] - x )
  df2<-as.data.frame(df2)
  df2[df2<0]<-0
  AncDesc_SV[,descendants]<-as.data.frame(df2)
  
}


SVaccum_df<-data.frame()
for (svtype in unique(AncDesc_SV$PE)) {
  m<-AncDesc_SV[AncDesc_SV$PE==svtype,]
  AncDesc_SV_stats<-data.frame(colSums(m[6:96]))
  AncDesc_SV_stats$sample<-rownames(AncDesc_SV_stats)
  colnames(AncDesc_SV_stats)<-c("counts","sample")
  AncDesc_SV_stats$Group<-gsub("15.*","15",AncDesc_SV_stats$sample)
  AncDesc_SV_stats$Group<-gsub("C.*","C",AncDesc_SV_stats$Group)
  AncDesc_SV_stats$svtype<-svtype
  SVaccum_df<-rbind(SVaccum_df,data.frame(AncDesc_SV_stats))
}

SVaccum_df<-SVaccum_df[which(SVaccum_df$sample!="PKS15.1.1.1"),]##Remove PKS15.1.1--Abandoned because not enough SNPs; might affect SV results
SVaccum_df$condition<-SVaccum_df$Group
SVaccum_df$condition[SVaccum_df$condition=="PKSC" | SVaccum_df$condition=="WTC"]<-"control"
SVaccum_df$condition[SVaccum_df$condition=="PKS15" | SVaccum_df$condition=="WT15"]<-"evolved"
SVaccum_df$condition[SVaccum_df$condition=="orig-PKSC" | SVaccum_df$condition=="orig-WTC"]<-"orig-control"
SVaccum_df$condition[SVaccum_df$condition=="orig-PKS15" | SVaccum_df$condition=="orig-WT15"]<-"orig-evolved"
SVaccum_df$condition <- factor(SVaccum_df$condition, levels = c("control","evolved","orig-control","orig-evolved"))
SVaccum_df<-SVaccum_df[which(SVaccum_df$condition=="evolved"|SVaccum_df$condition=="control"),]
SVaccum_df$lineage<-gsub('\\..$',"",SVaccum_df$sample)



colnames(SVaccum_df)<-c("accumulated_SVs","sample","Group","svtype","condition","lineage")

SVaccum_df_lineage<-aggregate(SVaccum_df$accumulated_SVs, by=list(SVaccum_df$lineage,SVaccum_df$svtype), FUN=mean)
colnames(SVaccum_df_lineage)<-c("lineage","svtype","accumulated_SVs")
SVaccum_df_lineage$Group<-gsub('\\..*$',"",SVaccum_df_lineage$lineage)
SVaccum_df_lineage$Group<-gsub('C.',"C",SVaccum_df_lineage$Group)
SVaccum_df_lineage$condition<-SVaccum_df_lineage$Group

SVaccum_df_lineage$condition<-gsub("PKSC.*|WTC.*","control",SVaccum_df_lineage$condition)
SVaccum_df_lineage$condition[SVaccum_df_lineage$condition=="PKS15" | SVaccum_df_lineage$condition=="WT15"]<-"evolved"
SVaccum_df_lineage$condition[SVaccum_df_lineage$condition=="orig-PKSC" | SVaccum_df_lineage$condition=="orig-WTC"]<-"orig-control"
SVaccum_df_lineage$condition[SVaccum_df_lineage$condition=="orig-PKS15" | SVaccum_df_lineage$condition=="orig-WT15"]<-"orig-evolved"
SVaccum_df_lineage<-data.frame(SVaccum_df_lineage)
SVaccum_df_lineage<-SVaccum_df_lineage[grep("^WT|^PKS",SVaccum_df_lineage$Group),]


SVaccum_df$lineage<-gsub("C.*","C",SVaccum_df$lineage)
SVaccum_df<-SVaccum_df[grep("^WT|^PKS",SVaccum_df$sample),]
kw.SVaccum_lineage<-SVaccum_df %>% 
  group_by(svtype) %>%
  kruskal_test(accumulated_SVs ~ lineage)
dunn.SVaccum_lineage<-SVaccum_df %>% 
  group_by(svtype) %>%
  dunn_test(accumulated_SVs ~ lineage)

dunn.SVaccum_lineage<-as.data.frame(dunn.SVaccum_lineage)
dunn.SVaccum_lineage<-dunn.SVaccum_lineage[(grepl("WTC",dunn.SVaccum_lineage$group1)&grepl("PKSC",dunn.SVaccum_lineage$group2))|
                                             (grepl("PKSC",dunn.SVaccum_lineage$group1)&grepl("WTC",dunn.SVaccum_lineage$group2))|
                                             (grepl("WTC",dunn.SVaccum_lineage$group1)&grepl("WT15",dunn.SVaccum_lineage$group2))|
                                             (grepl("WT15",dunn.SVaccum_lineage$group1)&grepl("WTC",dunn.SVaccum_lineage$group2))|
                                             (grepl("PKSC",dunn.SVaccum_lineage$group1)&grepl("PKS15",dunn.SVaccum_lineage$group2))|
                                             (grepl("PKS15",dunn.SVaccum_lineage$group1)&grepl("PKSC",dunn.SVaccum_lineage$group2))|
                                             (grepl("PKS15",dunn.SVaccum_lineage$group1)&grepl("PKS15",dunn.SVaccum_lineage$group2))|
                                             (grepl("WT15",dunn.SVaccum_lineage$group1)&grepl("WT15",dunn.SVaccum_lineage$group2)),]
dunn.SVaccum_lineage<-dunn.SVaccum_lineage[dunn.SVaccum_lineage$p.adj<=0.05,]
dunn.SVaccum_lineage$condition<-NA

SVaccum_8000Gy_lineage<-ggplot(SVaccum_df, aes(x=lineage, y=accumulated_SVs, fill=condition)) +
  geom_boxplot()+ylab("No. of SVs after 8000 Gy")+xlab("Lineages")+
  facet_wrap(~svtype, scales="free_y",ncol=2)+
  scale_fill_manual(breaks = c("control", "evolved"), 
                    values=c("blue", "red"))+
  stat_pvalue_manual(dunn.SVaccum_lineage, step.group.by="svtype",
                     y.position = c(2.2,2.35,2.50,2.65,
                                    5,5.5,6,6.5,7,
                                    35,37,39,41))


SVaccum_df[c("accumulated_SVs","sample","Group","svtype","condition")]
colnames(svcallers.df_orig)<-c("accumulated_SVs","sample","Group","svtype","condition")

sv_8000Gy_sample<-ggplot(rbind(SVaccum_df[c("accumulated_SVs","sample","Group","svtype","condition")],svcallers.df_orig), 
                         aes(x=sample, y=accumulated_SVs, fill=condition)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~svtype, scales="free_x",ncol=4)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background=element_blank())+
  scale_fill_manual(breaks = c("control", "evolved",
                               "orig-control","orig-evolved"), 
                    values=c("blue", "red","gray","black"))+
  coord_flip() +ylab("SNP counts") + scale_x_discrete(limits=rev)


####################################################################################
###################################### SNPeff ######################################
####################################################################################
### SNPeff statistics
mu_type<-unique(SNPeff$Type)

SNPeff_sample<-list()
for (i in mu_type) {
  temp<-SNPeff[c(2,6:97)]
  temp<- lapply(subset(temp, Type==i,select=-Type),table)
  temp <-plyr::ldply(temp, rbind)
  temp["1"][is.na(temp["1"])]<-0 
  temp$genotyped<-rowSums(temp[c("0","1")],na.rm=TRUE)
  temp$alt<-unlist(temp["1"])
  temp$freq_0<-unlist(temp["0"]/temp$genotyped)
  temp$freq_1<-unlist(temp["1"]/temp$genotyped)
  temp$freq_alt<-unlist((temp["1"])/temp$genotyped)
  temp$Type<-i
  SNPeff_sample[[i]]<-temp
}

SNPeff_sample=plyr::ldply(SNPeff_sample, rbind)
colnames(SNPeff_sample)[1:5]<-c("Sample","mis","count_0","count_1","genotyped")
SNPeff_sample$condition<-gsub("orig-.*C.*","orig-control",SNPeff_sample$Sample)
SNPeff_sample$condition<-gsub("orig-.*15.*","orig-evolved",SNPeff_sample$condition)
SNPeff_sample$condition<-gsub(".*C.*","control",SNPeff_sample$condition)
SNPeff_sample$condition<-gsub(".*15.*","evolved",SNPeff_sample$condition)

SNPeff_sample$Group<-gsub("^WT15.*","WT15",SNPeff_sample$Sample)
SNPeff_sample$Group<-gsub("^PKS15.*","PKS15",SNPeff_sample$Group)
SNPeff_sample$Group<-gsub("^WTC.*","WTC",SNPeff_sample$Group)
SNPeff_sample$Group<-gsub("^PKSC.*","PKSC",SNPeff_sample$Group)
SNPeff_sample$Group<-gsub("orig-WT15.*","orig-WT15",SNPeff_sample$Group)
SNPeff_sample$Group<-gsub("orig-PKS15.*","orig-PKS15",SNPeff_sample$Group)
SNPeff_sample$Group<-gsub("orig-WTC.*","orig-WTC",SNPeff_sample$Group)
SNPeff_sample$Group<-gsub("orig-PKSC.*","orig-PKSC",SNPeff_sample$Group)
SNPeff_sample$Group<-factor(SNPeff_sample$Group, levels = c("WTC", "WT15","PKSC","PKS15",
                                                            "orig-WTC", "orig-WT15","orig-PKSC","orig-PKS15"))
SNPeff_sample$Isolate<-gsub("C.$","C",SNPeff_sample$Sample)
SNPeff_sample<-SNPeff_sample[which(SNPeff_sample$condition=="evolved"|SNPeff_sample$condition=="control"),]
SNPeff_sample$Isolate<-gsub("\\..$","",SNPeff_sample$Isolate)

SNPeff_sample$accumulated_snp_generation<-SNPeff_sample$count_1


kw.SNPeff_lineage<-SNPeff_sample %>% 
  group_by(Type) %>%
  kruskal_test(alt ~ Isolate)
kw.SNPeff_lineage[which(kw.SNPeff_lineage$p<=0.05),]

dunn.SNPeff_lineage<-SNPeff_sample %>% 
  group_by(Type) %>%
  dunn_test(alt ~ Isolate)

dunn.SNPeff_lineage<-as.data.frame(dunn.SNPeff_lineage)
dunn.SNPeff_lineage<-dunn.SNPeff_lineage[(grepl("WTC",dunn.SNPeff_lineage$group1)&grepl("PKSC",dunn.SNPeff_lineage$group2))|
                                           (grepl("PKSC",dunn.SNPeff_lineage$group1)&grepl("WTC",dunn.SNPeff_lineage$group2))|
                                           (grepl("WTC",dunn.SNPeff_lineage$group1)&grepl("WT15",dunn.SNPeff_lineage$group2))|
                                           (grepl("WT15",dunn.SNPeff_lineage$group1)&grepl("WTC",dunn.SNPeff_lineage$group2))|
                                           (grepl("PKSC",dunn.SNPeff_lineage$group1)&grepl("PKS15",dunn.SNPeff_lineage$group2))|
                                           (grepl("PKS15",dunn.SNPeff_lineage$group1)&grepl("PKSC",dunn.SNPeff_lineage$group2))|
                                           (grepl("PKS15",dunn.SNPeff_lineage$group1)&grepl("PKS15",dunn.SNPeff_lineage$group2))|
                                           (grepl("WT15",dunn.SNPeff_lineage$group1)&grepl("WT15",dunn.SNPeff_lineage$group2)),]
dunn.SNPeff_lineage<-dunn.SNPeff_lineage[(!grepl("^orig",dunn.SNPeff_lineage$group1)),]
dunn.SNPeff_lineage<-dunn.SNPeff_lineage[(!grepl("^orig",dunn.SNPeff_lineage$group2)),]
dunn.SNPeff_lineage<-dunn.SNPeff_lineage[dunn.SNPeff_lineage$p.adj<=0.05,]
dunn.SNPeff_lineage$condition<-NA

SNPeff_isolate_plot<-ggplot(SNPeff_sample, aes(x=Isolate, y=alt, fill=condition)) +
  geom_boxplot()+ facet_wrap(~Type, scales="free_x")+
  scale_fill_manual(breaks = c("control", "evolved",
                               "orig-control","orig-evolved"), 
                    values=c("blue", "red","dark blue","dark red"))+
  coord_flip()+ scale_x_discrete(limits=rev)+
  stat_pvalue_manual(dunn.SNPeff_lineage,y.position=c(400))


kw.SNPeff_populations<-SNPeff_sample %>% 
  group_by(Type) %>%
  kruskal_test(alt ~ Group)
kw.SNPeff_populations[kw.SNPeff_populations$p<=0.05,]$Type

dunn.SNPeff_populations<-SNPeff_sample[SNPeff_sample$Type %in% c(kw.SNPeff_populations[kw.SNPeff_populations$p<=0.05,]$Type),] %>% 
  group_by(Type) %>%
  dunn_test(alt ~ Group)

dunn.SNPeff_populations<-as.data.frame(dunn.SNPeff_populations)
dunn.SNPeff_populations<-dunn.SNPeff_populations[(grepl("WTC",dunn.SNPeff_populations$group1)&grepl("PKSC",dunn.SNPeff_populations$group2))|
                                                   (grepl("PKSC",dunn.SNPeff_populations$group1)&grepl("WTC",dunn.SNPeff_populations$group2))|
                                                   (grepl("WTC",dunn.SNPeff_populations$group1)&grepl("WT15",dunn.SNPeff_populations$group2))|
                                                   (grepl("WT15",dunn.SNPeff_populations$group1)&grepl("WTC",dunn.SNPeff_populations$group2))|
                                                   (grepl("PKSC",dunn.SNPeff_populations$group1)&grepl("PKS15",dunn.SNPeff_populations$group2))|
                                                   (grepl("PKS15",dunn.SNPeff_populations$group1)&grepl("PKSC",dunn.SNPeff_populations$group2)),]
dunn.SNPeff_populations<-dunn.SNPeff_populations[(!grepl("orig",dunn.SNPeff_populations$group1)),]
dunn.SNPeff_populations<-dunn.SNPeff_populations[(!grepl("orig",dunn.SNPeff_populations$group2)),]
dunn.SNPeff_populations<-dunn.SNPeff_populations[dunn.SNPeff_populations$p.adj<=0.05,]
dunn.SNPeff_populations$condition<-NA


pdf("SNPeff_sample_isolate.pdf",width=25,height=50)
#SNPeff_sample_plot
SNPeff_isolate_plot
graphics.off()


AncDesc_SNPeff<-merge(SNPeff[,c("CHROM","POS","gene","Type")],AncDesc, by=c("CHROM","POS"),all=FALSE)
AncDesc_SNPeff[5:96]<-as.data.frame(sapply(AncDesc_SNPeff[5:96],as.numeric))

SNPeffaccum_df<-data.frame()
for (Type in unique(AncDesc_SNPeff$Type)) {
  m<-AncDesc_SNPeff[AncDesc_SNPeff$Type==Type,]
  AncDesc_SNPeff_stats<-data.frame(colSums(m[5:96],na.rm =TRUE))
  AncDesc_SNPeff_stats$sample<-rownames(AncDesc_SNPeff_stats)
  colnames(AncDesc_SNPeff_stats)<-c("counts","sample")
  AncDesc_SNPeff_stats$Group<-gsub("15.*","15",AncDesc_SNPeff_stats$sample)
  AncDesc_SNPeff_stats$Group<-gsub("C.*","C",AncDesc_SNPeff_stats$Group)
  AncDesc_SNPeff_stats$Type<-Type
  SNPeffaccum_df<-rbind(SNPeffaccum_df,data.frame(AncDesc_SNPeff_stats))
}

SNPeffaccum_df<-SNPeffaccum_df[which(SNPeffaccum_df$sample!="PKS15.1.1.1"),]##Remove PKS15.1.1--Abandoned because not enough SNPs; might affect SNPeff results
SNPeffaccum_df$condition<-SNPeffaccum_df$Group
SNPeffaccum_df$condition[SNPeffaccum_df$condition=="PKSC" | SNPeffaccum_df$condition=="WTC"]<-"control"
SNPeffaccum_df$condition[SNPeffaccum_df$condition=="PKS15" | SNPeffaccum_df$condition=="WT15"]<-"evolved"
SNPeffaccum_df$condition[SNPeffaccum_df$condition=="orig-PKSC" | SNPeffaccum_df$condition=="orig-WTC"]<-"orig-control"
SNPeffaccum_df$condition[SNPeffaccum_df$condition=="orig-PKS15" | SNPeffaccum_df$condition=="orig-WT15"]<-"orig-evolved"
SNPeffaccum_df$condition <- factor(SNPeffaccum_df$condition, levels = c("control","evolved","orig-control","orig-evolved"))
SNPeffaccum_df<-SNPeffaccum_df[which(SNPeffaccum_df$condition=="evolved"|SNPeffaccum_df$condition=="control"),]

SNPeffaccum_df$lineage<-gsub('\\..$',"",SNPeffaccum_df$sample)
colnames(SNPeffaccum_df)<-c("accumulated_SNPeffs","sample","Group","SNPefftype","condition","lineage")

SNPeffaccum_df_lineage<-aggregate(SNPeffaccum_df$accumulated_SNPeffs, by=list(SNPeffaccum_df$lineage,SNPeffaccum_df$SNPefftype), FUN=mean)
colnames(SNPeffaccum_df_lineage)<-c("lineage","SNPefftype","accumulated_SNPeffs")

SNPeffaccum_df_lineage[which(SNPeffaccum_df_lineage$condition=="evoled"|SNPeffaccum_df_lineage$condition=="control"),]
SNPeffaccum_df_lineage$Group<-gsub('\\..*$',"",SNPeffaccum_df_lineage$lineage)
SNPeffaccum_df_lineage$Group<-gsub('C.',"C",SNPeffaccum_df_lineage$Group)
SNPeffaccum_df_lineage$condition<-SNPeffaccum_df_lineage$Group

SNPeffaccum_df_lineage$condition<-gsub("PKSC.*|WTC.*","control",SNPeffaccum_df_lineage$condition)
SNPeffaccum_df_lineage$condition[SNPeffaccum_df_lineage$condition=="PKS15" | SNPeffaccum_df_lineage$condition=="WT15"]<-"evolved"
SNPeffaccum_df_lineage$condition[SNPeffaccum_df_lineage$condition=="orig-PKSC" | SNPeffaccum_df_lineage$condition=="orig-WTC"]<-"orig-control"
SNPeffaccum_df_lineage$condition[SNPeffaccum_df_lineage$condition=="orig-PKS15" | SNPeffaccum_df_lineage$condition=="orig-WT15"]<-"orig-evolved"
SNPeffaccum_df_lineage<-data.frame(SNPeffaccum_df_lineage)
SNPeffaccum_df_lineage<-SNPeffaccum_df_lineage[grep("^WT|^PKS",SNPeffaccum_df_lineage$Group),]
SNPeffaccum_df$lineage<-gsub("C.*","C",SNPeffaccum_df$lineage)



SNPeffaccum_df$lineage<-gsub("C.*","C",SNPeffaccum_df$lineage)
SNPeffaccum_df<-SNPeffaccum_df[grep("^WT|^PKS",SNPeffaccum_df$sample),]
kw.SNPeffaccum_lineage<-SNPeffaccum_df %>% 
  group_by(SNPefftype) %>%
  kruskal_test(accumulated_SNPeffs ~ lineage)
dunn.SNPeffaccum_lineage<-SNPeffaccum_df[which(SNPeffaccum_df$SNPefftype=="SPLICE_SITE_ACCEPTOR+INTRON"),] %>% 
  group_by(SNPefftype) %>%
  dunn_test(accumulated_SNPeffs ~ lineage)

dunn.SNPeffaccum_lineage<-as.data.frame(dunn.SNPeffaccum_lineage)
dunn.SNPeffaccum_lineage<-dunn.SNPeffaccum_lineage[(grepl("WTC",dunn.SNPeffaccum_lineage$group1)&grepl("PKSC",dunn.SNPeffaccum_lineage$group2))|
                                                     (grepl("PKSC",dunn.SNPeffaccum_lineage$group1)&grepl("WTC",dunn.SNPeffaccum_lineage$group2))|
                                                     (grepl("WTC",dunn.SNPeffaccum_lineage$group1)&grepl("WT15",dunn.SNPeffaccum_lineage$group2))|
                                                     (grepl("WT15",dunn.SNPeffaccum_lineage$group1)&grepl("WTC",dunn.SNPeffaccum_lineage$group2))|
                                                     (grepl("PKSC",dunn.SNPeffaccum_lineage$group1)&grepl("PKS15",dunn.SNPeffaccum_lineage$group2))|
                                                     (grepl("PKS15",dunn.SNPeffaccum_lineage$group1)&grepl("PKSC",dunn.SNPeffaccum_lineage$group2)),]
dunn.SNPeffaccum_lineage<-dunn.SNPeffaccum_lineage[dunn.SNPeffaccum_lineage$p.adj<=0.05,]
dunn.SNPeffaccum_lineage$condition<-NA

SNPeffaccum_8000Gy_lineage<-ggplot(SNPeffaccum_df, aes(x=lineage, y=accumulated_SNPeffs, fill=condition)) +
  geom_boxplot()+ylab("No. of SNPeffs after 8000 Gy")+xlab("Lineages")+
  facet_wrap(~SNPefftype, scales="free_y",ncol=2)+
  scale_fill_manual(breaks = c("control", "evolved"), 
                    values=c("blue", "red"))+
  stat_pvalue_manual(dunn.SNPeffaccum_lineage, step.group.by="SNPefftype",
                     y.position = c(2.2,2.4,2.6))

##########################################################################################
#### Parallel SNPs: mutations per Loci ###################################################
##########################################################################################
SNPs_loci<-t(SNPs_df0)
colnames(SNPs_loci)<-paste0(SNPs_loci[1,],":",SNPs_loci[2,])
SNPs_loci<-SNPs_loci[3:94,]
SNPs_loci<-as.data.frame(cbind(SNPs_loci))
#SNPs_loci_summary<-t(plyr::ldply(lapply(SNPs_loci,table),rbind))[c(2,3,5),]
SNPs_loci_summary<-(plyr::ldply(lapply(SNPs_loci,table),rbind))
SNPs_loci_summary[3,][is.na(SNPs_loci_summary[3,])]<-0
SNPs_loci_summary<-t(SNPs_loci_summary)
colnames(SNPs_loci_summary)<-colnames(SNPs_loci)
SNPs_loci<-t(rbind(SNPs_loci,SNPs_loci_summary))

###count frequency of SNPs per loci
colnames(SNPs_loci)[94:95]<-c("count_0","count_1")
SNPs_loci<-data.frame(SNPs_loci)
SNPs_loci$count_0<-as.numeric(as.character(SNPs_loci$count_0))
SNPs_loci$count_1<-as.numeric(as.character(SNPs_loci$count_1))
###SNPs_loci$count_2<-as.numeric(as.character(SNPs_loci$count_2))

SNPs_loci$genotyped<-SNPs_loci$count_0+SNPs_loci$count_1##+SNPs_loci$count_2
SNPs_loci$alt<-unlist(SNPs_loci$count_1)##+SNPs_loci$count_2)
SNPs_loci$freq_0<-unlist(SNPs_loci$count_0/SNPs_loci$genotyped)
SNPs_loci$freq_1<-unlist(SNPs_loci$count_1/SNPs_loci$genotyped)
##SNPs_loci$freq_2<-unlist(SNPs_loci$count_2/SNPs_loci$genotyped)
SNPs_loci$freq_alt<-unlist((SNPs_loci$count_1)/SNPs_loci$genotyped)


### Identify parallel SNPs at the locus level
Groups<-c("WTC","PKSC",unique(gsub("\\..$","",colnames(SNPs_loci)[18:79])))
##group isolate "reps" to identify SNPs that occurred in separate samples
for (i in Groups) {
  DF<-SNPs_loci[,grep(pattern=paste0(i,".*"), x=colnames(SNPs_loci))]
  ###check if there is a consensus between samples/replicates of an Group
  SNPs_loci[paste0(i,"_same")]<-transform(DF, same = apply(DF, 1, function(x) length(unique(x)) == 1))["same"]
  ###check if there are potential parallel Groups (get alternative allele "1", "2" is unlikely)
  SNPs_loci[paste0(i,"_max")]<-transform(apply(DF, 1, function(x) max(unique(x))))
}

max_allele<-SNPs_loci[grep("_max",colnames(SNPs_loci))]
###max_allele[max_allele$WTC_max>0 & max_allele$PKSC_max>0,] ### no parallelism between WTC and PKSC almost all isolates have the same mutation
max_allele<-cbind(max_allele,plyr::ldply(lapply(data.frame(t(max_allele)),table),rbind))
max_allele<-transform(max_allele, same = apply(DF, 1, function(x) length(unique(x)) == 1))
max_allele<-max_allele[- grep("1", max_allele$WTC_max),]
max_allele<-max_allele[- grep("1", max_allele$PKSC_max),]
max_allele$tot<-rowSums(as.data.frame(lapply(max_allele[,1:9],function(x){as.numeric(x)})))
#scaffold_4: 293156
#scaffold_4: 293183

write.csv(SNPs_loci,"SNPs_loci.csv")

### SNPeff statistics
SNPeff$tot_snps<-rowSums(SNPeff[,101:120])
write.csv(SNPeff,"SNPs_loci.snpeff.csv")

### Identify Parallelism at Gene level: SNPs occur at different sites of same gene 
##for each row, identify samples with GT="1", list their names but remove the replicate notation
SNPeff$Isolate<-apply(SNPeff[,grep("^WT|^PKS",colnames(SNPeff))]=='1', 1, FUN= function(x) toString(unique(gsub("\\..$","",names(x)[x]))))
SNPeff$Anc<-apply(SNPeff[,grep("^orig",colnames(SNPeff))]=='1', 1, FUN= function(x) toString(unique(names(x)[x])))

##remove loci where control has GT="1"
SNPeff$Isolate<-gsub("C",NA,SNPeff$Isolate)
## Has Genes; Isolates: PKS15.1.1, PKS15.3.2
parallel_loci<-SNPeff[grep(",",SNPeff$Isolate),]

####Function searches for Parallelism at Gene level but also occurrence in Ancestral sequences
Parallel_Genes<-function(variants,filename){
  SNPeff_coding<-SNPeff[variants,]
  parallel_genes<-table(SNPeff_coding$gene)
  ##genes with multiple SNPs (nonsyn, syn, start/stop/frameshift) from different samples
  parallel_genes<-names(parallel_genes[parallel_genes>1])
  SNPeff_coding<-SNPeff_coding[SNPeff_coding$gene %in% parallel_genes,]
  
  ##get list of Anc names:compress list of genes into comma delimited list
  SNPeff_coding_orig<-SNPeff_coding[grep("orig",SNPeff_coding$Anc),]
  SNPeff_coding_orig<-unique(SNPeff_coding_orig[,c("gene","Anc")])
  SNPeff_coding_orig_tmp<-reshape(SNPeff_coding[,c("gene","Anc","Anc")],
                                  idvar=c("gene"),timevar = "Anc",direction="wide")
  SNPeff_coding_orig_tmp$Anc<-apply(SNPeff_coding_orig_tmp[,grep("orig",colnames(SNPeff_coding_orig_tmp))],1,FUN=function(x) toString(na.omit(x)))
  SNPeff_coding_orig_tmp<-SNPeff_coding_orig_tmp[grep(",",SNPeff_coding_orig_tmp$Anc),][c("gene","Anc")]
  
  ##get list of Isolate names:compress list of genes into comma delimited list
  SNPeff_coding_tmp<-reshape(SNPeff_coding[,c("gene","Isolate","Isolate")],
                             idvar=c("gene"),timevar = "Isolate",direction="wide")
  SNPeff_coding_tmp$Isolate<-apply(SNPeff_coding_tmp[,grep("Isolate",colnames(SNPeff_coding_tmp))],1,FUN=function(x) toString(na.omit(x)))
  SNPeff_coding_tmp<-SNPeff_coding_tmp[grep(", $",SNPeff_coding_tmp$Isolate,invert=TRUE),]
  SNPeff_coding_tmp<-SNPeff_coding_tmp[grep(",",SNPeff_coding_tmp$Isolate),][c("gene","Isolate")]
  ##get list of Types
  SNPeff_coding_tmp2<-reshape(SNPeff_coding[,c("gene","Type","Type")],
                              idvar=c("gene"),timevar = "Type",direction="wide")
  SNPeff_coding_tmp2$Type<-apply(SNPeff_coding_tmp2[,grep("Type",colnames(SNPeff_coding_tmp2))],1,FUN=function(x) toString(na.omit(x)))
  SNPeff_coding_tmp2<-SNPeff_coding_tmp2[,c("gene","Type")]
  SNPeff_coding<-Reduce(function(...) merge(..., by='gene'), 
                        list(SNPeff_coding_tmp,SNPeff_coding_tmp2,SNPeff_coding_orig_tmp,annot))
  SNPeff_coding$events<-str_count(SNPeff_coding$Isolate, ",")+1
  SNPeff_coding<-merge(genenames,SNPeff_coding,by="gene")
  write.csv(SNPeff_coding,filename)
  
  table<-SNPeff_coding
  ego<-enricher(gene = table$gene, 
                universe = annot$gene,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                TERM2GENE = term2gene_go[,c("goName","gene")])
  a1<-dotplot(ego,orderBy="GeneRatio")+ggtitle("activated")
  a1<-tryCatch({print(a1)}, error = function(e){
    a1<-NA
    return(a1)})
  ekegg <- enricher(gene = table$gene,
                    universe = annot$gene,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    TERM2GENE = term2gene_kegg[,c("pathway","gene")])
  b1<-dotplot(ekegg,orderBy="GeneRatio")+ggtitle("activated")
  b1<-tryCatch({print(b1)}, error = function(e){
    b1<-NA
    return(b1)})
  pdf(paste0(filename,".go.pdf"))
  print(a1)
  graphics.off()
  
  pdf(paste0(filename,".kegg.pdf"))
  print(b1)
  graphics.off()
}

Parallel_Genes(variants = grep("STREAM|UTR|INTR|SPLICE|SYN",SNPeff$Type),filename="parallel_genes.aroundcoding_regions.csv") ### NO GO or KEGG
Parallel_Genes(variants = -grep("STREAM|UTR|SPLICE|INTR",SNPeff$Type),filename="parallel_genes.coding_regions.csv") ### NO GO or KEGG
Parallel_Genes(variants = grep("START|STOP|FRAME",SNPeff$Type),filename="parallel_genes.start_stop_frame.csv") ### *YES GO but NO KEGG*
Parallel_Genes(variants = grep("SYN",SNPeff$Type),filename="parallel_genes.nonsynonymous_synonymous.csv") ### NO GO or KEGG

Parallel_Genes(variants = grep("SYN|START|STOP|FRAME",SNPeff$Type),filename="parallel_genes.nonsyn_syn_start_stop_frame.csv") 
parallel_genes.coding<-read.csv("parallel_genes.nonsyn_syn_start_stop_frame.csv")
parallel_genes.aroundcoding<-read.csv("parallel_genes.aroundcoding_regions.csv")

### Top 10 Genes w/annotation with highest associated SNPs
Top10<-list()
variants<-colnames(SNPeff)[grep("^variant.*|tot_.*",colnames(SNPeff))]
for (i in variants) {
  Top10[[i]]<-unique(SNPeff[order(-SNPeff[,i]),c("gene","genename","goName", "pathway","kogdefline","iprDesc",i)])[1:10,]
  Top10[[i]]<-subset(Top10[[i]],Top10[[i]][,i]>=2)
}
write.csv(Top10[["variants_impact_HIGH"]],"variants_impact_HIGH.csv")
View(Top10[["variants_impact_MODERATE"]])
View(Top10[["variants_impact_MODIFIER"]])

write.csv(Top10[["variants_effect_NON_SYNONYMOUS_CODING"]],"variants_effect_NON_SYNONYMOUS_CODING.csv")
View(Top10[["variants_effect_SYNONYMOUS_CODING"]])

View(Top10[["variants_effect_FRAME_SHIFT"]])

View(Top10[["variants_effect_START_GAINED"]])
#View(Top10[["variants_effect_START_LOST"]])
View(Top10[["variants_effect_STOP_GAINED"]])
#View(Top10[["variants_effect_STOP_LOST"]])

View(Top10[["variants_effect_UTR_3_PRIME"]])
View(Top10[["variants_effect_UTR_5_PRIME"]])

#####################################
##### Identify Parallel SVs #########
#####################################
svcallers_parallel<-svcallers
svcallers_parallel$counts<-rowSums(svcallers_parallel[6:96]) ### count number of mutations
svcallers_parallel<-svcallers_parallel[svcallers_parallel$counts>1,] ### select mutations that occur more than once

svcallers_parallel$isolates<-0
###for each locus, list isolate name for every "1" genotype
for (i in rownames(svcallers_parallel)) {
  svcallers_parallel[i,"isolates"]<-paste0(colnames(svcallers_parallel[i,6:96][,svcallers_parallel[i,6:96]==1]),collapse = ",")
}
write.csv(svcallers_parallel, "SV_parallel.csv") ### writes a file that has loci with multiple SV mutations, must verify if truly parallel

##WTC and PKSC
WTCPKSC_in<-svcallers_parallel[grepl("PKSC.*WTC|WTC.*PKSC",svcallers_parallel$isolates),] ##has WTC and PKSC
#WTCPKSC_ex<-WTCPKSC_in[!grep("PKS",WTCPKSC_in$isolates),]
##WT15 and PKS15
WT15PKS15_in<-svcallers_parallel[grepl("PKS15.*WT15|WT15.*PKS15",svcallers_parallel$isolates),] ##has WT15 and PKS15
#WT15PKS15_ex<-WT15PKS15_in[!grep("PKS",WT15PKS15_in$isolates),]
##WTC and WT15
WTCWT15_in<-svcallers_parallel[grepl("WTC.*WT15|WT15.*WTC",svcallers_parallel$isolates),] ##has WT
#WTCWT15_ex<-WTCWT15_in[!grep("PKS",WTCWT15_in$isolates),]
##PKSC and PKS15
PKSCPKS15_in<-svcallers_parallel[grepl(c("PKSC.*PKS15|PKS15.*PKSC"),svcallers_parallel$isolates),] ##has PKS
#PKSCPKS15_ex<-PKSCPKS15_in[!grep("WT",PKSCPKS15_in$isolates),]

##WTC
WTC_in<-svcallers_parallel[grepl("WTC.*WTC",svcallers_parallel$isolates),] ##has WTC
#WTC_ex<-WTC_in[!grep("15|PKS",WTC_in$isolates),] ##only WTC
##WT15
WT15_in<-svcallers_parallel[grepl("WT15.*WT15",svcallers_parallel$isolates),] ##has WT15
#WT15_ex<-WT15_in[!grep("WTC|PKS",WT15_in$isolates),] ##only WT15
##PKSC
PKSC_in<-svcallers_parallel[grepl("PKSC.*PKSC",svcallers_parallel$isolates),] ##has PKSC
#PKSC_ex<-PKSC_in[!grep("15|WT",PKSC_in$isolates),] ##only PKS
##PKS15
PKS15_in<-svcallers_parallel[grepl("PKS15.*PKS15",svcallers_parallel$isolates),] ##has PKS15
#PKS15_ex<-PKS15_in[!grep("PKSC|WT",PKS15_in$isolates),] ##only PKS15

write.csv(WTCPKSC_in, "WTCPKSC_in.csv")
write.csv(WT15PKS15_in, "WT15PKS15_in.csv")
write.csv(WTCWT15_in, "WTCWT15_in.csv")
write.csv(PKSCPKS15_in, "PKSCPKS15_in.csv")
write.csv(WTC_in, "WTC_in.csv")
write.csv(WT15_in, "WT15_in.csv")
write.csv(PKSC_in, "PKSC_in.csv")
write.csv(PKS15_in, "PKS15_in.csv")

##################################################
#### Read in SV genes and annotate ###############
##################################################
C.group_genes<-read.table("WT_PKS_bedfiles/C.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
evol15.group_genes<-read.table("WT_PKS_bedfiles/15.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
WT.group_genes<-read.table("WT_PKS_bedfiles/WT.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
PKS.group_genes<-read.table("WT_PKS_bedfiles/PKS.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))

WTC.group_genes<-read.table("WT_PKS_bedfiles/WTC.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
PKSC.group_genes<-read.table("WT_PKS_bedfiles/PKSC.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
WT15.group_genes<-read.table("WT_PKS_bedfiles/WT15.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
PKS15.group_genes<-read.table("WT_PKS_bedfiles/PKS15.group_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))

C.parallel_genes<-read.table("WT_PKS_bedfiles/C.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
evol15.parallel_genes<-read.table("WT_PKS_bedfiles/15.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
WT.parallel_genes<-read.table("WT_PKS_bedfiles/WT.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
PKS.parallel_genes<-read.table("WT_PKS_bedfiles/PKS.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))

WTC.parallel_genes<-read.table("WT_PKS_bedfiles/WTC.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
PKSC.parallel_genes<-read.table("WT_PKS_bedfiles/PKSC.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
WT15.parallel_genes<-read.table("WT_PKS_bedfiles/WT15.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))
PKS15.parallel_genes<-read.table("WT_PKS_bedfiles/PKS15.parallel_genes.bed",col.names=c("CHROM","start","stop","gene","note","svtype"))

C.group_annot<-merge(C.group_genes,annot,by="gene")
evol15.group_annot<-merge(evol15.group_genes,annot,by="gene")
WT.group_annot<-merge(WT.group_genes,annot,by="gene")
PKS.group_annot<-merge(PKS.group_genes,annot,by="gene")

WTC.group_annot<-merge(WTC.group_genes,annot,by="gene")
PKSC.group_annot<-merge(PKSC.group_genes,annot,by="gene")
WT15.group_annot<-merge(WT15.group_genes,annot,by="gene")
PKS15.group_annot<-merge(PKS15.group_genes,annot,by="gene")
write.csv(WTC.group_annot, "WTC.group_annot.csv")
write.csv(WT15.group_annot, "WT15.group_annot.csv")
write.csv(PKSC.group_annot, "PKSC.group_annot.csv")
write.csv(PKS15.group_annot, "PKS15.group_annot.csv")

C.parallel_annot<-merge(C.parallel_genes,annot,by="gene")
evol15.parallel_annot<-merge(evol15.parallel_genes,annot,by="gene")
WT.parallel_annot<-merge(WT.parallel_genes,annot,by="gene")
PKS.parallel_annot<-merge(PKS.parallel_genes,annot,by="gene")

WTC.parallel_annot<-merge(WTC.parallel_genes,annot,by="gene")
PKSC.parallel_annot<-merge(PKSC.parallel_genes,annot,by="gene")
WT15.parallel_annot<-merge(WT15.parallel_genes,annot,by="gene")
PKS15.parallel_annot<-merge(PKS15.parallel_genes,annot,by="gene")
write.csv(WTC.parallel_annot, "WTC.parallel_annot.csv")
write.csv(WT15.parallel_annot, "WT15.parallel_annot.csv")
write.csv(PKSC.parallel_annot, "PKSC.parallel_annot.csv")
write.csv(PKS15.parallel_annot, "PKS15.parallel_annot.csv")

###############################################
###### Length of each SV type #################
###############################################
WTC.group_genes$Group<-"WTC"
PKSC.group_genes$Group<-"PKSC"
WT15.group_genes$Group<-"WT15"
PKS15.group_genes$Group<-"PKS15"

Groups.bed<-rbind(WTC.group_genes,PKSC.group_genes,WT15.group_genes,PKS15.group_genes)
Groups.bed<-unique(Groups.bed[,c("CHROM","start","stop","gene","svtype","Group")])
Groups.bed$condition<-gsub(".*C","control",Groups.bed$Group)
Groups.bed$condition<-gsub(".*15","evolved",Groups.bed$condition)
Groups.bed$length<-Groups.bed$stop-Groups.bed$start


summary(Groups.bed[Groups.bed$svtype=="DEL",]$length)
summary(Groups.bed[Groups.bed$svtype=="DUP",]$length)
summary(Groups.bed[Groups.bed$svtype=="INV",]$length)

summary(Groups.bed[Groups.bed$svtype=="INV" & Groups.bed$Group=="WTC",]$length)
summary(Groups.bed[Groups.bed$svtype=="INV" & Groups.bed$Group=="WT15",]$length)
summary(Groups.bed[Groups.bed$svtype=="INV" & Groups.bed$Group=="PKSC",]$length)
summary(Groups.bed[Groups.bed$svtype=="INV" & Groups.bed$Group=="PKS15",]$length)

my_comparisons <- list( c("WTC", "WT15"), c("PKSC", "PKS15"), 
                        c("WTC", "PKSC"), c("WT15", "PKS15"))

#####################################
#### Identify Inherited SVs #########
#####################################
###a file that has loci with multiple SV mutations, must verify if truly parallel
svcallers_inherit<-svcallers_parallel[grep("orig",svcallers_parallel$isolates),]
svcallers_inherit$anc<-gsub(".*orig-","",svcallers_inherit$isolates)
svcallers_inherit$desc<-gsub(",orig-.*","",svcallers_inherit$isolates)
matches<-sapply(1:nrow(svcallers_inherit),function(i)agrepl(svcallers_inherit$anc[i],svcallers_inherit$desc[i]))
svcallers_inherit<-svcallers_inherit[matches,]
###none are inherited and parallel, only inherited

###############################################################################
############# Annotate Gene files of SNPs, SVs, Transcriptome #################
###############################################################################
### Read SNP Data
WTC.snps<-read.csv(paste0(DNA_dir,"parallel_genes.WTC.csv"))
PKSC.snps<-read.csv(paste0(DNA_dir,"parallel_genes.PKSC.csv"))
WT15.snps<-read.csv(paste0(DNA_dir,"parallel_genes.WT15.csv"))
PKS15.snps<-read.csv(paste0(DNA_dir,"parallel_genes.PKS15.csv"))
### Read SV Data
WTC.svs<-read.csv(paste0(DNA_dir,"WTC.parallel_annot.csv"))
PKSC.svs<-read.csv(paste0(DNA_dir,"PKSC.parallel_annot.csv"))
WT15.svs<-read.csv(paste0(DNA_dir,"WT15.parallel_annot.csv"))
PKS15.svs<-read.csv(paste0(DNA_dir,"PKS15.parallel_annot.csv"))
### Read Transcriptome Data
WTC.irr<-read.csv(paste0(RNA_dir,"res_WT_irr_vs_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKSC.irr<-read.csv(paste0(RNA_dir,"res_PKS_vs_WT_irr_vs_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.irr<-read.csv(paste0(RNA_dir,"res_WT15common_vs_WT_irr.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.irr<-read.csv(paste0(RNA_dir,"res_PKS15common_vs_PKS_irr.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.1.2.irr<-read.csv(paste0(RNA_dir,"res_WT15.1.2_vs_WT_irr.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.2.2.irr<-read.csv(paste0(RNA_dir,"res_WT15.2.2_vs_WT_irr.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.4.2.irr<-read.csv(paste0(RNA_dir,"res_WT15.4.2_vs_WT_irr.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.1.1.irr<-read.csv(paste0(RNA_dir,"res_PKS15.1.1_vs_PKS_irr.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.2.2.irr<-read.csv(paste0(RNA_dir,"res_PKS15.2.2_vs_PKS_irr.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.3.2.irr<-read.csv(paste0(RNA_dir,"res_PKS15.3.2_vs_PKS_irr.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.4.2.irr<-read.csv(paste0(RNA_dir,"res_PKS15.4.2_vs_PKS_irr.GO_KEGG_KOG_IPR_annotation.csv"))
WTC.ctl<-read.csv(paste0(RNA_dir,"res_WT_irr_vs_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKSC.ctl<-read.csv(paste0(RNA_dir,"res_PKS_vs_WT_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.ctl<-read.csv(paste0(RNA_dir,"res_WT15common_vs_WT_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.ctl<-read.csv(paste0(RNA_dir,"res_PKS15common_vs_PKS_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.1.2.ctl<-read.csv(paste0(RNA_dir,"res_WT15.1.2_vs_WT_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.2.2.ctl<-read.csv(paste0(RNA_dir,"res_WT15.2.2_vs_WT_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
WT15.4.2.ctl<-read.csv(paste0(RNA_dir,"res_WT15.4.2_vs_WT_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.1.1.ctl<-read.csv(paste0(RNA_dir,"res_PKS15.1.1_vs_PKS_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.2.2.ctl<-read.csv(paste0(RNA_dir,"res_PKS15.2.2_vs_PKS_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.3.2.ctl<-read.csv(paste0(RNA_dir,"res_PKS15.3.2_vs_PKS_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
PKS15.4.2.ctl<-read.csv(paste0(RNA_dir,"res_PKS15.4.2_vs_PKS_ctl.GO_KEGG_KOG_IPR_annotation.csv"))
colnames(WTC.irr)[2]<-"genename"
colnames(PKSC.irr)[2]<-"genename"
colnames(WT15.irr)[2]<-"genename"
colnames(PKS15.irr)[2]<-"genename"
colnames(WTC.ctl)[2]<-"genename"
colnames(PKSC.ctl)[2]<-"genename"
colnames(WT15.ctl)[2]<-"genename"
colnames(PKS15.ctl)[2]<-"genename"
colnames(WT15.1.2.irr)[2]<-"genename"
colnames(WT15.2.2.irr)[2]<-"genename"
colnames(WT15.4.2.irr)[2]<-"genename"
colnames(PKS15.1.1.irr)[2]<-"genename"
colnames(PKS15.2.2.irr)[2]<-"genename"
colnames(PKS15.3.2.irr)[2]<-"genename"
colnames(PKS15.4.2.irr)[2]<-"genename"
colnames(WTC.ctl)[2]<-"genename"
colnames(PKSC.ctl)[2]<-"genename"
colnames(WT15.ctl)[2]<-"genename"
colnames(PKS15.ctl)[2]<-"genename"
colnames(WT15.1.2.ctl)[2]<-"genename"
colnames(WT15.2.2.ctl)[2]<-"genename"
colnames(WT15.4.2.ctl)[2]<-"genename"
colnames(PKS15.1.1.ctl)[2]<-"genename"
colnames(PKS15.2.2.ctl)[2]<-"genename"
colnames(PKS15.3.2.ctl)[2]<-"genename"
colnames(PKS15.4.2.ctl)[2]<-"genename"

WTC.irr<-merge(WTC.irr,annot)
PKSC.irr<-merge(PKSC.irr,annot)
WT15.irr<-merge(WT15.irr,annot)
PKS15.irr<-merge(PKS15.irr,annot)
WTC.ctl<-merge(WTC.ctl,annot)
PKSC.ctl<-merge(PKSC.ctl,annot)
WT15.ctl<-merge(WT15.ctl,annot)
PKS15.ctl<-merge(PKS15.ctl,annot)

WTC.irr<-merge(WTC.irr,annot)
PKSC.irr<-merge(PKSC.irr,annot)
WT15.irr<-merge(WT15.irr,annot)
PKS15.irr<-merge(PKS15.irr,annot)
WTC.ctl<-merge(WTC.ctl,annot)
PKSC.ctl<-merge(PKSC.ctl,annot)
WT15.ctl<-merge(WT15.ctl,annot)
PKS15.ctl<-merge(PKS15.ctl,annot)
WT15.1.2.irr<-merge(WT15.1.2.irr,annot)
WT15.2.2.irr<-merge(WT15.2.2.irr,annot)
WT15.4.2.irr<-merge(WT15.4.2.irr,annot)
PKS15.1.1.irr<-merge(PKS15.1.1.irr,annot)
PKS15.2.2.irr<-merge(PKS15.2.2.irr,annot)
PKS15.3.2.irr<-merge(PKS15.3.2.irr,annot)
PKS15.4.2.irr<-merge(PKS15.4.2.irr,annot)
WTC.ctl<-merge(WTC.ctl,annot)
PKSC.ctl<-merge(PKSC.ctl,annot)
WT15.ctl<-merge(WT15.ctl,annot)
PKS15.ctl<-merge(PKS15.ctl,annot)
WT15.1.2.ctl<-merge(WT15.1.2.ctl,annot)
WT15.2.2.ctl<-merge(WT15.2.2.ctl,annot)
WT15.4.2.ctl<-merge(WT15.4.2.ctl,annot)
PKS15.1.1.ctl<-merge(PKS15.1.1.ctl,annot)
PKS15.2.2.ctl<-merge(PKS15.2.2.ctl,annot)
PKS15.3.2.ctl<-merge(PKS15.3.2.ctl,annot)
PKS15.4.2.ctl<-merge(PKS15.4.2.ctl,annot)

###############################################################################
############ Intersect Gene files of SNPs, SVs, Transcriptome #################
###############################################################################
### SNPs, SVs
WTC.snps.svs<-merge(WTC.snps,WTC.svs, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.snps.svs<-merge(WT15.snps,WT15.svs, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.snps.svs<-merge(PKSC.snps,PKSC.svs, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.snps.svs<-merge(PKS15.snps,PKS15.svs, by=c("gene","goName","pathway","kogdefline","iprDesc"))
### SNPs, RNA
WTC.snps.irr<-merge(WTC.snps,WTC.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.snps.irr<-merge(WT15.snps,WT15.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.snps.irr<-merge(PKSC.snps,PKSC.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.snps.irr<-merge(PKS15.snps,PKS15.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WTC.snps.ctl<-merge(WTC.snps,WTC.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.snps.ctl<-merge(WT15.snps,WT15.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.snps.ctl<-merge(PKSC.snps,PKSC.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.snps.ctl<-merge(PKS15.snps,PKS15.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))

WT15.1.2.snps.irr<-merge(WT15.snps,WT15.1.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.2.2.snps.irr<-merge(WT15.snps,WT15.2.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.4.2.snps.irr<-merge(WT15.snps,WT15.4.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.1.1.snps.irr<-merge(PKS15.snps,PKS15.1.1.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.2.2.snps.irr<-merge(PKS15.snps,PKS15.2.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.3.2.snps.irr<-merge(PKS15.snps,PKS15.3.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.4.2.snps.irr<-merge(PKS15.snps,PKS15.4.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))

WT15.1.2.snps.ctl<-merge(WT15.snps,WT15.1.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.2.2.snps.ctl<-merge(WT15.snps,WT15.2.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.4.2.snps.ctl<-merge(WT15.snps,WT15.4.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.1.1.snps.ctl<-merge(PKS15.snps,PKS15.1.1.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.2.2.snps.ctl<-merge(PKS15.snps,PKS15.2.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.3.2.snps.ctl<-merge(PKS15.snps,PKS15.3.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.4.2.snps.ctl<-merge(PKS15.snps,PKS15.4.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))


### SVs, RNA
WTC.svs.irr<-merge(WTC.svs,WTC.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.svs.irr<-merge(WT15.svs,WT15.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.svs.irr<-merge(PKSC.svs,PKSC.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.svs.irr<-merge(PKS15.svs,PKS15.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.1.2.svs.irr<-merge(WT15.svs,WT15.1.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.2.2.svs.irr<-merge(WT15.svs,WT15.2.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.4.2.svs.irr<-merge(WT15.svs,WT15.4.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.1.1.svs.irr<-merge(PKS15.svs,PKS15.1.1.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.2.2.svs.irr<-merge(PKS15.svs,PKS15.2.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.3.2.svs.irr<-merge(PKS15.svs,PKS15.3.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.4.2.svs.irr<-merge(PKS15.svs,PKS15.4.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))

WTC.svs.ctl<-merge(WTC.svs,WTC.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.svs.ctl<-merge(WT15.svs,WT15.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.svs.ctl<-merge(PKSC.svs,PKSC.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.svs.ctl<-merge(PKS15.svs,PKS15.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.1.2.svs.ctl<-merge(WT15.svs,WT15.1.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.2.2.svs.ctl<-merge(WT15.svs,WT15.2.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.4.2.svs.ctl<-merge(WT15.svs,WT15.4.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.1.1.svs.ctl<-merge(PKS15.svs,PKS15.1.1.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.2.2.svs.ctl<-merge(PKS15.svs,PKS15.2.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.3.2.svs.ctl<-merge(PKS15.svs,PKS15.3.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.4.2.svs.ctl<-merge(PKS15.svs,PKS15.4.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))

### SNPs, SVs, RNA
WTC.snps.svs.irr<-merge(WTC.snps.svs,WTC.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.snps.svs.irr<-merge(WT15.snps.svs,WT15.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.snps.svs.irr<-merge(PKSC.snps.svs,PKSC.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.snps.svs.irr<-merge(PKS15.snps.svs,PKS15.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.1.2.snps.svs.irr<-merge(WT15.snps.svs,WT15.1.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.2.2.snps.svs.irr<-merge(WT15.snps.svs,WT15.2.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.4.2.snps.svs.irr<-merge(WT15.snps.svs,WT15.4.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.1.1.snps.svs.irr<-merge(PKS15.snps.svs,PKS15.1.1.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.2.2.snps.svs.irr<-merge(PKS15.snps.svs,PKS15.2.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.3.2.snps.svs.irr<-merge(PKS15.snps.svs,PKS15.3.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.4.2.snps.svs.irr<-merge(PKS15.snps.svs,PKS15.4.2.irr, by=c("gene","goName","pathway","kogdefline","iprDesc"))

WTC.snps.svs.ctl<-merge(WTC.snps.svs,WTC.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.snps.svs.ctl<-merge(WT15.snps.svs,WT15.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKSC.snps.svs.ctl<-merge(PKSC.snps.svs,PKSC.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.snps.svs.ctl<-merge(PKS15.snps.svs,PKS15.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.1.2.snps.svs.ctl<-merge(WT15.snps.svs,WT15.1.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.2.2.snps.svs.ctl<-merge(WT15.snps.svs,WT15.2.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
WT15.4.2.snps.svs.ctl<-merge(WT15.snps.svs,WT15.4.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.1.1.snps.svs.ctl<-merge(PKS15.snps.svs,PKS15.1.1.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.2.2.snps.svs.ctl<-merge(PKS15.snps.svs,PKS15.2.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.3.2.snps.svs.ctl<-merge(PKS15.snps.svs,PKS15.3.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))
PKS15.4.2.snps.svs.ctl<-merge(PKS15.snps.svs,PKS15.4.2.ctl, by=c("gene","goName","pathway","kogdefline","iprDesc"))

###############################################################################
############ Go Enrichment: SNPs, SVs, Transcriptome ##########################
###############################################################################

GOenrich<-function(table,file){
  write.csv(table,file = paste0(file,".csv"))
  table$gene<-as.character(table$gene)
  gene<-unique(table$gene)
  plots <- vector('list', 2)
  ego<-enricher(gene = gene, 
                universe = annot$gene,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                TERM2GENE = term2gene_go[,c("goName","gene")]) ###goName has 1000s and keytypes 3 categories (MF, BP, CC)
  a1<-dotplot(ego,orderBy="GeneRatio")+ggtitle("GO")
  a1<-tryCatch({print(a1)}, error = function(e){
    a1<-NA
    return(a1)})
  write.csv(ego,file = paste0(file,".go.csv"))
  
  ekegg <- enricher(gene = gene,
                    universe = annot$gene,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    TERM2GENE = term2gene_kegg[,c("pathway","gene")]) ###pathway and pathway class the same
  b1<-dotplot(ekegg,orderBy="GeneRatio")+ggtitle("KEGG")
  b1<-tryCatch({print(b1)}, error = function(e){
    b1<-NA
    return(b1)})
  write.csv(ekegg,file = paste0(file,".kegg.csv"))
  
  ekog <- enricher(gene = gene,
                   universe = annot$gene,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   TERM2GENE = term2gene_kog[,c("kogClass","gene")]) ###kogdefline 1000s kogClass 25
  c1<-dotplot(ekog,orderBy="GeneRatio")+ggtitle("KOG")
  c1<-tryCatch({print(c1)}, error = function(e){
    c1<-NA
    return(c1)})
  write.csv(ekog,file = paste0(file,".kog.csv"))
  #eipr <- enricher(gene = gene,
  #                 universe = annot$gene,
  #                 pAdjustMethod = "BH",
  #                 qvalueCutoff = 0.05,
  #                 TERM2GENE = term2gene_ipr[,c("domainDesc","gene")]) ### iprDesc and domainDesc 1000s but domainDesc more
  #d1<-dotplot(eipr,orderBy="GeneRatio")+ggtitle("ipr")
  #d1<-tryCatch({print(d1)}, error = function(e){
  #  d1<-NA
  #  return(d1)})
  
  plots<-annotate_figure(ggarrange(a1,b1,c1,common.legend = TRUE))
  pdf(paste0(file,".pdf"),width=32,height=44)
  print(plots)
  graphics.off()
}

GOenrich(table=WT15.snps.svs.irr,file="WT15.snps.svs.irr")
GOenrich(table=PKS15.snps.svs.irr,file="PKS15.snps.svs.irr")
GOenrich(table=WTC.snps.svs.irr,file="WTC.snps.svs.irr")
GOenrich(table=PKSC.snps.svs.irr,file="PKSC.snps.svs.irr")
GOenrich(table=WT15.snps.svs.ctl,file="WT15.snps.svs.ctl")
GOenrich(table=PKS15.snps.svs.ctl,file="PKS15.snps.svs.ctl")
GOenrich(table=WTC.snps.svs.ctl,file="WTC.snps.svs.ctl")
GOenrich(table=PKSC.snps.svs.ctl,file="PKSC.snps.svs.ctl")
GOenrich(table=WT15.1.2.snps.svs.irr,file="WT15.1.2.snps.svs.irr")
GOenrich(table=WT15.2.2.snps.svs.irr,file="WT15.2.2.snps.svs.irr")
GOenrich(table=WT15.4.2.snps.svs.irr,file="WT15.4.2.snps.svs.irr")
GOenrich(table=PKS15.1.1.snps.svs.irr,file="PKS15.1.1.snps.svs.irr")
GOenrich(table=PKS15.2.2.snps.svs.irr,file="PKS15.2.2.snps.svs.irr")
GOenrich(table=PKS15.3.2.snps.svs.irr,file="PKS15.3.2.snps.svs.irr")
GOenrich(table=PKS15.4.2.snps.svs.irr,file="PKS15.4.2.snps.svs.irr")
GOenrich(table=WT15.1.2.snps.svs.ctl,file="WT15.1.2.snps.svs.ctl")
GOenrich(table=WT15.2.2.snps.svs.ctl,file="WT15.2.2.snps.svs.ctl")
GOenrich(table=WT15.4.2.snps.svs.ctl,file="WT15.4.2.snps.svs.ctl")
GOenrich(table=PKS15.1.1.snps.svs.ctl,file="PKS15.1.1.snps.svs.ctl")
GOenrich(table=PKS15.2.2.snps.svs.ctl,file="PKS15.2.2.snps.svs.ctl")
GOenrich(table=PKS15.3.2.snps.svs.ctl,file="PKS15.3.2.snps.svs.ctl")
GOenrich(table=PKS15.4.2.snps.svs.ctl,file="PKS15.4.2.snps.svs.ctl")

GOenrich(table=WT15.snps.irr,file="WT15.snps.irr")
GOenrich(table=PKS15.snps.irr,file="PKS15.snps.irr")
GOenrich(table=WTC.snps.irr,file="WTC.snps.irr")
GOenrich(table=PKSC.snps.irr,file="PKSC.snps.irr")
GOenrich(table=WT15.snps.ctl,file="WT15.snps.ctl")
GOenrich(table=PKS15.snps.ctl,file="PKS15.snps.ctl")
GOenrich(table=WTC.snps.ctl,file="WTC.snps.ctl")
GOenrich(table=PKSC.snps.ctl,file="PKSC.snps.ctl")
GOenrich(table=WT15.1.2.snps.irr,file="WT15.1.2.snps.irr")
GOenrich(table=WT15.2.2.snps.irr,file="WT15.2.2.snps.irr")
GOenrich(table=WT15.4.2.snps.irr,file="WT15.4.2.snps.irr")
GOenrich(table=PKS15.1.1.snps.irr,file="PKS15.1.1.snps.irr")
GOenrich(table=PKS15.2.2.snps.irr,file="PKS15.2.2.snps.irr")
GOenrich(table=PKS15.3.2.snps.irr,file="PKS15.3.2.snps.irr")
GOenrich(table=PKS15.4.2.snps.irr,file="PKS15.4.2.snps.irr")
GOenrich(table=WT15.1.2.snps.ctl,file="WT15.1.2.snps.ctl")
GOenrich(table=WT15.2.2.snps.ctl,file="WT15.2.2.snps.ctl")
GOenrich(table=WT15.4.2.snps.ctl,file="WT15.4.2.snps.ctl")
GOenrich(table=PKS15.1.1.snps.ctl,file="PKS15.1.1.snps.ctl")
GOenrich(table=PKS15.2.2.snps.ctl,file="PKS15.2.2.snps.ctl")
GOenrich(table=PKS15.3.2.snps.ctl,file="PKS15.3.2.snps.ctl")
GOenrich(table=PKS15.4.2.snps.ctl,file="PKS15.4.2.snps.ctl")

#### Interesting terms: SVs vs transcriptome
GOenrich(table=WT15.svs.irr,file="WT15.svs.irr")
GOenrich(table=PKS15.svs.irr,file="PKS15.svs.irr")
GOenrich(table=WTC.svs.irr,file="WTC.svs.irr")
GOenrich(table=PKSC.svs.irr,file="PKSC.svs.irr")
GOenrich(table=WT15.svs.ctl,file="WT15.svs.ctl")
GOenrich(table=PKS15.svs.ctl,file="PKS15.svs.ctl")
GOenrich(table=WTC.svs.ctl,file="WTC.svs.ctl")
GOenrich(table=PKSC.svs.ctl,file="PKSC.svs.ctl")
GOenrich(table=WT15.1.2.svs.irr,file="WT15.1.2.svs.irr")
GOenrich(table=WT15.2.2.svs.irr,file="WT15.2.2.svs.irr")
GOenrich(table=WT15.4.2.svs.irr,file="WT15.4.2.svs.irr")
GOenrich(table=PKS15.1.1.svs.irr,file="PKS15.1.1.svs.irr")
GOenrich(table=PKS15.2.2.svs.irr,file="PKS15.2.2.svs.irr")
GOenrich(table=PKS15.3.2.svs.irr,file="PKS15.3.2.svs.irr")
GOenrich(table=PKS15.4.2.svs.irr,file="PKS15.4.2.svs.irr")
GOenrich(table=WT15.1.2.svs.ctl,file="WT15.1.2.svs.ctl")
GOenrich(table=WT15.2.2.svs.ctl,file="WT15.2.2.svs.ctl")
GOenrich(table=WT15.4.2.svs.ctl,file="WT15.4.2.svs.ctl")
GOenrich(table=PKS15.1.1.svs.ctl,file="PKS15.1.1.svs.ctl")
GOenrich(table=PKS15.2.2.svs.ctl,file="PKS15.2.2.svs.ctl")
GOenrich(table=PKS15.3.2.svs.ctl,file="PKS15.3.2.svs.ctl")
GOenrich(table=PKS15.4.2.svs.ctl,file="PKS15.4.2.svs.ctl")

GOenrich(table=WT15.snps.svs,file="WT15.snps.svs")
GOenrich(table=PKS15.snps.svs,file="PKS15.snps.svs")
GOenrich(table=WTC.snps.svs,file="WTC.snps.svs")
GOenrich(table=PKSC.snps.svs,file="PKSC.snps.svs")


################################################################################
####################### WT15 and PKS15 Manuscript Figures ######################
################################################################################

############################################# Figure 2 ################################################
SNPSV_8000Gy_lineage<-ggarrange(SNPaccum_8000Gy_lineage+
                                  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
                                        axis.title.x = element_blank(),
                                        panel.background = element_blank()),
                                SVaccum_8000Gy_lineage+
                                  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
                                        panel.background = element_blank()), 
                                ncol=1,nrow=2, common.legend = TRUE,legend = "bottom",heights = c(1,2))

pdf("evolFig2_SNPSVaccum_8000Gy_lineage.pdf",width=8.5, height=11)
SNPSV_8000Gy_lineage
graphics.off()
kw.SVaccum_lineage
dunn.SVaccum_lineage

############################################# Figure S1 ################################################
colnames(SNPaccum_df)<-c("alt","Sample","condition","Group","lineage")
SNPs_sample2<-rbind(SNPaccum_df,SNP_orig)
SNPs_sample2$Type<-"SNP"
SNPs_sample2$lineage<-NULL
colnames(SNPs_sample2)<-c("count_1","Sample","condition","Group","Type")

SVs_sample2<-rbind(SVaccum_df[c("accumulated_SVs","sample","Group","svtype","condition")],svcallers.df_orig)
SVs_sample2<-SVs_sample2[c("accumulated_SVs","sample","condition","Group","svtype")]
colnames(SVs_sample2)<-c("count_1","Sample","condition","Group","Type")

SNP_SV.df<-rbind(SNPs_sample2,SVs_sample2)
SNP_SV.df$Sample<-gsub("-",".",SNP_SV.df$Sample)
SNP_SV.df$Type<-factor(SNP_SV.df$Type, levels = c("SNP", "DEL", "DUP","INV","TRA"))

SNP_SV.df[SNP_SV.df$Sample=="orig.WTC.1.4",]$Sample<-"orig.WTC.15"
SNP_SV.df[SNP_SV.df$Sample=="orig.PKSC.1.4",]$Sample<-"orig.PKSC.15"

median(SNP_SV.df[SNP_SV.df$Group=="WT15" & SNP_SV.df$Type=="TRA",]$count_1)
median(SNP_SV.df[SNP_SV.df$Group=="PKS15" & SNP_SV.df$Type=="TRA",]$count_1)
### Counts of SNPs/SVs per Sample
mutation_counts_sample<-ggplot(SNP_SV.df, aes(x=Sample, y=count_1, fill=condition)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(~Type, scales="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.background = element_blank())+
  scale_fill_manual(breaks = c("control", "evolved",
                               "orig-control","orig-evolved"), 
                    values=c("blue", "red","gray","black"))+
  xlab("Isolate")+ylab("No. of mutations") + scale_x_discrete(limits=rev)+coord_flip()

pdf("evolFigS1_SNPsSVs_per_sample.pdf",width=8.5, height=11,onefile = FALSE)
mutation_counts_sample
graphics.off()


############################################# Table S3 ################################################

median(Groups.bed[Groups.bed$svtype=="DEL",]$length,na.rm = TRUE) #1762
median(Groups.bed[Groups.bed$svtype=="DUP",]$length,na.rm = TRUE) #1787
median(Groups.bed[Groups.bed$svtype=="INV",]$length,na.rm = TRUE) #1924

range(Groups.bed[Groups.bed$svtype=="DEL",]$length,na.rm = TRUE) #182-13,043
range(Groups.bed[Groups.bed$svtype=="DUP",]$length,na.rm = TRUE) #197-12,338
range(Groups.bed[Groups.bed$svtype=="INV",]$length,na.rm = TRUE) #182-18,910

############################################# Table 1 #################################################

### SNPs and RNA
WT15.snps.irr$profile<-"WT15common_vs_WT_irr"
PKS15.snps.irr$profile<-"PKS15common_vs_PKS_irr"
WT15.snps.ctl$profile<-"WT15common_vs_WT_ctl"
PKS15.snps.ctl$profile<-"PKS15common_vs_PKS_ctl"

WT15.1.2.snps.irr$profile<-"WT15.1.2_vs_WT_irr"
WT15.2.2.snps.irr$profile<-"WT15.2.2_vs_WT_irr"
WT15.4.2.snps.irr$profile<-"WT15.4.2_vs_WT_irr"
PKS15.1.1.snps.irr$profile<-"PKS15.1.2_vs_PKS_irr"
PKS15.2.2.snps.irr$profile<-"PKS15.2.2_vs_PKS_irr"
PKS15.3.2.snps.irr$profile<-"PKS15.3.2_vs_PKS_irr"
PKS15.4.2.snps.irr$profile<-"PKS15.4.2_vs_PKS_irr"

WT15.1.2.snps.ctl$profile<-"WT15.1.2_vs_WT_ctl"
WT15.2.2.snps.ctl$profile<-"WT15.2.2_vs_WT_ctl"
WT15.4.2.snps.ctl$profile<-"WT15.4.2_vs_WT_ctl"
PKS15.1.1.snps.ctl$profile<-"PKS15.1.1_vs_PKS_ctl"
PKS15.2.2.snps.ctl$profile<-"PKS15.2.2_vs_PKS_ctl"
PKS15.3.2.snps.ctl$profile<-"PKS15.3.2_vs_PKS_ctl"
PKS15.4.2.snps.ctl$profile<-"PKS15.4.2_vs_PKS_ctl"

#create ddataframe of parallel genes/regions with SNP data and transcriptome data
AncDesc_regions<-rbind(#WT15.snps.irr,
  #PKS15.snps.irr,
  #WT15.snps.ctl,
  #PKS15.snps.ctl,
  WT15.1.2.snps.irr,
  WT15.2.2.snps.irr,
  WT15.4.2.snps.irr,
  PKS15.1.1.snps.irr,
  PKS15.2.2.snps.irr,
  PKS15.3.2.snps.irr,
  PKS15.4.2.snps.irr,
  WT15.1.2.snps.ctl,
  WT15.2.2.snps.ctl,
  WT15.4.2.snps.ctl,
  PKS15.1.1.snps.ctl,
  PKS15.2.2.snps.ctl,
  PKS15.3.2.snps.ctl,
  PKS15.4.2.snps.ctl
)

###annotate 
AncDesc_regions$annotation<-paste(AncDesc_regions$goName,";",
                                  AncDesc_regions$pathway,";",
                                  AncDesc_regions$kogdefline,";",
                                  AncDesc_regions$iprDesc)

###match the parallel gene with the transcriptome profiles by isolate
AncDesc_regions<-AncDesc_regions[which(grepl("WT15.1.2",AncDesc_regions$Isolate.x) & grepl("WT15.1.2",AncDesc_regions$profile)|
                                         grepl("WT15.2.2",AncDesc_regions$Isolate.x) & grepl("WT15.2.2",AncDesc_regions$profile)|
                                         grepl("WT15.4.2",AncDesc_regions$Isolate.x) & grepl("WT15.4.2",AncDesc_regions$profile)|
                                         grepl("PKS15.1.1",AncDesc_regions$Isolate.x) & grepl("PKS15.1.1",AncDesc_regions$profile)|
                                         grepl("PKS15.2.2",AncDesc_regions$Isolate.x) & grepl("PKS15.2.2",AncDesc_regions$profile)|
                                         grepl("PKS15.3.2",AncDesc_regions$Isolate.x) & grepl("PKS15.3.2",AncDesc_regions$profile)|
                                         grepl("PKS15.4.2",AncDesc_regions$Isolate.x) & grepl("PKS15.4.2",AncDesc_regions$profile)),]

###subset for protein coding region
AncDesc_gene<-AncDesc_regions[grep("NON_SYN|STOP|START|FRAME|UTR",AncDesc_regions$Type),]

###make sure to get mutations that change the protein coding regions for all parallel mutations
AncDesc_gene0<-AncDesc_gene[which(!grepl(",",AncDesc_gene$Type)),]
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("NON",AncDesc_gene$Type)& grepl(",SYN",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("NON",AncDesc_gene$Type)& grepl("FRAME",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("NON",AncDesc_gene$Type)& grepl("STOP",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("NON",AncDesc_gene$Type)& grepl("START",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("NON",AncDesc_gene$Type)& grepl("UTR",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("SYN",AncDesc_gene$Type)& grepl("FRAME",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("SYN",AncDesc_gene$Type)& grepl("STOP",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("SYN",AncDesc_gene$Type)& grepl("START",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("SYN",AncDesc_gene$Type)& grepl("UTR",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("FRAME",AncDesc_gene$Type)& grepl("STOP",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("FRAME",AncDesc_gene$Type)& grepl("START",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("FRAME",AncDesc_gene$Type)& grepl("UTR",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("START",AncDesc_gene$Type)& grepl("STOP",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("START",AncDesc_gene$Type)& grepl("UTR",AncDesc_gene$Type)),])
AncDesc_gene0<-rbind(AncDesc_gene0,AncDesc_gene[which(grepl("STOP",AncDesc_gene$Type)& grepl("UTR",AncDesc_gene$Type)),])

###Ancestral,Parallel, around protein coding regions
AncParallel_regions<-parallel_genes.aroundcoding[c("genename.x","Anc")] ##reformat file containing parallel genes at ancestral sequences
colnames(AncParallel_regions)<-c("genename","Anc")
AncParallel_regions<-merge(AncParallel_regions,AncDesc_regions,by="genename")

###Ancestral,Parallel, around protein coding regions; but must have at least NS, S, Shift,Stop,Start,UTR
AncDesc_regions0<-AncDesc_regions[which(!grepl(",",AncDesc_regions$Type)),]
AncDesc_regions0<-AncDesc_regions0[grep("NON_SYN|STOP|START|FRAME|UTR",AncDesc_regions0$Type),]
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("NON",AncDesc_regions$Type)& grepl(",SYN",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("NON",AncDesc_regions$Type)& grepl("FRAME",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("NON",AncDesc_regions$Type)& grepl("STOP",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("NON",AncDesc_regions$Type)& grepl("START",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("NON",AncDesc_regions$Type)& grepl("UTR",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("SYN",AncDesc_regions$Type)& grepl("FRAME",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("SYN",AncDesc_regions$Type)& grepl("STOP",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("SYN",AncDesc_regions$Type)& grepl("START",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("SYN",AncDesc_regions$Type)& grepl("UTR",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("FRAME",AncDesc_regions$Type)& grepl("STOP",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("FRAME",AncDesc_regions$Type)& grepl("START",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("FRAME",AncDesc_regions$Type)& grepl("UTR",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("START",AncDesc_regions$Type)& grepl("STOP",AncDesc_regions$Type)),])
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("START",AncDesc_regions$Type)& grepl("UTR",AncDesc_regions$Type)),]) 
AncDesc_regions0<-rbind(AncDesc_regions0,AncDesc_regions[which(grepl("STOP",AncDesc_regions$Type)& grepl("UTR",AncDesc_regions$Type)),])

###Ancestral,Parallel,protein coding regions
AncParallel_genes.coding<-parallel_genes.coding[c("genename.x","Anc")]
colnames(AncParallel_genes.coding)<-c("genename","Anc")
AncDesc_gene0<-merge(AncParallel_genes.coding,AncDesc_gene0,by="genename")
###Ancestral,Parallel,protein coding regions related to DNA
AncDesc_DNAgene<-AncDesc_gene0[grep("DNA",AncDesc_gene0$annotation),]
#AncDesc_DNAgene<-AncDesc_DNAgene[c("genename","Type","Isolate","events","padj","regulation","annotation","profile" )]

###Check if Isolates with SNPs match RNAseq profile
#AncDesc_DNAgene<-rbind(AncDesc_DNAgene[grepl("WT",AncDesc_DNAgene$Isolate) & grepl("WT",AncDesc_DNAgene$profile),],
#                    AncDesc_DNAgene[grepl("PKS",AncDesc_DNAgene$Isolate) & grepl("PKS",AncDesc_DNAgene$profile),])

write.csv(unique(AncDesc_regions),"Table_ancSNPparallel_aroundcoding.csv") ###1000s of genes
write.csv(unique(AncDesc_DNAgene),"Table_ancSNPparallel_RNA_DNArepair.csv")
AncDesc_gene0<-AncDesc_gene0[c("genename","Type","Isolate","Anc","events","padj","regulation","annotation","profile" )]
write.csv(unique(AncDesc_gene0),"Table1_ancSNPparallel_RNA.csv")

############################################# Table S2 ################################################
####Parallel SNPs at the gene level that occur on different ancestral lineages of evolved strains #####
wt15_wt15<-parallel_genes.coding[grep("^(?=.*WT15.*WT15.*)(?!.*PKS15.*)",parallel_genes.coding$Isolate,perl=TRUE),]
pks15_pks15<-parallel_genes.coding[grep("^(?=.*PKS15.*PKS15.*)(?!.*WT15.*)",parallel_genes.coding$Isolate,perl=TRUE),]
wt15_pks15<-parallel_genes.coding[grep("^(?=.*WT15.*PKS15.*)(?=.*PKS15.*WT15.*)",parallel_genes.coding$Isolate,perl=TRUE),]
write.csv(rbind(wt15_wt15,pks15_pks15,wt15_pks15),"TableS4_ParallelSNPsGeneLevel_codingregion.csv")
