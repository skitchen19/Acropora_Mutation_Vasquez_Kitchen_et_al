########################
### Bar Plot of SMs ###
########################
library(ggplot2)
library(reshape2)

GOH <- c(217,86)
LOH <- c(50,47)

tab<-rbind(GOH,LOH)

tab2<-melt(tab)

colnames(tab)<-c("SM","Inherited SM")

ggplot(tab2, aes(Var1, value, fill=as.factor(Var2)))+
  geom_bar(stat="identity",position=position_dodge(0.65), color="white", width=0.6) +
  theme_classic()+
  ylab("Number of SMs") +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values=c("#004166","#57ABDB"))

##################
### Import VCF ###
##################

#Required R-packages for multi-locus genotype calling
library(data.table)
library(dplyr)
library(tidyr)
library(vcfR)
library(cowplot)
library(zoo)
library(ggpubr)

# Read in VCF input file.
#vcf <- read.vcfR("E:/PSU/NOAA/PRO100175_PSU175_SAX_b04/PRO100175_PSU175_SAX_b04/dup_removed_DB_snps_P1-5.vcf")

setwd("C:/Users/Sheila's Comp/Desktop")
vcf <- read.vcfR("STAGdb_04Jan2022.vcf")

# Read in population information
population_info_data_table <- read.table("popInfo_STAGdb_04Jan2022.txt",
                                         check.names=FALSE, header=F, 
                                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region")

# Read in mutation probe IDs
mut<-read.table("E:/PSU/NOAA/somatic_mutation/smutations.txt", sep="\t", header=T)
head(mut)

# number of mutations per A.digitifera scaffold
cnt<- mut %>% group_by(CHROM) %>% count()

#########################################
### Extract BAF and LRR from VCF file ###
#########################################

#extract metadata from VCF
fix<-as.data.frame(getFIX(vcf))

#B-allele frequencies
baf <- extract.gt(vcf, element = 'BAF', as.numeric = TRUE)
baf2<-cbind(fix[,1:2],baf)

#reformat BAF
dfm <- reshape2::melt(baf2,id.vars = c("CHROM","POS"))
dfm$POS<-as.numeric(paste(dfm$POS))
head(dfm)

#logR ratios
lrr<- extract.gt(vcf, element = 'LRR', as.numeric = TRUE)
lrr2<-cbind(fix[,1:2],lrr)

#reformat LRR
lfm <- reshape2::melt(lrr2,id.vars = c("CHROM","POS"))
lfm$POS<-as.numeric(paste(lfm$POS))
head(lfm)

#####################################################
### Subset dataframes for SELF and Parent samples ###
#####################################################

parent<-c("CHROM", "POS","a550962-4393310-052921-062_A01.CEL",
          "a550962-4393310-052921-062_C01.CEL",
          "a550962-4393310-052921-062_E01.CEL",
          "a550962-4393310-052921-062_G01.CEL",
          "a550962-4393310-052921-062_I01.CEL",
          "a550962-4393310-052921-062_K01.CEL",
          "a550962-4393310-052921-062_M01.CEL",
          "a550962-4393310-052921-062_O01.CEL",
          "a550962-4393310-052921-062_A03.CEL",
          "a550962-4393310-052921-062_C03.CEL",
          "a550962-4393310-052921-062_E03.CEL",
          "a550962-4393310-052921-062_G03.CEL",
          "a550962-4393310-052921-062_I03.CEL",
          "a550962-4393310-052921-062_K03.CEL",
          "a550962-4393310-052921-062_M03.CEL"
)

self<-c("CHROM", "POS","a550962-4393310-052921-062_O03.CEL","a550962-4393310-052921-062_E05.CEL",
        "a550962-4393310-052921-062_I05.CEL","a550962-4393310-052921-062_K05.CEL",
        "a550962-4393310-052921-062_M05.CEL","a550962-4393310-052921-062_O05.CEL",
        "a550962-4393310-052921-062_A07.CEL","a550962-4393310-052921-062_G07.CEL",
        "a550962-4393310-052921-062_I07.CEL","a550962-4393310-052921-062_C09.CEL",
        "a550962-4393310-052921-062_I09.CEL","a550962-4393310-052921-062_K09.CEL",
        "a550962-4393310-052921-062_M09.CEL","a550962-4393310-052921-062_O09.CEL",
        "a550962-4393310-052921-062_A11.CEL","a550962-4393310-052921-062_C11.CEL",
        "a550962-4393310-052921-062_E11.CEL","a550962-4393310-052921-062_G11.CEL",
        "a550962-4393310-052921-062_I11.CEL","a550962-4393310-052921-062_K11.CEL",
        "a550962-4393310-052921-062_M11.CEL","a550962-4393310-052921-062_O11.CEL",
        "a550962-4393310-052921-062_A13.CEL","a550962-4393310-052921-062_C13.CEL",
        "a550962-4393310-052921-062_E13.CEL","a550962-4393310-052921-062_G13.CEL",
        "a550962-4393310-052921-062_I13.CEL","a550962-4393310-052921-062_K13.CEL",
        "a550962-4393310-052921-062_M13.CEL","a550962-4393310-052921-062_O13.CEL")

# subset offspring
lrr3<-lrr2[,colnames(lrr2) %in% self ]
lfm3 <- melt(lrr3,id.vars = c("CHROM","POS"))
lfm3$POS<-as.numeric(paste(lfm3$POS))
lfm3$group<-"offspring"
head(lfm3)

baf3<-baf2[,colnames(baf2) %in% self ]
dfm3 <- melt(baf3,id.vars = c("CHROM","POS"))
dfm3$POS<-as.numeric(paste(dfm3$POS))
dfm3$group<-"offspring"
head(dfm3)

#subset parent ramets
baf4<-baf2[,colnames(baf2) %in% parent ]
dfm4 <-gather(baf4,key, value, -CHROM, -POS)
dfm4$POS<-as.numeric(paste(dfm4$POS))
dfm4$group<-"parent"
colnames(dfm4)<-c("CHROM","POS","variable","value","group")
str(dfm4)
head(dfm4)

lrr4<-lrr2[,colnames(lrr2) %in% parent ]
lfm4 <- melt(lrr4,id.vars = c("CHROM","POS"))
lfm4$POS<-as.numeric(paste(lfm4$POS))
lfm4$group<-"parent"
head(lfm4)

###########################################
### Combine offspring and parent tables ###
###########################################

#combine LRR tables
lfm4<-rbind(lfm3,lfm4)

lfm5<-lfm4 %>% 
  left_join(population_info_data_table %>% select("user_specimen_id","affy_id"), by=c("variable"="affy_id")) %>%
  left_join(mut %>% select("probe","POS"), by=c("POS")) %>%
  mutate(type="LRR")

#combine BAF tables
dfm4<-rbind(dfm3,dfm4)

dfm5<-dfm4 %>% 
  left_join(population_info_data_table %>% dplyr::select("user_specimen_id","affy_id"), by=c("variable"="affy_id"))%>%
  left_join(mut %>% dplyr::select("probe","POS"), by=c("POS"))%>%
  mutate(type="BAF")

#combine BAF and LRRs
lfm6<-rbind(lfm5,dfm5)

###################################################
### Extract A. digitifera Scaffolds with SMs ###
###################################################

#loop through chroms with SMs
chr<-lfm6 %>% drop_na(probe) %>% select(CHROM) %>% unique()
chr<-unlist(chr)

#############################
### Plotting BAF/LRR loop ###
#############################

for(i in chr){
  print(i)
  
  lfm7<-lfm6 %>% drop_na(probe) %>% subset(CHROM==i)
  lfm11 <-lfm6 %>% subset(CHROM==i)
  
  parent_samples<-c("a550962-4393310-052921-062_A01.CEL",
             "a550962-4393310-052921-062_C01.CEL",
             "a550962-4393310-052921-062_E01.CEL",
             "a550962-4393310-052921-062_G01.CEL",
             "a550962-4393310-052921-062_I01.CEL",
             "a550962-4393310-052921-062_K01.CEL",
             "a550962-4393310-052921-062_M01.CEL",
             "a550962-4393310-052921-062_O01.CEL",
             "a550962-4393310-052921-062_A03.CEL",
             "a550962-4393310-052921-062_C03.CEL",
             "a550962-4393310-052921-062_E03.CEL",
             "a550962-4393310-052921-062_G03.CEL",
             "a550962-4393310-052921-062_I03.CEL",
             "a550962-4393310-052921-062_K03.CEL",
             "a550962-4393310-052921-062_M03.CEL")
  
  lfm8<-lfm7[lfm7$variable %in% parent_samples,]
  lfm12<-lfm11[lfm11$variable %in% parent_samples,]

  p <- ggplot(lfm12,aes(x = POS/1e6, y = value)) 
  p <- p + geom_point(lfm12 %>% subset(type == "BAF"),
                    mapping=aes(x = POS/1e6, y = value),color=alpha("grey",0.5), pch=16, size=0.8) 
  p <- p + geom_point(lfm12 %>% subset(type == "LRR"),
                    mapping=aes(x = POS/1e6, y = value),color=alpha("grey",0.5), pch=16,size=0.8, inherit.aes=F)
  p <- p + geom_hline(yintercept = 0, color="grey")
  p <- p + geom_hline(yintercept = 1, color="grey")
  p <- p + geom_point(lfm8,mapping=aes(x = POS/1e6, y = value, color=probe, shape=probe, fill=probe), size=1.5) +
    facet_grid(type ~ user_specimen_id +group, scales='free_y')+
    #scale_color_manual(values=c("#93CEF1","#57ABDB","#157BB7","#004166", "#57ABDB"))+
    scale_color_manual(values=c("black","black","black","black", "black"))+
    scale_fill_manual(values=c("black","black","black","black", "black"))+
    scale_shape_manual(values=c(25,16,17,18,15))
  p <- p + scale_x_continuous(paste(i," position (Mb)"),breaks=seq(0, max(lfm6$POS)/1e6, 0.2)) +
    scale_y_continuous(paste())+
    theme_classic() +
    theme(legend.position = 'none', strip.background = element_rect(fill="grey95",colour = NA),
        panel.spacing.y = unit(0.3, "lines"),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 6))+
    coord_cartesian(expand = TRUE) + panel_border(size = 0.3)
  
  #offspring first 15 samples
  off1_samples<-c("a550962-4393310-052921-062_O03.CEL","a550962-4393310-052921-062_E05.CEL",
                  "a550962-4393310-052921-062_I05.CEL","a550962-4393310-052921-062_K05.CEL",
                  "a550962-4393310-052921-062_M05.CEL","a550962-4393310-052921-062_O05.CEL",
                  "a550962-4393310-052921-062_A07.CEL","a550962-4393310-052921-062_G07.CEL",
                  "a550962-4393310-052921-062_I07.CEL","a550962-4393310-052921-062_C09.CEL",
                  "a550962-4393310-052921-062_I09.CEL","a550962-4393310-052921-062_K09.CEL",
                  "a550962-4393310-052921-062_M09.CEL","a550962-4393310-052921-062_O09.CEL",
                  "a550962-4393310-052921-062_A11.CEL")
  
  lfm9<-lfm7[lfm7$variable %in% off1_samples,]
  lfm13<-lfm11[lfm11$variable %in% off1_samples,]
  
  p2 <- ggplot(lfm13,aes(x = POS/1e6, y = value)) 
  p2 <- p2 + geom_point(lfm13 %>% subset(type == "BAF"),
                      mapping=aes(x = POS/1e6, y = value),color=alpha("grey",0.5), pch=16, size=0.8) 
  p2 <- p2 + geom_point(lfm13 %>% subset(type == "LRR"),
                      mapping=aes(x = POS/1e6, y = value),color=alpha("grey",0.5), pch=16,size=0.8, inherit.aes=F)
  p2 <- p2 + geom_hline(yintercept = 0, color="grey")
  p2 <- p2 + geom_hline(yintercept = 1, color="grey")
  p2 <- p2 + geom_point(lfm9, mapping=aes(x = POS/1e6, y = value, color=probe, shape=probe, fill=probe), size=1.5) +
    facet_grid(type ~ user_specimen_id +group, scales='free_y')+
    # scale_color_manual(values=c("#93CEF1","#57ABDB","#157BB7","#004166", "#57ABDB"))+
    scale_color_manual(values=c("black","black","black","black", "black"))+
    scale_fill_manual(values=c("black","black","black","black", "black"))+
    scale_shape_manual(values=c(25,16,17,18,15))
  p2 <- p2 + scale_x_continuous(paste(i," position (Mb)"),breaks=seq(0, max(lfm6$POS)/1e6, 0.2)) +
    scale_y_continuous(paste())+
    theme_classic() +
    theme(legend.position = 'none', strip.background = element_rect(fill="grey95",colour = NA),
          panel.spacing.y = unit(0.3, "lines"),
          strip.text.x = element_text(size = 6),
          axis.text = element_text(size = 6))+
    coord_cartesian(expand = TRUE) + panel_border(size = 0.3)
  
  #offspring second 15 samples
  off2_samples<-c("a550962-4393310-052921-062_C11.CEL",
                  "a550962-4393310-052921-062_E11.CEL","a550962-4393310-052921-062_G11.CEL",
                  "a550962-4393310-052921-062_I11.CEL","a550962-4393310-052921-062_K11.CEL",
                  "a550962-4393310-052921-062_M11.CEL","a550962-4393310-052921-062_O11.CEL",
                  "a550962-4393310-052921-062_A13.CEL","a550962-4393310-052921-062_C13.CEL",
                  "a550962-4393310-052921-062_E13.CEL","a550962-4393310-052921-062_G13.CEL",
                  "a550962-4393310-052921-062_I13.CEL","a550962-4393310-052921-062_K13.CEL",
                  "a550962-4393310-052921-062_M13.CEL","a550962-4393310-052921-062_O13.CEL")
  
  lfm10 <-lfm7[lfm7$variable %in% off2_samples,]
  lfm14<-lfm11[lfm11$variable %in% off2_samples,]
  
  p3 <- ggplot(lfm14,aes(x = POS/1e6, y = value)) 
  p3 <- p3 + geom_point(lfm14 %>% subset(type == "BAF"),
                        mapping=aes(x = POS/1e6, y = value),color=alpha("grey",0.5), pch=16, size=0.8) 
  p3 <- p3 + geom_point(lfm14 %>% subset(type == "LRR"),
                        mapping=aes(x = POS/1e6, y = value),color=alpha("grey",0.5), pch=16,size=0.8, inherit.aes=F)
  p3 <- p3 + geom_hline(yintercept = 0, color="grey")
  p3 <- p3 + geom_hline(yintercept = 1, color="grey")
  p3 <- p3 + geom_point(lfm10, mapping=aes(x = POS/1e6, y = value, color=probe, shape=probe, fill=probe), size=1.5) +
    facet_grid(type ~ user_specimen_id +group, scales='free_y')+
    # scale_color_manual(values=c("#93CEF1","#57ABDB","#157BB7","#004166", "#57ABDB"))+
    scale_color_manual(values=c("black","black","black","black", "black"))+
    scale_fill_manual(values=c("black","black","black","black", "black"))+
    scale_shape_manual(values=c(25,16,17,18,15))
  p3 <- p3 + scale_x_continuous(paste(i," position (Mb)"),breaks=seq(0, max(lfm6$POS)/1e6, 0.2)) +
    scale_y_continuous(paste())+
    theme_classic() +
    theme(legend.position = 'bottom', legend.box = 'horizontal', strip.background = element_rect(fill="grey95",colour = NA),
          panel.spacing.y = unit(0.3, "lines"),
          strip.text.x = element_text(size = 6),
          axis.text = element_text(size = 6)) +
    coord_cartesian(expand = TRUE) + panel_border(size = 0.3)
  
  pdf(paste0("E:/PSU/NOAA/somatic_mutation/b-allele_freq/SM_plots_all_samples/",i,".pdf"), paper = "letter", width = 8.5, height = 11)
  print(ggarrange(p + rremove("x.title"),p2 + rremove("x.title"),p3,nrow = 3, align="v"))
  dev.off()
}

######################
## Import Mixes VCF ##
######################

# Read in VCF input file.
vcf <- read.vcfR("E:/PSU/NOAA/somatic_mutation/figures/OneDrive_1_5-16-2021/Mixed_Samples.vcf")

#########################################
### Extract BAF and LRR from VCF file ###
#########################################

#extract metadata from VCF
fix<-as.data.frame(getFIX(vcf))

#B-allele frequencies
baf <- extract.gt(vcf, element = 'BAF', as.numeric = TRUE,mask = FALSE)
baf2<-cbind(fix[,1:2],baf)

#reformat BAF
dfm <- reshape2::melt(baf2,id.vars = c("CHROM","POS"))
dfm$POS<-as.numeric(paste(dfm$POS))
head(dfm)

dfm2<-dfm %>% mutate(type="BAF") %>% group_by(variable) %>% mutate(id=row_number())

#logR ratios
lrr<- extract.gt(vcf, element = 'LRR', as.numeric = TRUE)
lrr2<-cbind(fix[,1:2],lrr)

#reformat LRR
lfm <- reshape2::melt(lrr2,id.vars = c("CHROM","POS"))
lfm$POS<-as.numeric(paste(lfm$POS))
head(lfm)

lfm2<-lfm %>%  mutate(type="LRR") %>% group_by(variable) %>% mutate(id=row_number())

#combine BAF and LRRs
df<-rbind(lfm2,dfm2)

###################################################
### Extract A. digitifera Scaffolds with SMs ###
###################################################

#loop through chroms with SMs
sample<-lfm %>% select(variable) %>% unique()
sample<-unlist(sample)

#############################
### Plotting BAF/LRR loop ###
#############################
library("ggExtra")

for(i in sample){
  print(i)
  
  lfm3<-df %>% subset(variable==i)
  
  p <- ggplot(lfm3,aes(x = id, y = value)) 
  p <- p + geom_point(lfm3 %>% subset(type == "BAF"),
                      mapping=aes(x = id, y = value,color = cut(value, c(Inf,0.95,0.55,0.45,0.05,-Inf))), pch=16, size=0.3) 
  p <- p + geom_point(lfm3 %>% subset(type == "LRR"),
                      mapping=aes(x = id, y = value),color=alpha("darkgrey",0.8), pch=16,size=0.3, inherit.aes=F)
  p <- p + geom_hline(yintercept = 0, color="black")+
    facet_grid(type ~ ., scales='free_y')+
    scale_color_manual(values=c("black","grey80","black","grey80", "black"))
  p <- p + scale_x_discrete(paste("Ordered Position"),breaks=seq(0, max(lfm3$id), 500)) +
    scale_y_continuous(paste())+
    theme_classic() +
    theme(legend.position = 'none', strip.background = element_rect(fill="grey95",colour = NA),
          panel.spacing.y = unit(0.3, "lines"),
          strip.text.x = element_text(size = 6),
          axis.text = element_text(size = 6))+
    coord_cartesian(expand = TRUE) + panel_border(size = 0.3)
  
  pdf(paste0("E:/PSU/NOAA/somatic_mutation/b-allele_freq/dna_mixes/",i,".pdf"), paper = "letter", width = 4.25, height = 2.75)
  print(p)
  dev.off()
}


ggplot(df %>% subset(type=="LRR"), aes(x=value,fill=variable)) + 
  geom_histogram(alpha=.5, bins=30,position="identity") +facet_wrap(~variable) +
  xlim(-1.5,1.5)+ theme(legend.position = 'none')+
  scale_fill_manual(values=c("#8D4C9E","#8D4C9E","#F79028","#F79028","black","black"))

ggplot(df %>% subset(type=="LRR"), aes(x=value,y=variable, fill=variable)) + 
  geom_violin(alpha=.5)+theme_bw()+theme(legend.position = 'none')+xlab("LRR")+
  scale_fill_manual(values=c("#8D4C9E","#8D4C9E","#F79028","#F79028","black","black"))

ggplot(df %>% subset(type=="BAF"), aes(x=value,fill=variable)) + 
  geom_histogram(alpha=.5, bins=50,position="identity") +facet_wrap(~variable)+theme_bw()+
  scale_y_continuous(expand=c(0,0), limits=c(0,7000)) +scale_x_continuous(expand=c(0.01,0.01)) +xlab("BAF")+
  theme(legend.position = 'none')+ geom_hline(yintercept = 5000, color="black") + 
  scale_fill_manual(values=c("#8D4C9E","#8D4C9E","#F79028","#F79028","black","black"))

######################################################
## Calculate number of alleles in each range of BAF ##
######################################################

df2<-df %>% subset(type=="BAF") %>% group_by(gr=cut(value, c(Inf,0.99,0.55,0.45,0.01,-Inf)),variable,type) %>% 
  summarise(n= n()) %>%
  arrange(as.numeric(gr))

df3<-df %>% subset(type=="LRR") %>% group_by(gr=cut(value, c(Inf,1,0.05,-0.05,-1,-Inf)),variable,type) %>% 
  summarise(n= n()) %>%
  arrange(as.numeric(gr))
