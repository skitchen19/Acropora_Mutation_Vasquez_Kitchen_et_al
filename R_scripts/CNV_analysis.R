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
vcf <- read.vcfR("E:/PSU/NOAA/PRO100175_PSU175_SAX_b04/PRO100175_PSU175_SAX_b04/dup_removed_DB_snps_P1-5.vcf")

# Read in population information
population_info_data_table <- read.table("E:/PSU/NOAA/PRO100175_PSU175_SAX_b04/popInfo_P1-5.txt",
                                         check.names=FALSE, header=F, 
                                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region")

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


####################
### Sample lists ###
####################

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

############################
### copy Number Analysis ###
############################

# load in R packages
library(DNAcopy)
library(stringr)
library(CGHcall)

# start with offspring
sub<-lfm5 %>% select(CHROM, POS,variable, value) %>%
  subset(variable %in% self) %>%
  spread(.,variable,value) %>%
  mutate(id=rownames(.), POS_end=POS+1) %>%
  relocate(id,.before=CHROM) %>%
  relocate(POS_end,.before="a550962-4393310-052921-062_A07.CEL")

sub$CHROM<-str_replace_all(sub$CHROM, "[a-zA-Z_]", " ")
sub$CHROM<-str_replace_all(sub$CHROM, "\\.1", " ")
head(sub)

#CGHcall program on offspring
cgh <- make_cghRaw(sub)

norm.cghdata <- normalize(cgh, method="median", smoothOutliers=TRUE)

seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy",undo.splits="sdundo",undo.SD=3,
                           clen=10, relSDlong=5)

postseg.cghdata <- postsegnormalize(seg.cghdata)
seg<-postseg.cghdata@assayData[["segmented"]]

#write.table(seg,"Segmented_offspring.txt",sep="\t", quote=F)

result <- CGHcall(postseg.cghdata,nclass=5,robustsig="yes", organism = "other")

result[,1]
result <- ExpandCGHcall(result,postseg.cghdata)

gain<-result@assayData[["probgain"]]
write.table(gain,"CNV_gain_offspring.txt", sep="\t",quote=F)

loss<-result@assayData[["probloss"]]
write.table(loss,"CNV_loss_offspring2.txt", sep="\t",quote=F)

write.table(sub,"CNV_LRR_offspring.txt", sep="\t",quote=F)

# now subset parentss
sub<-lfm5 %>% select(CHROM, POS,variable, value) %>%
  subset(variable %in% parent) %>%
  spread(.,variable,value) %>%
  mutate(id=rownames(.), POS_end=POS+1) %>%
  relocate(id,.before=CHROM) %>%
  relocate(POS_end,.before="a550962-4393310-052921-062_A01.CEL")

sub$CHROM<-str_replace_all(sub$CHROM, "[a-zA-Z_]", " ")
sub$CHROM<-str_replace_all(sub$CHROM, "\\.1", " ")
head(sub)

#CGHcall program on parents
cgh <- make_cghRaw(sub)

norm.cghdata <- normalize(cgh, method="median", smoothOutliers=TRUE)

seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy",undo.splits="sdundo",undo.SD=3,
                           clen=10, relSDlong=5)

postseg.cghdata <- postsegnormalize(seg.cghdata)
seg<-postseg.cghdata@assayData[["segmented"]]

#write.table(seg,"Segmented_parent.txt",sep="\t", quote=F)

result <- CGHcall(postseg.cghdata,nclass=5,robustsig="yes", organism = "other")

result <- ExpandCGHcall(result,postseg.cghdata)

gain<-result@assayData[["probgain"]]
write.table(gain,"CNV_gain_parent.txt", sep="\t",quote=F)

loss<-result@assayData[["probloss"]]
write.table(loss,"CNV_loss_parent.txt", sep="\t",quote=F)

write.table(sub,"CNV_LRR_parent.txt", sep="\t",quote=F)


##################################################
### Stats on CNV gains and losses ###
##################################################

library(chisq.posthoc.test)

# predicted losses between parent and offspring
dat_loss <- data.frame(Parent = c(1254,15,9), 
                       Offspring = c(1840,27,16), 
                       row.names = c("no change", "PEM", "PEM_inherited"))
dat_loss

chi.t <- chisq.test(dat_loss)
fisher.test(dat_loss)
mosaicplot(chi.t$observed, cex.axis =1 , main = "Observed counts")
chisq.posthoc.test(dat_loss)

# predicted gains between parent and offspring
dat_gain <- data.frame(Parent = c(556,6,4), 
                       Offspring = c(2442,37,20), 
                       row.names = c("no change", "PEM", "PEM_inherited"))
dat_gain

chi.t <- chisq.test(dat_gain)
mosaicplot(chi.t$observed, cex.axis =1 , main = "Observed counts")
chisq.posthoc.test(dat_gain)

# differences between parent and offspring gains and losses
dat_gl <-data.frame(loss = c(1867,1269), 
                    gains = c(2479,562), 
                    row.names = c("offspring", "parent"))

chi.t <- chisq.test(dat_gl)
mosaicplot(chi.t$observed, cex.axis =1 , main = "Observed counts")
chisq.posthoc.test(dat_gl)

chi.t$observed
chi.t$expected

# parent only gains vs. losses
dat_parent <- data.frame(loss = c(1254,15,9), 
                         gain = c(556,6,4), 
                         row.names = c("no change", "PEM", "PEM_inherited"))
chi.t <- chisq.test(dat_parent)
chisq.posthoc.test(dat_parent)

# offspring only gains vs. losses
dat_off <- data.frame(loss = c(1840,27,16), 
                      gain = c(2442,37,20), 
                      row.names = c("no change", "PEM", "PEM_inherited"))
chi.t <- chisq.test(dat_off)
chisq.posthoc.test(dat_off)

###########################################
### Plot CNV density by scaffold length ###
###########################################

# load in files of scaffold size and predicted CNVs, sorted in Excel from output above
scaff_len<-read.table("E:/PSU/NOAA/somatic_mutation/b-allele_freq/CNV/Adig_scaff_number_loss.txt",header=T)
scaff_len2<-read.table("E:/PSU/NOAA/somatic_mutation/b-allele_freq/CNV/Adig_scaff_number_gain.txt",header=T)

#plot density of gains and losses
ggplot(scaff_len, aes(x = Size/1e6)) +
  theme_classic()+ xlim(0,2)+
  geom_density(alpha=0.5, fill='blue')+
  geom_density(data=scaff_len2, aes(x = Size/1e6), alpha=0.5, fill="orange")

# load in file that includes in silico and  visual inspection, sorted in Excel from output above and visual inspection
scaff_len3<-read.table("E:/PSU/NOAA/somatic_mutation/b-allele_freq/CNV/Adig_scaff_CNV_visuals.txt",header=T, sep="\t")

p<-ggplot(scaff_len3, aes(x = Size/1e6)) +
  theme_classic()+ xlim(0,3)+
  geom_density(data=scaff_len3[, c(1:4,7)] %>% subset(InSilico=="1"), aes(x = Size/1e6), alpha=0.5, fill="pink")

p2<-ggplot(scaff_len3, aes(x = Size/1e6)) +
  theme_classic()+ xlim(0,3)+
  geom_density(data=scaff_len3[, c(1:4,8)] %>% subset(Visual=="1"), aes(x = Size/1e6), alpha=0.5, fill="orange")

ggarrange(p + rremove("x.title"),p2,nrow = 2, align="v")

# identify peak length of each density plot
which.max(density(scaff_len$Size)$y)
density(scaff_len$Size)$x[58]

which.max(density(scaff_len2$Size)$y)
density(scaff_len2$Size)$x[160]

which.max(density(scaff_len3$Size)$y)
density(scaff_len$Size)$x[58]

l<-scaff_len3 %>%
  subset(InSilico=="1") 

which.max(density(l$Size)$y)
density(l$Size)$x[131]

l<-scaff_len3 %>%
  subset(Visual=="1") 

which.max(density(l$Size)$y)
density(l$Size)$x[156]

# calculate median size with predicted CNV
scaff_len %>%
  filter(loss >= 0) %>%
  summarize(median(Size))

scaff_len2 %>%
  filter(gain >= 0) %>%
  summarize(median(Size))

# statistical tests for gains and losses of predicted CNVs
m1<-lm(l_s ~ Size, data=scaff_len)
summary(m1)
plot(m1$residuals, pch = 16, col = "red")

m2<-lm(g_s ~ Size, data=scaff_len)
summary(m2)
plot(m2$residuals, pch = 16, col = "red")

##################################
### UpSetR Plot of shared CNVs ###
##################################

#install.packages("UpSetR")
library("UpSetR")

Input <- read.table("CNV_freq_PA.txt",header=T)
head(Input)
Input$PEM<-as.factor(Input$PEM)

upset(
  Input, 
  sets=c("gain_frequency_offspring","loss_frequency_offspring",
         "gain_frequency_parent","loss_frequency_parent"),
  sets.x.label = "Number of P(CNV) > 0.5", 
  mainbar.y.label = "Intersection of CNV gains/losses",
  mb.ratio = c(0.7, 0.3),
  order.by="freq", 
  keep.order = T,
  line.size = 0.7,
  point.size = 3,
  number.angles = 25,
  text.scale = c(2, 1.3, 1, 1, 1.5, 1.3),
  query.legend = "top",
  sets.bar.color = c("#D55E00","#0072B2","#D55E00","#0072B2"),
  queries = list(
    list(query = elements, 
         params = list("PEM", c("Normal","PEM","PEM_inherited")),
         color = "grey", active = T,query.name="No change"),
    list(query = elements, 
         params = list("PEM", c("PEM","PEM_inherited")),
         color = "#004166", active = T,query.name="PEM"),
    list(query = elements, 
         params = list("PEM", "PEM_inherited"),
         color = "#57ABDB", active = T,query.name="PEM_inherited")))


Input2 <- read.table("CNV_freq_PA_PEM.txt",header=T)
head(Input2)
Input2$PEM<-as.factor(Input2$PEM)

upset(
  Input2, 
  sets=c("gain_frequency_offspring","loss_frequency_offspring",
         "gain_frequency_parent","loss_frequency_parent"),
  sets.x.label = "Number of P(CNV) > 0.5", 
  mainbar.y.label = "Number of shared CNV gains/losses",
  mb.ratio = c(0.8, 0.2),
  order.by="freq", 
  keep.order = T,
  line.size = 0.7,
  point.size = 3,
  number.angles = 25,
  text.scale = c(2, 1.3, 1, 1, 1.5, 1.3),
  sets.bar.color = c("#D55E00","#0072B2","#D55E00","#0072B2"),
  queries = list(
    list(query = elements, 
         params = list("PEM", c("PEM","PEM_inherited")),
         color = "#004166", active = T),
    list(query = elements, 
         params = list("PEM", c("PEM_inherited")),
         color = "#57ABDB", active = T))
)

###################
## Plot heatmap ###
###################

library(ComplexHeatmap)

myCol<-c("#0072B2","grey98","#D55E00")

#load in table of combined offspring and parent predictions from above
mat <- read.table("CNV_heatmap_input.txt",header=T,check.names=FALSE)

# different sets of samples of interests
parentColony<-c("a550962-4393310-052921-062_A01.CEL",
         "a550962-4393310-052921-062_C01.CEL",
         "a550962-4393310-052921-062_E01.CEL",
         "a550962-4393310-052921-062_G01.CEL",
         "a550962-4393310-052921-062_I01.CEL",
         "a550962-4393310-052921-062_K01.CEL",
         "a550962-4393310-052921-062_M01.CEL",
         "a550962-4393310-052921-062_O01.CEL",
         "a550962-4393310-052921-062_A03.CEL",
         "a550962-4393310-052921-062_C03.CEL")

neighborRamet<-c(
  "a550962-4393310-052921-062_E03.CEL",
  "a550962-4393310-052921-062_G03.CEL",
  "a550962-4393310-052921-062_I03.CEL",
  "a550962-4393310-052921-062_K03.CEL",
  "a550962-4393310-052921-062_M03.CEL")

automictic<-c("a550962-4393310-052921-062_K09.CEL",
              "a550962-4393310-052921-062_M13.CEL")

self<-c("a550962-4393310-052921-062_O03.CEL","a550962-4393310-052921-062_E05.CEL",
        "a550962-4393310-052921-062_I05.CEL","a550962-4393310-052921-062_K05.CEL",
        "a550962-4393310-052921-062_M05.CEL","a550962-4393310-052921-062_O05.CEL",
        "a550962-4393310-052921-062_A07.CEL","a550962-4393310-052921-062_G07.CEL",
        "a550962-4393310-052921-062_I07.CEL","a550962-4393310-052921-062_C09.CEL",
        "a550962-4393310-052921-062_I09.CEL",
        "a550962-4393310-052921-062_M09.CEL","a550962-4393310-052921-062_O09.CEL",
        "a550962-4393310-052921-062_A11.CEL","a550962-4393310-052921-062_C11.CEL",
        "a550962-4393310-052921-062_E11.CEL","a550962-4393310-052921-062_G11.CEL",
        "a550962-4393310-052921-062_I11.CEL","a550962-4393310-052921-062_K11.CEL",
        "a550962-4393310-052921-062_M11.CEL","a550962-4393310-052921-062_O11.CEL",
        "a550962-4393310-052921-062_A13.CEL","a550962-4393310-052921-062_C13.CEL",
        "a550962-4393310-052921-062_E13.CEL","a550962-4393310-052921-062_G13.CEL",
        "a550962-4393310-052921-062_I13.CEL","a550962-4393310-052921-062_K13.CEL",
        "a550962-4393310-052921-062_O13.CEL")

t2<- as.data.frame(t(mat[,8:52])) %>%
  mutate(Sample=row.names(as.data.frame(t(mat[,8:52])))) %>%
  mutate(group=ifelse(Sample %in% parentColony, "parentColony", 
                      ifelse(Sample %in% neighborRamet, "neighborRamet", 
                             ifelse(Sample %in% automictic, "automictic","selfed")))))


ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 6),
       heatmap_row_title_gp = gpar(fontsize = 6),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)

column_ha = HeatmapAnnotation(mutations = factor(mat$PEM, levels=c("Normal","PEM","PEM_inherited")),
                              col = list(mutations = c("Normal" = "grey", "PEM" = "#004166", "PEM_inherited" = "#57ABDB")),
                              annotation_legend_param = list(
                                mutations = list(direction = "horizontal")))

hm<-Heatmap(t(mat[,8:52]), name = "Parental haplotype", col=myCol, border = TRUE,
            row_split = factor(t2$group, levels=c("parentColony","neighborRamet","selfed", "automictic")),cluster_row_slices = FALSE, 
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, show_column_names =FALSE,row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = FALSE, 
            heatmap_legend_param = list(
              direction = "horizontal",
              labels = c("< 2", "2", "> 2")),
            row_title_rot = 0,row_gap = unit(1, "mm"),width = unit(15, "cm"), height = unit(8, "cm"),top_annotation = column_ha)
draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
