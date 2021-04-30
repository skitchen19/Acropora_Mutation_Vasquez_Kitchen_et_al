##################
### Import VCF ###
##################
library(vcfR)

vcf <- read.vcfR("E:/PSU/NOAA/PRO100175_PSU175_SAX_b04/PRO100175_PSU175_SAX_b04/dup_removed_DB_snps_P1-5.vcf")

##########################
### Convert to Genind  ###
##########################
library(poppr)

# Convert VCF file into a genind for the Poppr package.
genind_obj <- vcfR2genind(vcf)

# Add population information to the genind object.
population_info_data_table <- read.table("E:/PSU/NOAA/PRO100175_PSU175_SAX_b04/popInfo_P1-5.txt",
                                         check.names=FALSE, header=F, 
                                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
colnames(population_info_data_table) <- c("row_id", "affy_id", "user_specimen_id", "region")

genind_obj@pop <- as.factor(population_info_data_table$region)
strata(genind_obj) <- data.frame(pop(genind_obj))

# Convert genind object to a genclone object.
genind_clone <- as.genclone(genind_obj)

###########################################
### Calculate Absolute Genetic Distance ###
###########################################
# Calculate the bitwise distance between individuals.
bitwise_distance <- bitwise.dist(genind_clone)

genDist<-as.matrix(bitwise_distance)
#write.table(genDist,"geneticDistance2.txt",sep="\t")

#########################################################
## Histogram of Genetic Distance of A. palmata samples ##
## Kitchen et al. 2020, Scientific Reports             ##
#########################################################

#convert names of columns to match rownames of matrix
genDist2<-rownames_to_column(as.data.frame(genDist),"affy_id")

# read in list of A. palmata samples
AP_list<-read.table("aplam_list.txt",sep="\t", header=F,stringsAsFactors=F)

#subset matrix to only A. palmata samples in the list
gendist_AP<-genDist2[,colnames(genDist2) %in% unlist(AP_list)]
gendist_AP2<-genDist2[rownames((genDist2) %in% unlist(AP_list),]

# get upper genetic distance estimate
upper <- round(max(gendist_AP2), digits = 1)

#plot histogram
hist(gendist_AP2[gendist_AP2 > 0 & gendist_AP2 < 1],
     breaks=20,
     axes = TRUE,
     xlab = "Pairwise Genetic Distance",
     ylab = "Frequency of Pairs",
     xlim = c(0, upper), ylim=c(0,20000),main = "")

abline(v=0.032,lty=2, col="black", lwd=2) #genet threshold
abline(v=0.12773133,lty=1, col="#90C640", lwd=2) #Apalm genet avg

rect(0.005658987-0.005440774,0,0.005658987+0.005440774,20000, col= rgb(0,0,1.0,alpha=0.5),border = NA)#avg apalm ramet distance

#################################
### Subset Matrix for Samples ###
#################################
library(tibble)
library(dplyr)
library(tidyr)


# different sets of samples of interests
ramet<-c("a550962-4393310-052921-062_A01.CEL",
          "a550962-4393310-052921-062_C01.CEL",
          "a550962-4393310-052921-062_E01.CEL",
          "a550962-4393310-052921-062_G01.CEL",
          "a550962-4393310-052921-062_I01.CEL",
          "a550962-4393310-052921-062_K01.CEL",
          "a550962-4393310-052921-062_M01.CEL",
          "a550962-4393310-052921-062_O01.CEL",
          "a550962-4393310-052921-062_A03.CEL",
          "a550962-4393310-052921-062_C03.CEL")

rametNeighbor<-c(
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

cu_fl<-c("a550962-4393310-052921-062_E09.CEL",
         "a550962-4393310-052921-062_M07.CEL",
         "a550962-4393310-052921-062_O07.CEL",
         "a550962-4393310-052921-062_E07.CEL",
         "a550962-4393310-052921-062_C07.CEL",
         "a550962-4393310-052921-062_K07.CEL",
         "a550962-4393310-052921-062_A09.CEL",
         "a550962-4393310-052921-062_G09.CEL",
         "a550962-4393310-052921-062_G05.CEL",
         "a550962-4393310-052921-062_A05.CEL",
         "a550962-4393310-052921-062_C05.CEL")

samples<-c(ramet,rametNeighbor,automictic,self,cu_fl)

genDist3<-genDist2[genDist2$affy_id %in% samples,colnames(genDist2) %in% samples]

rownames(genDist3)<-colnames(genDist3)

genDist4<-genDist3 %>%
  mutate(affyID=rownames(.)) %>%
  mutate(group=ifelse(affyID %in% ramet, "ramet", 
                      ifelse(affyID %in% rametNeighbor, "rametNeighbor", 
                             ifelse(affyID %in% automictic, "automictic",
                                    ifelse(affyID %in% self, "selfed","CUxFL")))))
                  
#####################################
### Convert Matrix to Long Format ###
#####################################

library(reshape2)
meltGenDist<-melt(genDist4)

meltGenDist2<-meltGenDist %>%
  mutate(group2=ifelse(variable %in% ramet, "ramet", 
                      ifelse(variable %in% rametNeighbor, "rametNeighbor", 
                             ifelse(variable %in% automictic, "automictic",
                                    ifelse(variable %in% self, "selfed","CUxFL"))))) %>%
  mutate(group3=paste(group,"x",group2)) %>%
  filter(value !=0) %>%
  select(group3, value) %>%
  filter(group3 != "selfed x selfed") %>%
  filter(group3 != "automictic x automictic") %>%
  filter(group3 != "CUxFL x CUxFL") %>%
  filter(group3 != "rametNeighbor x rametNeighbor") %>%
  mutate(group3=replace(group3, group3 == "rametNeighbor x ramet", "ramet x rametNeighbor")) %>%
  mutate(group3=replace(group3, group3 == "rametNeighbor x selfed" , "parent x selfed")) %>%
  mutate(group3=replace(group3, group3 == "selfed x rametNeighbor" , "parent x selfed")) %>%
  mutate(group3=replace(group3, group3 == "ramet x selfed" , "parent x selfed")) %>%
  mutate(group3=replace(group3, group3 == "selfed x ramet" , "parent x selfed")) %>%
  mutate(group3=replace(group3, group3 == "rametNeighbor x automictic" , "parent x automictic")) %>%
  mutate(group3=replace(group3, group3 == "automictic x rametNeighbor" , "parent x automictic")) %>%
  mutate(group3=replace(group3, group3 == "ramet x automictic" , "parent x automictic")) %>%
  mutate(group3=replace(group3, group3 == "automictic x ramet" , "parent x automictic"))  %>%
  mutate(group3=replace(group3, group3 == "rametNeighbor x CUxFL" , "parent x CUxFL")) %>%
  mutate(group3=replace(group3, group3 == "CUxFL x rametNeighbor" , "parent x CUxFL")) %>%
  mutate(group3=replace(group3, group3 == "ramet x CUxFL" , "parent x CUxFL")) %>%
  mutate(group3=replace(group3, group3 == "CUxFL x ramet" , "parent x CUxFL")) %>%
  filter(group3 == "parent x CUxFL" | group3 == "parent x selfed"| group3 =="parent x automictic"| group3 == "ramet x ramet"
         | group3 == "ramet x rametNeighbor")

cu_cu<-read.csv("E:/PSU/NOAA/somatic_mutation/CURCUR_GenDIst.csv",
                check.names=FALSE, header=T, 
                na.strings=c("", "NA"), stringsAsFactors=FALSE,quote="")


cu_cu_filter<-cu_cu %>% select(Type_of_Relationship,dist)

colnames(cu_cu_filter)<-c("group3","value")

meltGenDist3<-rbind(meltGenDist2, cu_cu_filter)

##############################
### Plot Genetic Distances ###
##############################
library(ggplot2)
#library(ggbeeswarm)

#set names of factor levels
meltGenDist3$group3<-factor(meltGenDist3$group3, levels=c("ramet x ramet",
                                            "ramet x rametNeighbor",
                                            "parent x automictic",
                                            "parent x selfed",
                                            "parent x CUxFL",
                                            "parent x CUxCU"))

# plot genetic distance violin plot
p <- ggplot(meltGenDist3, aes(x=group3, y=value, color=group3, fill=group3)) + 
  geom_violin(lwd=0.8,trim = FALSE)+
  #geom_beeswarm(size=1.5, alpha=0.8,dodge.width = 0.3, cex=1.4)+
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  theme_classic() +
  ylab("Pairwise Genetic Distance")+
  xlab("")+
  theme(legend.position = "none")+
  scale_color_manual(values=c("black","grey","#E69F00","#D55E00","#0072B2", "#0072B2"))+
  scale_fill_manual(values=c("#00000020","#8F8F8F20","#E69F0020","#D55E0020","#0072B220","#0072B220"))+
  labs(color='Species') 
p 

#############################################################
### Summary Statistics of Genetic Distance for each Group ###
#############################################################

#Average and standard deviations
A <- mean(meltGenDist3$value[meltGenDist3$group3== "ramet x ramet"])
#0.004097279
As <- sd(meltGenDist3$value[meltGenDist3$group3== "ramet x ramet"])
#[1] 0.001712921
B <- mean(meltGenDist3$value[meltGenDist3$group3== "ramet x rametNeighbor"])
#0.005398558
Bs <- sd(meltGenDist3$value[meltGenDist3$group3== "ramet x rametNeighbor"])
#[1]0.002990064
C <- mean(meltGenDist3$value[meltGenDist3$group3== "parent x automictic"])
#0.02142652
Cs <- sd(meltGenDist3$value[meltGenDist3$group3== "parent x automictic"])
#[1]0.002351426
D <- mean(meltGenDist3$value[meltGenDist3$group3== "parent x selfed"])
#0.03713704
Ds <- sd(meltGenDist3$value[meltGenDist3$group3== "parent x selfed"])
#[1] 0.005153534
E <- mean(meltGenDist3$value[meltGenDist3$group3== "parent x CUxFL"])
#[1] 0.08248021
Es <- sd(meltGenDist3$value[meltGenDist3$group3== "parent x CUxFL"])
#[1] 0.002675624
F <- mean(meltGenDist3$value[meltGenDist3$group3== "parent x CUxCU"])
#[1] 0.06634868
Fs <- sd(meltGenDist3$value[meltGenDist3$group3== "parent x CUxCU"])
#[1] 0.003342121

# ANOVA tests on genetic distance by group
mod<-aov(meltGenDist3$value~meltGenDist3$group3)
summary(mod)

# post-hoc tests
tukey.test <- TukeyHSD(mod)

######################################################################
### Neighbor-joining Tree of Parent Colony and Neighboring Samples ###
######################################################################

samples2<-c(ramet,rametNeighbor)

genind_clone_Parents<-genind_clone[rownames(genind_clone@strata) %in% samples2,]

# Construct tree and plot.
set.seed(20210311)
nj_phylogeny_tree <- genind_clone_Parents  %>%
  aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=FALSE, showtree = FALSE, root=FALSE, mc.cores = 4)
nj_phylogeny_tree$tip.label <- population_info_data_table$user_specimen_id[match(nj_phylogeny_tree$tip.label, population_info_data_table$affy_id)]
plot.phylo(nj_phylogeny_tree, label.offset=0, type = "unrooted",
           cex=0.6, font=2, lwd=4, align.tip.label=F, no.margin=T)
# Add a scale bar showing 5% difference.
add.scale.bar(0, 0, length=0.001, cex=0.65, lwd=2)
nodelabels(nj_phylogeny_tree$node.label, cex=.5,  frame="n", font=3, xpd=TRUE)

#write.tree(nj_phylogeny_tree, file = "nj_parent_tree.nwk")

# Tree based on only the SMs
mut<-read.table("E:/PSU/NOAA/somatic_mutation/smutations.txt", sep="\t", header=T)

#subset to the SMs
genind_clone_Parents_mutOnly<-genind_clone_Parents[,genind_clone_Parents@loc.fac %in% mut$probe]

# Construct tree and plot.
set.seed(20210311)
nj_phylogeny_tree_mut <- genind_clone_Parents_mutOnly  %>%
  aboot(dist=provesti.dist, sample=100, tree="nj", cutoff=50, quiet=FALSE, showtree = FALSE, root=FALSE, mc.cores = 4)
nj_phylogeny_tree_mut$tip.label <- population_info_data_table$user_specimen_id[match(nj_phylogeny_tree_mut$tip.label, population_info_data_table$affy_id)]
plot.phylo(nj_phylogeny_tree_mut, label.offset=0, type = "unrooted",
           cex=0.6, font=2, lwd=4, align.tip.label=F, no.margin=T)
# Add a scale bar showing 5% difference.
add.scale.bar(0, 0, length=0.01, cex=0.65, lwd=2)
nodelabels(nj_phylogeny_tree_mut$node.label, cex=.5,  frame="n", font=3, xpd=TRUE)

#write.tree(nj_phylogeny_tree_mut, file = "nj_parent_tree_mutOnly.nwk")

###############################
### Calculate Heteozygosity ###
###############################

#subset samples from the genind object, 12 for each region

CU<-c("a550962-4368120-060520-252_I11.CEL",
      "a550962-4368120-060520-252_I07.CEL",
      "a550962-4368120-060520-252_I03.CEL",
      "a550962-4368120-060520-252_I15.CEL",
      "a550962-4368120-060520-252_I01.CEL",
      "a550962-4368120-060520-252_I09.CEL",
      "a100000-4368120-060520-256_C17.CEL",
      "a550962-4368120-060520-252_I17.CEL",
      "a100000-4368120-060520-256_I17.CEL",
      "a100000-4368120-060520-256_A17.CEL",
      "a100000-4368120-060520-256_O17.CEL",
      "a100000-4368120-060520-256_M17.CEL"
)

FL<-c("a100000-4368120-060520-256_C01.CEL",
      "a100000-4368120-060520-256_C03.CEL",
      "a550962-4368120-060520-256_K11.CEL",
      "a550962-4368120-060520-256_K17.CEL",
      "a550962-4368120-060520-256_M01.CEL",
      "a550962-4368120-060520-256_M11.CEL",
      "a550962-4368120-060520-256_M17.CEL",
      "a550962-4368120-060520-256_M21.CEL",
      "a550962-4368120-060520-256_O11.CEL",
      "a550962-4368120-060520-256_O17.CEL",
      "a550962-4393310-052921-062_A15.CEL",
      "a550962-4393310-052921-062_A23.CEL"
)

BE<-c("a550962-4368120-060520-252_G19.CEL",
      "a550962-4368120-060520-252_G21.CEL",
      "a550962-4368120-060520-252_K21.CEL",
      "a550962-4368120-060520-252_K23.CEL",
      "a550962-4368120-060520-252_M01.CEL",
      "a550962-4368120-060520-252_M15.CEL",
      "a550962-4368120-060520-252_M17.CEL",
      "a550962-4368120-060520-252_M21.CEL",
      "a550962-4368120-060520-252_O19.CEL",
      "a550962-4368120-060520-253_A03.CEL",
      "a550962-4368120-060520-253_E07.CEL",
      "a550962-4368120-060520-253_I01.CEL"
)

samples3<-c("FORMAT",ramet,rametNeighbor,automictic,self,cu_fl,BE,CU,FL)

vcf@gt<-vcf@gt[,colnames(vcf@gt) %in% samples3]

# Heterozygous alleles of all SNPs.
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE)
heterozygous_alleles_all <- apply(gt, MARGIN=2, function(x) {sum(lengths(regmatches(x, gregexpr("0/1", x))))})
heterozygous_alleles_all <- (heterozygous_alleles_all / nrow(gt)) * 100
heterozygous_alleles_all_data_frame <- data.frame(heterozygous_alleles_all)
heterozygous_alleles_all_data_table <- setDT(heterozygous_alleles_all_data_frame, keep.rownames=TRUE)[]

# Rename the rn column.
setnames(heterozygous_alleles_all_data_table, c("rn"), c("affy_id"))
# Rename the heterozygous_alleles column.
setnames(heterozygous_alleles_all_data_table, c("heterozygous_alleles_all"), c("percent_heterozygous_coral_all"))
# Round data to two digits.
heterozygous_alleles_all_data_table$percent_heterozygous_coral_all <- round(heterozygous_alleles_all_data_table$percent_heterozygous_coral_all, digits=2)

# write out table of percent heteozygosity
#write.csv(heterozygous_alleles_all_data_table, "heterozygosity.csv", row.names=F)

###########################
### Plot Heterozygosity ###
###########################

# set sequential id of samples
het<-heterozygous_alleles_all_data_table %>%
  mutate(group=ifelse(affy_id %in% CU, "CU",
                      ifelse(affy_id %in% BE, "BE",
                             ifelse(affy_id %in% FL, "FL",
                                    ifelse(affy_id %in% ramet, "ramet",
                                           ifelse(affy_id %in% rametNeighbor, "rametNeighbor",
                                                  ifelse(affy_id %in% cu_fl, "CUxFL", 
                                                         ifelse(affy_id %in% automictic, "automictic","selfed")))))))) %>%
  arrange(factor(group, levels=c("BE","FL","CUxFL","CU","ramet","rametNeighbor","selfed","automictic")),percent_heterozygous_coral_all) %>%
  mutate(idu= seq(1,length(percent_heterozygous_coral_all)))

# heterozygosity plot
ggplot(het, aes(idu,percent_heterozygous_coral_all, color=group, fill=group))+
  geom_point(size=1.7, pch=21) + ylim(5,20)+
  theme_classic() +
  xlab("")+
  ylab("Heterozygosity")+
  scale_color_manual(values=c(BE="black",CU="grey80",CUxFL="#0072B2",FL="grey40",selfed="#D55E00",
                              ramet="#FF801F",rametNeighbor="#F0E442",automictic="#E69F00"))+
  scale_fill_manual(values=c(BE="black",CU="grey80",CUxFL="#0072B250",FL="grey40",selfed="#D55E0080",
                             ramet="#FF801F80",rametNeighbor="#F0E44280",automictic="#E69F0080"))



###########################################
### Statistical Tests on Heterozygosity ###
###########################################

# ANOVA tests of subset of samples above
mod<-aov(percent_heterozygous_coral_all ~ group, data=het)
summary(mod)

# post-hoc tests
TukeyHSD(mod)

# read in heterozygosity estimates of all samples in vcf file
het_tab <- read.csv("E:/PSU/NOAA/somatic_mutation/heterozygosity.csv")

het_tab2<-het_tab %>% filter(species=="A.palmata")

# mixed-effects model of the effect of location (fixed) and plate (random)
library(lme4)
library(lmerTest)
library(emmeans)

mod2<-lmer(percent_heterozygous_coral_all~ location + (1|plate), data = het_tab2)
summary(mod2)

# pairwise comparisons with post-hoc multiple test correction
emmeans(mod2, list(pairwise ~ location), adjust = "tukey")

###########################
## Plots of DNA mixtures ##
###########################

# genetic distance estimates of DNA mixtures to donor genets
dat<-cbind.data.frame(pair=c("GAxGB","M1xGA","M1xGB","M2xGA","M2xGB"),
           gen.dist=c("0.070866342","0.001929718","0.067159252","0.013025594","0.046541743"))

# plot DNA mixture genetic distances
ggplot(dat, aes(x=pair,y=as.numeric(gen.dist))) + 
  geom_point(size=6)+ 
  geom_hline(yintercept = 0.005398558, color="grey")+ # ramet distance
  coord_flip() + xlab("")+ ylab("Pairwise Genetic Distance")+
  theme_classic()

# heterozygosity and missing data of DNA mixtures
het2<-cbind.data.frame(sample=c("GA1","GA2","M1","M2","GB1","GB2"),
                       heterozygosity=c("13.08895207","13.03310317","13.0077173","14.6831844","13.90129976","13.30239643"),
                       missing=c("0.51","0.59","0.61","2.57","0.95","0.53"),
                       group=c("GA","GA","M1","M2","GB","GB"))
het2$sample <- factor(het2$sample, levels = het2$sample)

# heterozygosity plot of DNA mixtures
ggplot(het2, aes(x=sample,y=as.numeric(heterozygosity))) + 
  geom_point(size=3)+ ylim(12.5,15)+
  xlab("")+ ylab("% Heterozygosity") +
  theme_classic()

# barplot of fixed sites between genet A and genet B of the DNA mixes
fixed_SNPs<-cbind.data.frame(sample=c("Genet A","M1","M2","Genet B"),sharedA=c("100","66.5","3.7","0"),
                             sharedB=c("0","0","0","100"),notCalled=c("0","23.1","2.8","0"),
                             heterozygous=c("0","10.4","93.5","0"))
library(reshape2)
meltFixed<-melt(fixed_SNPs,id=c("sample"))
meltFixed$value<-as.numeric(meltFixed$value)
meltFixed$sample <- factor(meltFixed$sample, levels = fixed_SNPs$sample)

ggplot(meltFixed, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="fill", stat="identity") + theme_classic()+
  scale_y_continuous(expand = c(0, 0)) + ylab("Proportion of Fixed SNPs")+
  xlab("") +  scale_fill_manual(values=c("#F79028", "#8D4C9E","grey","black"))


######################################
### Calculate Recombination Events ###
######################################
library(hsphase)

#extract genotypes and convert format of alleles
geno <- extract.gt(vcf)

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

sum(is.na(G)) # should not contain NAs
o3 <- G
o3[is.na(o3)] <- 9
sum(is.na(o3))
dim(o3)

SNPdata <- t(o3)
colnames(SNPdata)<-row.names(geno)
row.names(SNPdata)<-colnames(geno)
SNPdata<-SNPdata[row.names(SNPdata) %in% c(self,automictic),]

# calculate blocks
recombinationBlocks <- bmh(SNPdata)
recombinationBlocks[1:5,1:3]

# calculate number of recombination events
rec<-cbind(row.names(SNPdata),recombinations(bmh(SNPdata))) 

colnames(rec)<-c("affy_id","recombinations")

#write.table(rec, "recombination_events_offspring.txt",sep="\t",row.names=F)

#################################
### Plot Recombination Events ###
#################################

myCol<-c("#0072B2","#D55E00","#FFA35C")

mat<-as.matrix(recombinationBlocks)

uniquelength <- sapply(as.data.frame(mat), function(x) length(unique(x)))
mat <- subset(mat, select=uniquelength > 1)

mat<-mat[, colSums(is.na(mat))< 3 ]

t2<- as.data.frame(mat) %>%
  mutate(Sample=row.names(mat)) %>%
  left_join(het %>% dplyr::select(affy_id, group), by=c("Sample"="affy_id"))

t3<-t2%>%
  dplyr::left_join(as.data.frame(rec) %>% dplyr::select(affy_id, recombinations), by=c("Sample"="affy_id"))

order<-c("automictic","selfed")


## Plot heat map 
library(ComplexHeatmap)

ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 6),
       heatmap_row_title_gp = gpar(fontsize = 6),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)

ha = rowAnnotation(RR = anno_barplot(as.numeric(t3$recombinations), height = unit(1.5, "cm"), gp = gpar(fill = "black")))

hm<-Heatmap(mat, name = "Parental haplotype", col=myCol, row_split = factor(t2$group, levels=order),cluster_row_slices = FALSE, 
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, show_column_names =FALSE,row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = FALSE, heatmap_legend_param = list(direction = "horizontal"),
            row_title_rot = 0,row_gap = unit(1, "mm"),width = unit(15, "cm"), height = unit(8, "cm"),
            right_annotation = ha)
draw(hm)