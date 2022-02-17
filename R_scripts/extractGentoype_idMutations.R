##################
### Import VCF ###
##################
library(vcfR)
library(dplyr)

#vcf <- read.vcfR("E:/PSU/NOAA/PRO100175_PSU175_SAX_b04/PRO100175_PSU175_SAX_b04/dup_removed_DB_snps_P1-5.vcf")

setwd("C:/Users/Sheila's Comp/Desktop")
vcf <- read.vcfR("STAGdb_04Jan2022.vcf")

##################
### Subset VCF ###
##################

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


offspring<-c("a550962-4393310-052921-062_O03.CEL","a550962-4393310-052921-062_E05.CEL",
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
         "a550962-4393310-052921-062_C05.CEL"
)

sire<-c("a100000-4368120-060520-256_C03.CEL",
        "120207031_(Axiom_AcropSNP)_C15.CEL",
        "120207034_(Axiom_AcropSNP)_C21.CEL")

samples<-c("FORMAT",ramet,rametNeighbor,offspring, cu_fl,sire)

vcf@gt<-vcf@gt[,colnames(vcf@gt) %in% samples]

###################
### Filter SNPs ###
###################

probAnnotation<-read.csv("E:/PSU/NOAA/PRO100175_PSU175_SAX_b01/PRO100175_PSU175_SAX_b01/AxiomReference/Axiom_AcropSNP_coral_Annotation.r1.csv/Axiom_AcropSNP_coral_Annotation.r1.csv",
                         skip=18)

gt <- extract.gt(vcf, element="GT", as.numeric=FALSE)
gt[1:3,1:3]

# extract parent genotypes (ramet and rametNeighbor samples)
parent_gt<-gt[,colnames(gt) %in% c(ramet,rametNeighbor)]

noDups_parent<-as.data.frame(parent_gt) %>%
  mutate(same = ifelse(rowSums(.[,c(1:15)] == "1/1") == ncol(.), "allSame", 
                       ifelse(rowSums(.[,c(1:15)] == "0/0") == ncol(.), "allSame",
                              ifelse(rowSums(.[,c(1:15)] == "0/1") == ncol(.), "allSame","diff")))) %>%
  mutate(probe=row.names(parent_gt)) %>%
  relocate(probe,.before="a550962-4393310-052921-062_A01.CEL") %>%
  left_join(probAnnotation %>% dplyr::select(Probe.Set.ID, Chromosome, Physical.Position), by=c("probe"="Probe.Set.ID"))

probes<-noDups_parent %>%
  filter(same == "diff")

# probes with mutations in parental samples
#write.table(as.data.frame(probes[,1]),"smutations.txt", sep="\t", quote=F,row.names=F)

#####
# extract sire genotypes that differ from dam genet
parentSire_gt<-gt[,colnames(gt) %in% c(sire)]

# compare to parental genotypes
noDups_parentSire<-as.data.frame(parentSire_gt) %>%
  mutate(probe=row.names(.)) %>%
  left_join(noDups_parent,by=c("probe")) %>%
  mutate(sire2dam_alleles = ifelse(rowSums(.[,c(1:3,5:19)] == "1/1") == ncol(.[,c(1:3,5:19)]), "allSame", 
                       ifelse(rowSums(.[,c(1:3,5:19)] == "0/0") == ncol(.[,c(1:3,5:19)]), "allSame",
                              ifelse(rowSums(.[,c(1:3,5:19)] == "0/1") == ncol(.[,c(1:3,5:19)]), "allSame","diff"))),
         sire2dam_alleles=ifelse(is.na(sire2dam_alleles),"missing",sire2dam_alleles)) %>%
  filter(!sire2dam_alleles == "missing") %>%
  filter(!sire2dam_alleles == "allSame") %>%
  mutate(num_00_sire=rowSums(.[,c(1:3)] == "0/0"), num_01_sire=rowSums(.[,c(1:3)] == "0/1"),
         num_11_sire=rowSums(.[,c(1:3)] == "1/1")) %>%
  mutate(num_00_dam=rowSums(.[,c(5:19)] == "0/0"), num_01_dam=rowSums(.[,c(5:19)] == "0/1"),num_11_dam=rowSums(.[,c(5:19)] == "1/1"))

####################################
### Calculate Number of Changes  ###
###    for each Mutation         ###
####################################

# in parent samples and neighboring colonies
cnt_parents<-probes[,colnames(probes) %in% c(ramet,rametNeighbor)] %>%
  mutate(num_00=rowSums(. == "0/0"), num_01=rowSums(. == "0/1"),num_11=rowSums(. == "1/1"), probeID=probes$probe)

# in offspring
allSample_mut<-gt[row.names(gt) %in% probes$probe,]

cnt_offspring <- as.data.frame(allSample_mut[,colnames(allSample_mut) %in% offspring]) %>%
  mutate(num_00=rowSums(. == "0/0",na.rm = TRUE), num_01=rowSums(. == "0/1",na.rm = TRUE),
         num_11=rowSums(. == "1/1",na.rm = TRUE), num_NA=30-(num_00 + num_01 +num_11), probeID=rownames(allSample_mut))

#####################################
### Calculate Number of SM Changes ###
###      in CUxFL offspring       ###
#####################################

# CUxFL offspring with somatic mutations
cuFL_mut<-gt[row.names(gt) %in% probes$probe,colnames(gt) %in% cu_fl]

cnt_cu_fl <- as.data.frame(cuFL_mut) %>%
  mutate(num_00=as.numeric(rowSums(. == "0/0",na.rm = TRUE)), num_01=as.numeric(rowSums(. == "0/1",na.rm = TRUE)),
         num_11=as.numeric(rowSums(. == "1/1",na.rm = TRUE)), num_NA=as.numeric(11-(num_00 + num_01 +num_11)), probe=rownames(cuFL_mut)) %>%
  left_join(noDups_parentSire %>% select(num_00_sire, num_01_sire, num_11_sire,
                                         num_00_dam, num_01_dam, num_11_dam,probe),by=c("probe")) %>%
  mutate(dam_major_allele=ifelse(num_11_dam > num_01_dam & num_11_dam > num_00_dam, "11",
                                 ifelse(num_01_dam > num_11_dam & num_01_dam > num_00_dam, "01", "00")),
         sire_major_allele=ifelse(num_11_sire > num_01_sire & num_11_sire > num_00_sire, "11",
                                  ifelse(num_01_sire > num_11_sire & num_01_sire > num_00_sire, "01", "00")),
         offspring_prob_00=ifelse(sire_major_allele == "00" & dam_major_allele == "00", (1*11), 
                                  ifelse(sire_major_allele == "11" & dam_major_allele == "11",0,
                                         ifelse(sire_major_allele == "01" & dam_major_allele == "01", (.25*11),
                                                ifelse(sire_major_allele == "00" & dam_major_allele == "01", (.5*11),
                                                       ifelse(sire_major_allele == "01" & dam_major_allele == "00", (.5*11),0))))),
         offspring_prob_01=ifelse(sire_major_allele == "01" & dam_major_allele == "01", (.5*11),
                                  ifelse(sire_major_allele == "11" & dam_major_allele == "01", (.5*11),
                                         ifelse(sire_major_allele == "01" & dam_major_allele == "11", (.5*11),
                                                ifelse(sire_major_allele == "00" & dam_major_allele == "01", (.5*11),
                                                       ifelse(sire_major_allele == "01" & dam_major_allele == "00", (.5*11),
                                                              ifelse(sire_major_allele == "11" & dam_major_allele == "00", (1*11),
                                                                     ifelse(sire_major_allele == "00" & dam_major_allele == "11", (1*11),0))))))),
         offspring_prob_11=ifelse(sire_major_allele == "00" & dam_major_allele == "00", 0, 
                                  ifelse(sire_major_allele == "11" & dam_major_allele == "11",(1*11),
                                         ifelse(sire_major_allele == "01" & dam_major_allele == "01", (.25*11),
                                                ifelse(sire_major_allele == "11" & dam_major_allele == "01", (.5*11),
                                                       ifelse(sire_major_allele == "01" & dam_major_allele == "11", (.5*11),0)))))) %>%
  select(probe,num_00,num_01,num_11,num_NA,offspring_prob_00,offspring_prob_01,offspring_prob_11, 
         dam_major_allele, num_00_dam,num_01_dam,num_11_dam, sire_major_allele,num_00_sire,num_01_sire,num_11_sire)

#write.table(cnt_cu_fl,"cufl_somatic_mutations_alleles.txt", sep="\t",row.names = FALSE)

#####################################
### Calculate Number of Mutations ###
###     unique to offspring       ###
#####################################

# extract offspring genotypes
offspring_gt<-gt[,colnames(gt) %in% c(offspring)]

# add column to identify those that differ between offspring
noDups_offspring<-as.data.frame(offspring_gt) %>%
  mutate(same = ifelse(rowSums(.[,c(1:30)] == "1/1") == ncol(.), "allSame", 
                       ifelse(rowSums(.[,c(1:30)] == "0/0") == ncol(.), "allSame",
                              ifelse(rowSums(.[,c(1:30)] == "0/1") == ncol(.), "allSame","diff")))) %>%
  mutate(probe=row.names(offspring_gt)) %>%
  relocate(probe,.before="a550962-4393310-052921-062_A07.CEL") %>%
  left_join(probAnnotation %>% dplyr::select(Probe.Set.ID, Chromosome, Physical.Position), by=c("probe"="Probe.Set.ID")) %>%
  left_join(noDups_parent %>% select("probe","same", "a550962-4393310-052921-062_O01.CEL"), by=c("probe")) %>%
  rename("same_offspring"=same.x)%>%
  rename("same_parent"=same.y)

# subset table to only those that differ in offspring but not in parent samples
probes_offspring<-noDups_offspring %>%
  filter(same_offspring == "diff" & same_parent=="allSame")

# probes with mutations in offspring samples
write.table(as.data.frame(probes_offspring),"uniq_offspring.txt", sep="\t", quote=F,row.names=F)

# tabulate the number of mutations for each probe
cnt_offspring_uniq<-probes_offspring[,colnames(probes_offspring) %in% offspring] %>%
  mutate(num_00=rowSums(. == "0/0"), num_01=rowSums(. == "0/1"),num_11=rowSums(. == "1/1"), probeID=probes_offspring$probe)

write.table(as.data.frame(cnt_offspring_uniq),"uniq_offspring_cnt.txt", sep="\t", quote=F,row.names=F)
