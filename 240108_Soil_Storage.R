----------------Begin--------------------------
library(ggplot2)
library(plyr)
library(gridExtra)
library(car)
library(visreg)
library(lme4)
library(MASS) 
library(reshape) 
library(ade4)
library(hillR)
library(lmerTest)
library(emmeans)
library(vegan)
library(reshape2)
library(pbkrtest)
library(svglite)
library(mvabund)
library(egg)




#------read data short-term -------
setwd("~/Storage_repository")

ITS_data_raw <- read.table("Shortterm_ITS_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")

ITS_data <- t(subset(ITS_data_raw, rowSums(ITS_data_raw) > 0))

ITS_meta <- read.csv("Shortterm_ITS_env.csv", header = T, row.names = 1)


AMF_data_raw <- read.table("Shortterm_AMF_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")

AMF_data <- t(subset(AMF_data_raw, rowSums(AMF_data_raw) > 0))

AMF_meta <- read.csv("Shortterm_AMF_env.csv", header = T, row.names = 1)


BAC_data_raw <- read.table("Shortterm_16S_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")

BAC_data <- t(subset(BAC_data_raw, rowSums(BAC_data_raw) > 0))

BAC_meta <- read.csv("Shortterm_16S_env.csv", header = T, row.names = 1)





#----------- alpha diversity short-term-----------

#q=0 (richness)
#Create vector of average alpha diversity (q=0) for each sample (row)
ITS_meta$alphaq0<-apply(ITS_data,1,hill_taxa,q=0) #Average richness (q=0) in each sample

AMF_meta$alphaq0<-apply(AMF_data,1,hill_taxa,q=0) #Average richness (q=0) in each sample

BAC_meta$alphaq0<-apply(BAC_data,1,hill_taxa,q=0) #Average richness (q=0) in each sample


#q=1 (Shannon entropy; weighted according to rel. abund)
ITS_meta$alphaq1<-apply(ITS_data,1,hill_taxa,q=1) #Average Shannon entropy (q=1) in each sample

AMF_meta$alphaq1<-apply(AMF_data,1,hill_taxa,q=1) #Average Shannon entropy (q=1) in each sample

BAC_meta$alphaq1<-apply(BAC_data,1,hill_taxa,q=1) #Average Shannon entropy (q=1) in each sample


#q=2 (Inverse simpsons; downweighting rare taxa)
ITS_meta$alphaq2<-apply(ITS_data,1,hill_taxa,q=2) #Average inverse Simpson's (q=2) in each sample

AMF_meta$alphaq2<-apply(AMF_data,1,hill_taxa,q=2) #Average inverse Simpson's (q=2) in each sample

BAC_meta$alphaq2<-apply(BAC_data,1,hill_taxa,q=2) #Average inverse Simpson's (q=2) in each sample





ITS_meta$name <- as.numeric(row.names(ITS_meta))
AMF_meta$name <- as.numeric(row.names(AMF_meta))
BAC_meta$name <- as.numeric(row.names(BAC_meta))


ITS_ordered <- ITS_meta[order(ITS_meta$name),]
AMF_ordered <- AMF_meta[order(AMF_meta$name),]
BAC_ordered <- BAC_meta[order(BAC_meta$name),]


ITS_mol <- subset(ITS_ordered, sampleStorage == "molecular")

ITS_bulk <- subset(ITS_ordered, sampleStorage == "bulk")

AMF_mol <- subset(AMF_ordered, sampleStorage == "molecular")

AMF_bulk <- subset(AMF_ordered, sampleStorage == "bulk")

BAC_mol <- subset(BAC_ordered, sampleStorage == "molecular")

BAC_bulk <- subset(BAC_ordered, sampleStorage == "bulk")


#Perform t-test for differences in raw reads per storage treatment
t.test(ITS_mol$read.tot,ITS_bulk$read.tot,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_mol$read.tot,AMF_bulk$read.tot,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(BAC_mol$read.tot,BAC_bulk$read.tot,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

#Perform t-test for differences in % read retained per storage treatment
t.test(ITS_mol$input,ITS_bulk$input,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_mol$input,AMF_bulk$input,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(BAC_mol$input,BAC_bulk$input,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364



#Perform t-test for differences in % read retained per storage treatment
t.test(ITS_mol$final_perc_reads_retained,ITS_bulk$final_perc_reads_retained,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_mol$final_perc_reads_retained,AMF_bulk$final_perc_reads_retained,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(BAC_mol$final_perc_reads_retained,BAC_bulk$final_perc_reads_retained,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364



#---------- beta diversity short-term ----------

#Create distance matrix for beta at q=0, q=1, q=2
#Output is a number between 1 and 2 (effective number of communities; -1 for proportional turnover)
ITS_ddatq0 <-  (as.dist(hill_taxa_parti_pairwise(ITS_data,q=0,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq0 <-  (as.dist(hill_taxa_parti_pairwise(AMF_data,q=0,output="matrix",pairs="full")$TD_beta)) -1

BAC_ddatq0 <-  (as.dist(hill_taxa_parti_pairwise(BAC_data,q=0,output="matrix",pairs="full")$TD_beta)) -1


ITS_ddatq1 <-  (as.dist(hill_taxa_parti_pairwise(ITS_data,q=1,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq1 <-  (as.dist(hill_taxa_parti_pairwise(AMF_data,q=1,output="matrix",pairs="full")$TD_beta)) -1

BAC_ddatq1 <-  (as.dist(hill_taxa_parti_pairwise(BAC_data,q=1,output="matrix",pairs="full")$TD_beta)) -1


ITS_ddatq2 <-  (as.dist(hill_taxa_parti_pairwise(ITS_data,q=2,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq2 <-  (as.dist(hill_taxa_parti_pairwise(AMF_data,q=2,output="matrix",pairs="full")$TD_beta)) -1

BAC_ddatq2 <-  (as.dist(hill_taxa_parti_pairwise(BAC_data,q=2,output="matrix",pairs="full")$TD_beta)) -1



#dbRDA for storage effects 

set.seed(420)

ITS_meta$pair <- as.factor(ITS_meta$pair)
AMF_meta$pair <- as.factor(ITS_meta$pair)
BAC_meta$pair <- as.factor(ITS_meta$pair)


ITS_modq0 <- dbrda(ITS_ddatq0~sampleStorage  + Condition(pair), data=ITS_meta)
anova(ITS_modq0,strata=ITS_meta$pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


AMF_modq0 <- dbrda(AMF_ddatq0~sampleStorage  + Condition(pair), data=AMF_meta)
anova(AMF_modq0,strata=AMF_meta$pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

BAC_modq0 <- dbrda(BAC_ddatq0~sampleStorage  + Condition(pair), data=BAC_meta)
anova(BAC_modq0,strata=BAC_meta$pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


ITS_modq1 <- dbrda(ITS_ddatq1~sampleStorage  + Condition(pair), data=ITS_meta)
anova(ITS_modq1,strata=ITS_meta$pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

AMF_modq1 <- dbrda(AMF_ddatq1~sampleStorage  + Condition(pair), data=AMF_meta)
anova(AMF_modq1,strata=AMF_meta$pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

BAC_modq1 <- dbrda(BAC_ddatq1~sampleStorage  + Condition(pair), data=BAC_meta)
anova(BAC_modq1,strata=BAC_meta$pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



ITS_modq2 <- dbrda(ITS_ddatq2~sampleStorage  + Condition(pair), data=ITS_meta)
anova(ITS_modq2,strata=ITS_meta$pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

AMF_modq2 <- dbrda(AMF_ddatq2~sampleStorage  + Condition(pair), data=AMF_meta)
anova(AMF_modq2,strata=AMF_meta$pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

BAC_modq2 <- dbrda(BAC_ddatq2~sampleStorage  + Condition(pair), data=BAC_meta)
anova(BAC_modq2,strata=BAC_meta$pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    





#------ read data long-term ---------


ITS_data_raw <- read.table("Longterm_ITS_counts.csv", header = T, row.names = 1, check.names = F, sep = ",")

ITS_data <- t(subset(ITS_data_raw, rowSums(ITS_data_raw) > 0))

ITS_meta <- read.csv("Longterm_ITS_env.csv", header = T, row.names = 1)


AMF_data_raw <- read.table("Longterm_AMF_counts.csv", header = T, row.names = 1, check.names = F, sep = ",")

AMF_data <- t(subset(AMF_data_raw, rowSums(AMF_data_raw) > 0))


AMF_meta <- read.csv("Longterm_AMF_env.csv", header = T, row.names = 1)


#---- alpha diversity long-term--------

ITS_meta$alphaq0<-apply(ITS_data,1,hill_taxa,q=0) #Average richness (q=0) in each sample

AMF_meta$alphaq0<-apply(AMF_data,1,hill_taxa,q=0) #Average richness (q=0) in each sample

#q=1 (Shannon entropy; weighted according to rel. abund)
ITS_meta$alphaq1<-apply(ITS_data,1,hill_taxa,q=1) #Average Shannon entropy (q=1) in each sample

AMF_meta$alphaq1<-apply(AMF_data,1,hill_taxa,q=1) #Average Shannon entropy (q=1) in each sample


#q=2 (Inverse simpsons; downweighting rare taxa)
ITS_meta$alphaq2<-apply(ITS_data,1,hill_taxa,q=2) #Average inverse Simpson's (q=2) in each sample

AMF_meta$alphaq2<-apply(AMF_data,1,hill_taxa,q=2) #Average inverse Simpson's (q=2) in each sample


ITS_ordered <- ITS_meta[order(ITS_meta$Pair),]
AMF_ordered <- AMF_meta[order(AMF_meta$Pair),]


ITS_dry <- subset(ITS_ordered, DryFrozen == "Dry")

ITS_frozen <- subset(ITS_ordered, DryFrozen == "Frozen")

AMF_dry <- subset(AMF_ordered, DryFrozen == "Dry")

AMF_frozen <- subset(AMF_ordered, DryFrozen == "Frozen")


#Perform t-test for differences in raw reads per storage treatment
t.test(ITS_dry$input,ITS_frozen$input,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_dry$input,AMF_frozen$input,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

#Perform t-test for differences in raw reads per storage treatment
t.test(ITS_dry$final_perc_reads_retained,ITS_frozen$final_perc_reads_retained,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_dry$final_perc_reads_retained,AMF_frozen$final_perc_reads_retained,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

#Perform t-test for differences in raw reads per storage treatment
t.test(ITS_dry$nonchim,ITS_frozen$nonchim,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_dry$nonchim,AMF_frozen$nonchim,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364




#Perform t-test for differences in alpha diversity per storage treatment

#q=0
t.test(ITS_dry$alphaq0,ITS_frozen$alphaq0,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_dry$alphaq0,AMF_frozen$alphaq0,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364


#q=1
t.test(ITS_dry$alphaq1,ITS_frozen$alphaq1,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364


t.test(AMF_dry$alphaq1,AMF_frozen$alphaq1,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364



#q=2
t.test(ITS_dry$alphaq2,ITS_frozen$alphaq2,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364

t.test(AMF_dry$alphaq2,AMF_frozen$alphaq2,paired=TRUE) #t = 1.7294, df = 15, p-value = 0.1043 #t = -0.47622, df = 41, p-value = 0.6364



#------  beta diversity long-term -----

#Create distance matrix for beta at q=0, q=1, q=2
#Output is a number between 1 and 2 (effective number of communities; -1 for proportional turnover)
ITS_ddatq0 <-  (as.dist(hill_taxa_parti_pairwise(ITS_data,q=0,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq0 <-  (as.dist(hill_taxa_parti_pairwise(AMF_data,q=0,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq1 <-  (as.dist(hill_taxa_parti_pairwise(ITS_data,q=1,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq1 <-  (as.dist(hill_taxa_parti_pairwise(AMF_data,q=1,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq2 <-  (as.dist(hill_taxa_parti_pairwise(ITS_data,q=2,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq2 <-  (as.dist(hill_taxa_parti_pairwise(AMF_data,q=2,output="matrix",pairs="full")$TD_beta)) -1



#dbRDA for storage effects 

set.seed(420)

AMF_meta$Pair <- as.factor(AMF_meta$Pair)

ITS_modq0 <- dbrda(ITS_ddatq0~DryFrozen  + Condition(Pairs), data=ITS_meta)
anova(ITS_modq0,strata=ITS_meta$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

AMF_modq0 <- dbrda(AMF_ddatq0~DryFrozen  + Condition(Pair), data=AMF_meta)
anova(AMF_modq0,strata=AMF_meta$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

RsquareAdj(ITS_modq0)
RsquareAdj(AMF_modq0)


ITS_modq1 <- dbrda(ITS_ddatq1~DryFrozen  + Condition(Pairs), data=ITS_meta)
anova(ITS_modq1,strata=ITS_meta$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

AMF_modq1 <- dbrda(AMF_ddatq1~DryFrozen  + Condition(Pair), data=AMF_meta)
anova(AMF_modq1,strata=AMF_meta$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

RsquareAdj(ITS_modq1)
RsquareAdj(AMF_modq1)


ITS_modq2 <- dbrda(ITS_ddatq2~DryFrozen  + Condition(Pairs), data=ITS_meta)
anova(ITS_modq2,strata=ITS_meta$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

AMF_modq2 <- dbrda(AMF_ddatq2~DryFrozen  + Condition(Pair), data=AMF_meta)
anova(AMF_modq2,strata=AMF_meta$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

RsquareAdj(ITS_modq2)
RsquareAdj(AMF_modq2)


#-------------- 2013 vs 2020 --------


ITS_all <- data.frame(ITS_meta, ITS_data)


ITS_all <- ITS_all[order(ITS_all$Pairs),]

old <- subset(ITS_all, SampleYear == 2013)


ITS_meta_old <- old[,0:20]
ITS_data_raw_old <- t(old[,21:ncol(old)])

ITS_data_old <- t(subset(ITS_data_raw_old, rowSums(ITS_data_raw_old) > 0))

AMF_all <- data.frame(AMF_meta, AMF_data)
AMF_all <- AMF_all[order(AMF_all$Pair),]
old <- subset(AMF_all, SampleYear == 2013)

AMF_meta_old <- old[,0:17]
AMF_data_raw_old <- t(old[,18:ncol(old)])
AMF_data_old <- t(subset(AMF_data_raw_old, rowSums(AMF_data_raw_old) > 0))


#Create distance matrix for beta at q=0, q=1, q=2
#Output is a number between 1 and 2 (effective number of communities; -1 for proportional turnover)
ITS_ddatq0_old <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_old,q=0,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq0_old <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_old,q=0,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq1_old <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_old,q=1,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq1_old <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_old,q=1,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq2_old <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_old,q=2,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq2_old <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_old,q=2,output="matrix",pairs="full")$TD_beta)) -1



ITS_modq0_old <- dbrda(ITS_ddatq0_old~DryFrozen  + Condition(Pairs), data=ITS_meta_old)
anova(ITS_modq0_old,strata=ITS_meta_old$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


AMF_modq0_old <- dbrda(AMF_ddatq0_old~DryFrozen  + Condition(Pair), data=AMF_meta_old)
anova(AMF_modq0_old,strata=AMF_meta_old$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


ITS_modq1_old <- dbrda(ITS_ddatq1_old~DryFrozen  + Condition(Pairs), data=ITS_meta_old)
anova(ITS_modq1_old,strata=ITS_meta_old$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


AMF_modq1_old <- dbrda(AMF_ddatq1_old~DryFrozen  + Condition(Pair), data=AMF_meta_old)
anova(AMF_modq1_old,strata=AMF_meta_old$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


ITS_modq2_old <- dbrda(ITS_ddatq2_old~DryFrozen  + Condition(Pairs), data=ITS_meta_old)
anova(ITS_modq2_old,strata=ITS_meta_old$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



AMF_modq2_old <- dbrda(AMF_ddatq2_old~DryFrozen  + Condition(Pair), data=AMF_meta_old)
anova(AMF_modq2_old,strata=AMF_meta_old$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


#new

new <- subset(ITS_all, SampleYear == 2020)


ITS_meta_new <- new[,0:20]
ITS_data_raw_new <- t(new[,21:ncol(new)])

ITS_data_new <- t(subset(ITS_data_raw_new, rowSums(ITS_data_raw_new) > 0))

AMF_all <- data.frame(AMF_meta, AMF_data)
AMF_all <- AMF_all[order(AMF_all$Pair),]
new <- subset(AMF_all, SampleYear == 2020)

AMF_meta_new <- new[,0:17]
AMF_data_raw_new <- t(new[,18:ncol(new)])
AMF_data_new <- t(subset(AMF_data_raw_new, rowSums(AMF_data_raw_new) > 0))


#Create distance matrix for beta at q=0, q=1, q=2
#Output is a number between 1 and 2 (effective number of communities; -1 for proportional turnover)
ITS_ddatq0_new <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_new,q=0,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq0_new <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_new,q=0,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq1_new <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_new,q=1,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq1_new <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_new,q=1,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq2_new <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_new,q=2,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq2_new <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_new,q=2,output="matrix",pairs="full")$TD_beta)) -1



ITS_modq0_new <- dbrda(ITS_ddatq0_new~DryFrozen  + Condition(Pairs), data=ITS_meta_new)
anova(ITS_modq0_new,strata=ITS_meta_new$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    

#anova(ITS_modq0,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



AMF_modq0_new <- dbrda(AMF_ddatq0_new~DryFrozen  + Condition(Pair), data=AMF_meta_new)
anova(AMF_modq0_new,strata=AMF_meta_new$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    
#anova(AMF_modq0,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    




ITS_modq1_new <- dbrda(ITS_ddatq1_new~DryFrozen  + Condition(Pairs), data=ITS_meta_new)
anova(ITS_modq1_new,strata=ITS_meta_new$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    
#anova(ITS_modq1,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



AMF_modq1_new <- dbrda(AMF_ddatq1_new~DryFrozen  + Condition(Pair), data=AMF_meta_new)
anova(AMF_modq1_new,strata=AMF_meta_new$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    
#anova(AMF_modq1,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


ITS_modq2_new <- dbrda(ITS_ddatq2_new~DryFrozen  + Condition(Pairs), data=ITS_meta_new)
anova(ITS_modq2_new,strata=ITS_meta_new$Pairs,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    
#anova(ITS_modq2,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



AMF_modq2_new <- dbrda(AMF_ddatq2_new~DryFrozen  + Condition(Pair), data=AMF_meta_new)
anova(AMF_modq2_old,strata=AMF_meta_old$Pair,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    
#anova(AMF_modq2,permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



#frozen vs dry - environmental variables 


#--------  frozen v dry ------------


froze <- subset(ITS_all, DryFrozen == "Frozen")
dry <- subset(ITS_all, DryFrozen == "Dry")


ITS_meta_froze <- froze[,0:20]
ITS_data_raw_froze <- t(froze[,21:ncol(froze)])


ITS_meta_dry <- dry[,0:20]
ITS_data_raw_dry <- t(dry[,21:ncol(dry)])



ITS_data_froze <- t(subset(ITS_data_raw_froze, rowSums(ITS_data_raw_froze) > 0))
ITS_data_dry <- t(subset(ITS_data_raw_dry, rowSums(ITS_data_raw_dry) > 0))


ITS_ddatq0_f <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_froze,q=0,output="matrix",pairs="full")$TD_beta)) -1
ITS_ddatq0_d <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_dry,q=0,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq1_f <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_froze,q=1,output="matrix",pairs="full")$TD_beta)) -1
ITS_ddatq1_d <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_dry,q=1,output="matrix",pairs="full")$TD_beta)) -1

ITS_ddatq2_f <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_froze,q=2,output="matrix",pairs="full")$TD_beta)) -1
ITS_ddatq2_d <-  (as.dist(hill_taxa_parti_pairwise(ITS_data_dry,q=2,output="matrix",pairs="full")$TD_beta)) -1



mantel(ITS_ddatq0_f, ITS_ddatq0_d, method = "pearson", permutations = 999)
mantel(ITS_ddatq1_f, ITS_ddatq1_d, method = "pearson", permutations = 999)
mantel(ITS_ddatq2_f, ITS_ddatq2_d, method = "pearson", permutations = 999)

colnames(AMF_all)[4] = "pctECM"

froze <- subset(AMF_all, DryFrozen == "Frozen")
dry <- subset(AMF_all, DryFrozen == "Dry")


AMF_meta_froze <- froze[,0:17]
AMF_data_raw_froze <- t(froze[,18:ncol(froze)])


AMF_meta_dry <- dry[,0:17]
AMF_data_raw_dry <- t(dry[,18:ncol(dry)])



AMF_data_froze <- t(subset(AMF_data_raw_froze, rowSums(AMF_data_raw_froze) > 0))
AMF_data_dry <- t(subset(AMF_data_raw_dry, rowSums(AMF_data_raw_dry) > 0))


AMF_ddatq0_f <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_froze,q=0,output="matrix",pairs="full")$TD_beta)) -1
AMF_ddatq0_d <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_dry,q=0,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq1_f <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_froze,q=1,output="matrix",pairs="full")$TD_beta)) -1
AMF_ddatq1_d <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_dry,q=1,output="matrix",pairs="full")$TD_beta)) -1

AMF_ddatq2_f <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_froze,q=2,output="matrix",pairs="full")$TD_beta)) -1
AMF_ddatq2_d <-  (as.dist(hill_taxa_parti_pairwise(AMF_data_dry,q=2,output="matrix",pairs="full")$TD_beta)) -1



mantel(AMF_ddatq0_f, AMF_ddatq0_d, method = "pearson", permutations = 999)
mantel(AMF_ddatq1_f, AMF_ddatq1_d, method = "pearson", permutations = 999)
mantel(AMF_ddatq2_f, AMF_ddatq2_d, method = "pearson", permutations = 999)


ITS_modq0_f <- dbrda(ITS_ddatq0_f~SampleYear+ Depth + Site + pctECM , data=ITS_meta_froze, add = TRUE)
anova(ITS_modq0_f, type = 3, by = "margin",permutations = how(nperm=9999)) 


AMF_modq0_f <- dbrda(AMF_ddatq0_f~SampleYear + Depth + Site + pctECM , data=AMF_meta_froze)
anova(AMF_modq0_f,type =3, by = "margin",permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    


ITS_modq0_d <- dbrda(ITS_ddatq0_d~SampleYear+ Depth + Site + pctECM , data=ITS_meta_dry, add = TRUE)
anova(ITS_modq0_d, type = 3, by = "margin",permutations = how(nperm=9999)) 


AMF_modq0_d <- dbrda(AMF_ddatq0_d~SampleYear + Depth + Site + pctECM , data=AMF_meta_dry)
anova(AMF_modq0_d,type =3, by = "margin",permutations = how(nperm=9999)) #SumofSqs = 0.21021, F = 1.0567, Pr(>F) =  0.0399*,          Df SumOfSqs      F Pr(>F)    



ITS_modq1_f <- dbrda(ITS_ddatq1_f~SampleYear+ Depth + Site + pctECM , data=ITS_meta_froze)
anova(ITS_modq1_f, by = "margin",permutations = how(nperm=9999)) 


AMF_modq1_f <- dbrda(AMF_ddatq1_f~Depth + Site + pctECM +SampleYear, data=AMF_meta_froze)
anova(AMF_modq1_f, by = "margin",permutations = how(nperm=9999))


ITS_modq1_d <- dbrda(ITS_ddatq1_d~SampleYear+ Depth + Site + pctECM , data=ITS_meta_dry)
anova(ITS_modq1_f, by = "margin",permutations = how(nperm=9999)) 


AMF_modq1_d <- dbrda(AMF_ddatq1_d~Depth + Site + pctECM +SampleYear, data=AMF_meta_dry)
anova(AMF_modq1_d, by = "margin",permutations = how(nperm=9999))


ITS_modq2_f <- dbrda(ITS_ddatq2_f~SampleYear + Depth + Site + pctECM, data=ITS_meta_froze)
anova(ITS_modq2_f,by = "margin",permutations = how(nperm=9999)) 


AMF_modq2_f <- dbrda(AMF_ddatq2_f~SampleYear + Depth + Site + pctECM , data=AMF_meta_froze)
anova(AMF_modq2_f, by = "margin",permutations = how(nperm=9999)) 

ITS_modq2_d <- dbrda(ITS_ddatq2_d~SampleYear + Depth + Site + pctECM, data=ITS_meta_dry)
anova(ITS_modq2_d,by = "margin",permutations = how(nperm=9999)) 

AMF_modq2_d <- dbrda(AMF_ddatq2_d~SampleYear + Depth + Site + pctECM , data=AMF_meta_dry)
anova(AMF_modq2_d, by = "margin",permutations = how(nperm=9999)) 
