# download data and make train and test sets 
# set working directory with files
#######################
setwd("~/Desktop/texas workshop")
load("phenos.Rdata")
load("snps.Rdata")
phenos<-subset(phenos,YEAR == "2013"& LOCATION == "B" & REP == 1 )
phenos1<-phenos[!duplicated(phenos$CLONE),]


# TAKE A SUBSET OF GENOTYPES AS "TRAINING DATA" --------------------------------------------------------------------
genos<-as.character(unique(phenos1$CLONE)) 
set.seed(123)
train<-sample(genos,length(genos)*.8) # Set 4/5 as training
test<-setdiff(genos,train) # Set the rest as test

# training data
phenos.trn=phenos1[which(phenos1$CLONE %in% train),]
phenos.trn$CLONE<-droplevels(phenos.trn$CLONE) 
row.names(phenos.trn)<-phenos.trn$CLONE
train.data<-merge(phenos.trn,snps,by=0)
snps.train<-train.data[11:1010]

#validation-test data
phenos.tst=phenos1[which(phenos1$CLONE %in% test),]
phenos.tst$CLONE<-droplevels(phenos.tst$CLONE) 
row.names(phenos.tst)<-phenos.tst$CLONE

test.data<-merge(phenos.tst,snps,by=0)
snps.test<-test.data[11:1010]

# start testing GWAS: fit CBSDRS training data with first SNP
# You can get the SNP with the SNP names or with columns numbers on train.data[,1] etc.
fit = lm(train.data$CBSDRS ~ train.data[,11])
#m1<-lmer(CBSDRS ~ (1|CLONE) + LOCATION + (REP %in% LOCATION), data=phenos1)


summary(fit) # gives more extensive summary
summary(fit)$coefficients # only gives table with estimates and p-values

# Function to automate GWAS on training data;
# The function is made so that 'x' can be passed as one SNP column and then
# only the line with estimate and p-value is reported.
gwas_CBSDRS = function(x) {
  fit = lm(train.data$CBSDRS ~ x)
  summary(fit)$coefficients[2, ]
}
# test it for instance by re-doing the first SNP
gwas_CBSDRS(train.data[,11])

# Complete GWAS in training data: 
# apply the gwas_CBSDRS function to the whole table of genotypes.
# Note: the apply function organizes the results in columns, but it looks nicer
# in rows, therefore t() (transpose) is applied.
# The results table should be the number rows for all SNPs, and may run a few seconds. 
gwas_CBSDRS_results = t(apply(train.data[,11:1000], 2, gwas_CBSDRS))

# For a Manhattan plot you can plot the -10log(pvalues),
# the p-values are in column 4 of the results
plot(-log10(gwas_CBSDRS_results[, 4]), cex = 0.5)

# Identification of most significant SNPs;
# Start with a low correction 0.05 threshold
snp_select = which(gwas_CBSDRS_results[, 4] < 0.05)

# extract from results table with the significant SNPs
gwas_CBSDRS_results[snp_select,]


# Prediction in test data with GWAS significant SNPs
# We can directly apply the estimates from gwas_CBSDRS_results (column 1) to make
# a prediction by multiplying with the genotypes in the test data.
# We just sum the effects of all significant SNPs.
# Extra: you could try if it makes a difference when using a cleaned-up list
# after use of the multi-SNP model to remove redundant ones.

gebv_test = as.matrix(snps.train[,snp_select]) %*% gwas_CBSDRS_results[snp_select, 1]

cor(train.data$CBSDRS, gebv_test, use = "pair")

# Computing SNP explained variance in training
# For an explained variance per SNP this can be done by computing 2pq(effect^2).
# However, as Heffner warns about QTL mapping, these effects will be (very) over-estimated
SNPp = colMeans(test.data[,11:1000], na.rm = T)/2
snp_select_var = 2 * SNPp[snp_select] * (1 - SNPp[snp_select]) * (gwas_CBSDRS_results[snp_select, 1])^2
# However, summing this up as the explained variance is not a good idea,if there are redundant SNPs




#Genomic relationship matrices
#1 ------- VanRaden method 2 G matrix
# This version is based on computing the SNP means, the SNP frequency 'p' (half the means),
# and the SNP variance as 2p(1-p). snpmean and snpvar should all be vectors of size 1000.
# It should also be possible to use the scale() function specifying scale=sqrt(snpvar).
# NA's are not removed at first; after centering NAs are replaced by zero, which are the
# SNP means after centering.

snps.ALL=snps[genos,]
snpmean = colMeans(snps.ALL,na.rm=T)
SNPp = snpmean/2
snpvar = 2*SNPp*(1-SNPp)
geno_m2 = t((t(snps.ALL)-snpmean)/sqrt(snpvar))
geno_m2[is.na(geno_m2)]=0
G2 = geno_m2 %*% t(geno_m2) / (ncol(snps.ALL))# this divides by the total number of snps

#check distribution of diagonal values and the clustering of groups
hist(diag(G2))
heatmap(G2,symm=T)

# 2 -------- A "modified" VanRaden method 3 G-matrix.
# I call it method 3 because like VanRaden method 3 it is based on a global centering
# and global scaling. The modification is that it doesn't use pedigree to obtain the
# centering and scaling coefficients, but this version takes mean(MM') for the
# global centering sum(2pq) for global scaling.
# A small bit of extra work is to handle the missing genotypes.
# We can't center first and then replace all missings by zero (the mean after centering);
# we have to replace NAs before centering. To do that, we have to go through all columns
# of the SNP data and replace the missings in each column with the SNP-mean of that column.
# This is done in a loop (there may be smarter ways).

snps.ALL=snps[genos,]
snpmeans = colMeans(snps.ALL,na.rm=T)
SNPp = snpmeans/2
for(i in 1:1000){ snps.ALL[is.na(snps.ALL[,i]),i] = snpmeans[i]}
MMtr = snps.ALL %*% t(snps.ALL) # MM' in VanRaden
global_centering = mean(MMtr)
global_scaling = sum(2*SNPp*(1-SNPp))
G3 = (MMtr - global_centering)/global_scaling
hist(diag(G3))
heatmap(G3,symm=T)






