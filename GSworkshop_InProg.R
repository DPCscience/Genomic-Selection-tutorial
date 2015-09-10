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
# test by re-doing the first SNP
gwas_CBSDRS(train.data[,11])

# Complete GWAS in training data: 
# apply the function to the whole table of genotypes.
# The apply function organizes the results in columns in rows, therefore t() (transpose) is applied
# The results table should be the number rows for all SNPs
gwas_CBSDRS_results = t(apply(train.data[,11:1000], 2, gwas_CBSDRS))

# For a Manhattan plot you can plot the -10log(pvalues),
# the p-values are in column 4 of the results, we cant see the chromosome numbers in needed use qqman R package
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
# Version is based on computing SNP means, the SNP frequency 'p' (half the means),
# and the SNP variance as 2p(1-p) = 2pq. snpmean and snpvar should all be vectors of size 1000 (#markers).
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
# Like VanRaden method 3 it is based on a global centering
# and global scaling. The modification is that it doesn't use pedigree to obtain the
# centering and scaling coefficients, but this version takes mean(MM') for the
# global centering sum(2pq) for global scaling.
# we replace NAs before centering. To do that, we have to go through all columns
# of the SNP data and replace the missings in each column with the SNP-mean of that column.
# This is done in a loop 

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


# LOAD DATA different subsetting
setwd("~/Google Drive/TAMU_Workshop_2015/R")
load("phenos.Rdata")
load("snps.Rdata")
phenos1<-subset(phenos,YEAR == "2013")

# CENTERING PREDICTORS (MARKER DOSAGES) ------------------
# If not corrected for, markers with higher allele frequency will contribute differentially to prediction models
#marker dosages subtracted 2P (P = freq. of alt. allele)

frq <- apply(snps, 2, function(x) { mean(x,na.rm=T) })/2
P <- matrix(rep(frq,nrow(snps)),byrow=T,ncol=ncol(snps))
M <- snps-2*P # The rescaled predictors now effectively have mean of zero

# TAKE A SUBSET OF GENOTYPES AS "TRAINING DATA" --------------------------------------------------------------------
genos<-as.character(unique(phenos1$CLONE)) # 408 clones
set.seed(123)
train<-sample(genos,length(genos)*.8) # Set 4/5 as training
test<-setdiff(genos,train) # Set the rest as test

phenos.trn=phenos1[which(phenos1$CLONE %in% train),]
phenos.trn$CLONE<-droplevels(phenos.trn$CLONE) 
M.trn=M[train,]
snps.trn=snps[train,]

phenos.tst=phenos1[which(phenos1$CLONE %in% test),]
phenos.tst$CLONE<-droplevels(phenos.tst$CLONE) 
M.tst=M[test,]
snps.tst=snps[test,]

# RANDOM REGRESSION FUNCTIONS IN THE RRBLUP R PACKAGE ----------------------------------------------------------

library(rrBLUP)

# INCIDENCE MATRICES ----------------------------------------------------------------------------------------
# Replications within and among locations will be included as fixed effects
# Mixed solve requires that we supply the incidence matrices (a.k.a. model matrices) for fixed and random effects
# We will create these with the model.matrix function
# The matrix for fixed effects will be N_obs by N_levels
X=model.matrix(~LOCATION+(factor(REP) %in% LOCATION),data=phenos.trn) # Location effect plus Rep nested in location
# The matrix for random effects will be N_obs x N_snps
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) # -1 because we don't want an intercept here
# Vector of phenotypes
y=phenos.trn$CBSDRS

# RR-BLUP variance in marker effect, calculate the additive genetic variance in GBLUP 
ans.rr <- mixed.solve(y=y,Z=Z%*%M.trn,X=X)
barplot(ans.rr$u^2)
barplot(abs(ans.rr$u))
2*sum(freq*(1-freq))*ans.rr$Vu # 2*sum(pq)*Vu total variance explained by the markers

# G-BLUP (mixed.solve) : use A.mat function to create relationship matrices
M=X
M= X[,1:ncol(X)]-1
p1=round((apply(X,2,sum)+nrow(X))/(nrow(X)*2),3)
K <- (M.trn%*%t(M.trn))/2*(sum(p1*(1-p1)) # OPTION 1 GRM
K <- A.mat(snps.trn-1) # OPTION 2 GRM
Z=model.matrix(~factor(CLONE,levels=rownames(K))-1,data=phenos.trn)
ans.g <- mixed.solve(y=y,Z=Z,K=K,X=X)
ans.g$Vu #differences in blup same as traits, values centered on zero 
#genetic variance explained matches RR-BLUP total variance

# G-BLUP (kin.blup) ---------------------------------------------------------------------------------------------------
# Wrapper for mixed.solve
phenos2<-phenos1
phenos2[which(phenos2$CLONE %in% phenos.tst$CLONE),6:9]<-NA
phenos2$LOC.REP<-paste(phenos2$LOCATION,phenos2$REP,sep=".")
ans.g.kin<-kin.blup(data=phenos2,geno="CLONE",pheno="CBSDRS",K=A.mat(snps-1),fixed=c("LOCATION","LOC.REP"))

# Comparison RR-BLUP and GBLUP -------------------------------------------------------------------------------------------
pred.rr=(M.trn%*%as.matrix(ans.rr$u))
pred.g=ans.g$u
plot(pred.rr[,1],pred.g,xlab="RR-BLUP Prediction",ylab="GBLUP Prediction")
abline(a=0,b=1,col='red')

# PREDICTING UN-PHENOTYPED INDIVIDUALS --------------------------------------------------------------------------------
# With rr-BLUP mixed solve you obtain marker effects, which can be used to predict any indidual genotyped with those markeres
# With GBLUP you obtain a GEBV for every individual represented in the kinship matrix used, with or without phenotypes
# When we fit the GBLUP model with mixed.solve, we used only marker data for training set to create a kinship matrix (ans.g$u)
# Use full SNP matrix with kin.blup, to obtain GEBV (ans.g.kin$g) for all genotyped individuals

# PREDICTION ACCURACY --------------------------------------------------------------------------------------------------
# Correlation between GEBV and raw phenotypes
# Correlation between GEBV and EGV (estimated genotypic value)
# EGV is BLUP obtained by fitting the same mixed model but without a kinship matrix (K = I, levels of genotype are i.i.d.) 
# prediction of validation-test sets
pred.tst.rr=as.matrix(M.tst%*%as.matrix(ans.rr$u))
pred.tst.g=as.matrix(ans.g.kin$g[test])

# Correlate prediction with raw data
phenos.tst1<-merge(phenos.tst[,c("CLONE","CBSDRS")],pred.tst.g,by.x="CLONE",by.y="row.names")
plot(phenos.tst1$V1+mean(phenos.tst1$CBSDRS,na.rm=T),phenos.tst1$CBSDRS,xlab="GEBV",ylab="Raw CBSDRS Score")
abline(a=0,b=1,col='red')
cor(phenos.tst1$CBSDRS,phenos.tst1$V1,use='complete.obs')

# Correlate prediction with estimated genetic value (EGV)
library(lme4)
phenos1$REP<-factor(phenos1$REP)
m1<-lmer(CBSDRS ~ (1|CLONE) + LOCATION + (REP %in% LOCATION), data=phenos1)
egv<-as.matrix(ranef(m1)$CLONE)
egv.tst<-as.matrix(egv[rownames(egv) %in% test,])
tst=merge(egv.tst,pred.tst.g,by="row.names")

plot(tst[,3],tst[,2],ylab="Estimated Genetic Value",xlab="GEBV (GBLUP)")
abline(a=0,b=1,col='red')
cor(tst[,2],tst[,3])

# BAYESIAN PREDICTION TOOLS (BGLR Package) -----------------------------------------------------------------------------
library(BGLR)

# Bayesian RR-BLUP 
# For BGLR, first make a list with the model elements to fit.
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) 
ETA.rr = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BRR"))
fit.rr = BGLR(y=y,ETA=ETA.rr)
# Output from BGLR is complicated. You can look at the contents with:
str(fit.rr)
fit.rr$varE # residual variance
# The list object ETA has one level for each term in the ETA list originally specified 
# (one for fixed, one for random/markers in this case)
# Let's look at the SNP effects specifially
fit.rr$ETA[[2]]$varB    # the SNP variance
barplot(fit.rr$ETA[[2]]$b^2)      # the SNP effects
# Compare marker effects from Bayesian vs. Frequentist RR-BLUP
plot(fit.rr$ETA[[2]]$b,ans.rr$u)
# Accuracy of predicting test set
pred.tst.brr=as.matrix(M.tst%*%as.matrix(fit.rr$ETA[[2]]$b))
tst=merge(egv.tst,pred.tst.brr,by="row.names")
plot(tst[,3],tst[,2],ylab="Estimated Genetic Value",xlab="GEBV (Bayesian RR-BLUP)")
abline(a=0,b=1,col='red')
cor(tst[,2],tst[,3]) # Accuracy is almost identical to GBLUP / RR-BLUP

# BayesA
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) 
ETA.BayesA = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BayesA"))
fit.BayesA = BGLR(y=y,ETA=ETA.BayesA)
str(fit.BayesA)

barplot(fit.BayesA$ETA[[2]]$b^2)
plot(fit.BayesA$ETA[[2]]$b,ans.rr$u)
abline(a=0,b=1,col='red')

pred.tst.BayesA=as.matrix(M.tst%*%as.matrix(fit.BayesA$ETA[[2]]$b))
tst=merge(egv.tst,pred.tst.BayesA,by="row.names")
cor(tst[,2],tst[,3]) 

# BayesB
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) 
ETA.BayesB = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BayesB"))
fit.BayesB = BGLR(y=y,ETA=ETA.BayesB)
str(fit.BayesB)

barplot(fit.BayesB$ETA[[2]]$b^2)
plot(fit.BayesB$ETA[[2]]$b,ans.rr$u)
abline(a=0,b=1,col='red')

pred.tst.BayesB=as.matrix(M.tst%*%as.matrix(fit.BayesB$ETA[[2]]$b))
tst=merge(egv.tst,pred.tst.BayesB,by="row.names")
cor(tst[,2],tst[,3])

# Bayesian LASSO
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) 
ETA.BL = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BL"))
fit.BL = BGLR(y=y,ETA=ETA.BL)
str(fit.BL)

barplot(fit.BL$ETA[[2]]$b^2)
plot(fit.BL$ETA[[2]]$b,ans.rr$u)
abline(a=0,b=1,col='red')

pred.tst.BL=as.matrix(M.tst%*%as.matrix(fit.BL$ETA[[2]]$b))
tst=merge(egv.tst,pred.tst.BL,by="row.names")
cor(tst[,2],tst[,3])

# Using more iterations / burnIn to solve Bayesian models
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) 
ETA.BL = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BL"))
fit.BL1 = BGLR(y=y,ETA=ETA.BL,nIter=6000,burnIn=1000)
str(fit.BL1)

barplot(fit.BL1$ETA[[2]]$b^2)
plot(fit.BL1$ETA[[2]]$b,ans.rr$u)
abline(a=0,b=1,col='red')

pred.tst.BL1=as.matrix(M.tst%*%as.matrix(fit.BL1$ETA[[2]]$b))
tst=merge(egv.tst,pred.tst.BL1,by="row.names")
cor(tst[,2],tst[,3]) # More iterations did improve accuracy

# Flexibility of BGRL in response type
# CMD3S is a disease severity score 1-5 (it is an ordinal variable with a distribution skewed to 1 [resistant])
y=phenos.trn$CMD3S

# BayesB Ordinal Response for CMD3S
Z=model.matrix(~factor(CLONE,levels=rownames(M.trn))-1,data=phenos.trn) 
ETA.BBord = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BayesB"))
fit.BBord = BGLR(y=y,ETA=ETA.BayesB,response_type="ordinal")
barplot(fit.BBord$ETA[[2]]$b^2)
# BayesB Continuous Resposne for CMD3S
ETA.BB = list(list(X = X, model = "FIXED"), list(X = Z%*%M.trn, model = "BayesB"))
fit.BB = BGLR(y=y,ETA=ETA.BayesB)
barplot(fit.BB$ETA[[2]]$b^2)
# GEBV from ordinal vs. continuous model
pred.tst.BBord=as.matrix(M.tst%*%as.matrix(fit.BBord$ETA[[2]]$b))
pred.tst.BB=as.matrix(M.tst%*%as.matrix(fit.BB$ETA[[2]]$b))
# EGV of CMD3S (continuous)
m2<-lmer(CMD3S ~ (1|CLONE) + LOCATION + (REP %in% LOCATION), data=phenos1)
egv.cmd<-as.matrix(ranef(m2)$CLONE)
egv.cmd.tst<-as.matrix(egv.cmd[rownames(egv.cmd) %in% test,])
# Accuracy of ordinal
tst=merge(egv.cmd.tst,pred.tst.BBord,by="row.names")
cor(tst[,2],tst[,3])
# Accuracy of continuous
tst=merge(egv.cmd.tst,pred.tst.BB,by="row.names")
cor(tst[,2],tst[,3])



