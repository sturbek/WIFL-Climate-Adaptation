
library(vegan)
library(tidyverse)
library(LEA)
library(psych)

gen.imp = read.lfmm("WIFL_contemp_samples_221_imputed_gl.lfmm")
dim(gen.imp) # # 221 individuals, 128,147 loci

pred = read.table("WIFL_contemp_samples_221_top3_environmental_data.txt",header=T) # precip tmax tdmean pc1
dim(pred)

cor(pred,method="pearson")
#           precip       tmax     tdmean        pc1
# precip  1.0000000 -0.2751677  0.6155092 -0.3173534
# tmax   -0.2751677  1.0000000 -0.2476562  0.5180017
# tdmean  0.6155092 -0.2476562  1.0000000 -0.4739768
# pc1    -0.3173534  0.5180017 -0.4739768  1.0000000

WIFL.rda <- rda(gen.imp ~ precip + tmax + tdmean + Condition(pc1), data=pred, scale=T)

RsquareAdj(WIFL.rda)

summary(eigenvals(WIFL.rda, model = "constrained"))
# Importance of components:
#                          RDA1     RDA2     RDA3
# Eigenvalue            1170.3784 787.2025 644.0944
# Proportion Explained     0.4499   0.3026   0.2476
# Cumulative Proportion    0.4499   0.7524   1.0000

screeplot(WIFL.rda)

signif.full <- anova.cca(WIFL.rda, parallel=getOption("mc.cores"))
# Model: rda(formula = gen.imp ~ precip + tmax + tdmean + Condition(pc1), data = pred, scale = T)
# Df Variance      F Pr(>F)    
# Model      3     2602 1.5534  0.001 ***
# Residual 216   120591

signif.axis.tmax <- anova.cca(WIFL.rda, by="axis", parallel=getOption("mc.cores"))
# Model: rda(formula = gen.imp ~ precip + tmax + tdmean + Condition(pc1), data = pred, scale = T)
# Df Variance      F Pr(>F)    
# RDA1       1     1170 2.0964  0.001 ***
# RDA2       1      787 1.4100  0.001 ***
# RDA3       1      644 1.1537  0.001 ***
# Residual 216   120591                  

vif.cca(WIFL.rda)
# pc1   precip     tmax   tdmean 
# 1.657880 1.655861 1.403825 1.884750 

intersetcor(WIFL.rda)[,1:3] # axis 1 associated with tmax, axis 2 associated with tdmean, axis 3 associated with precip
#           RDA1       RDA2       RDA3
# precip 0.1562985  0.3941304 -0.8611545
# tmax   0.7782275 -0.2061215  0.2447087
# tdmean 0.2075102  0.9115485 -0.1207737


load.rda <- scores(WIFL.rda, choices=c(1:4), display="species")  # Species scores for the first three constrained axes

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3)
length(cand1)  # 245
cand2 <- outliers(load.rda[,2],3)
length(cand2) # 344
cand3 <- outliers(load.rda[,3],3)
length(cand3) # 392

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand # 981 SNPs

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=3)  # 3 columns for 4 predictors
colnames(foo) <- c("precip","tmax","tdmean")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred[,-4],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  # 9 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # 0 duplicates on axis 1
table(foo[foo[,1]==2,2]) # 2 duplicates on axis 2
table(foo[foo[,1]==3,2]) # 7 duplicates on axis 3

cand <- cand[!duplicated(cand$snp),]
length(cand$snp) # 972 SNPs

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
}

colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)
# precip tdmean   tmax 
# 338    340    294 

write.csv(cand,"WIFL.RDA.candidates.csv",row.names = F,quote=F)

cand_precip <- subset(cand, predictor == "precip")
cand_tmax <- subset(cand, predictor == "tmax")
cand_tdmean <- subset(cand, predictor == "tdmean")

cand_precip$snp<-gsub("V","",as.character(cand_precip$snp))
cand_tmax$snp<-gsub("V","",as.character(cand_tmax$snp))
cand_tdmean$snp<-gsub("V","",as.character(cand_tdmean$snp))

cand_precip_lfmm <- read.csv("WIFL_contemp_samples_221_imputed_gl_precip_candidates.csv",header=T)
cand_precip_lfmm$x <- as.character(cand_precip_lfmm$x)
cand_tmax_lfmm <- read.csv("WIFL_contemp_samples_221_imputed_gl_tmax_candidates.csv",header=T)
cand_tmax_lfmm$x <- as.character(cand_tmax_lfmm$x)
cand_tdmean_lfmm <- read.csv("WIFL_contemp_samples_221_imputed_gl_tdmean_candidates.csv",header=T)
cand_tdmean_lfmm$x <- as.character(cand_tdmean_lfmm$x)

precip_intersect <- intersect(cand_precip$snp, cand_precip_lfmm$x) ## found by both LFMM and RDA
length(precip_intersect) # 72 snps

tmax_intersect <- intersect(cand_tmax$snp, cand_tmax_lfmm$x) ## found by both LFMM and RDA
length(tmax_intersect) # 56 snps

tdmean_intersect <- intersect(cand_tdmean$snp, cand_tdmean_lfmm$x) ## found by both LFMM and RDA
length(tdmean_intersect) # 104 snps
