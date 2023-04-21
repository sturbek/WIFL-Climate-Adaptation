
library(LEA)
library(tidyverse)

# convert ped file to lfmm
ped2geno("WIFL_contemp_samples_221_imputed_gl.ped")
geno2lfmm("WIFL_contemp_samples_221_imputed_gl.geno")
# 221 individuals 128147 loci

# read in lfmm file and environmental predictors
gen.imp <- read.lfmm("WIFL_contemp_samples_221_imputed_gl.lfmm")
dim(gen.imp)

pred <- read.table("WIFL_contemp_samples_221_top3_environmental_data.txt",header=F)
dim(pred)

project = load.lfmmProject("WIFL_contemp_samples_221_imputed_gl.lfmmProject")

### VARIABLE 1 (Precip) ###
#Record z-scores from the 5 runs in the zs matrix
zs.d1 = z.scores(project, K = 3, d = 1)
#Combine z-scores using the median
zs.median.d1 = apply(zs.d1, MARGIN = 1, median)

# compute adjusted p-values from the combined z-scores
lambda = median(zs.median.d1^2) / 0.456
lambda # 0.989096
adj.p.values = pchisq(zs.median.d1^2 / lambda, df = 1, lower = F)

#histogram of p-values
hist(adj.p.values, col = "red")

pdf(file = "Histogram_AdjustedPvalues_Precip.pdf",width = 8,height = 4) 
hist(adj.p.values, col = "red")
dev.off()

## FDR control: Benjamini-Hochberg at level q
L = 128147 # number of loci
q = 0.1 # fdr level
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
length(candidates.bh) # 536

write.table(candidates.bh,file=paste("WIFL_contemp_samples_221_imputed_gl_precip_candidates.csv",sep=""),row.names=FALSE)

### VARIABLE 2 (tmax) ###
#Record z-scores from the 5 runs in the zs matrix
zs.d2 = z.scores(project, K = 3, d = 2)
#Combine z-scores using the median
zs.median.d2 = apply(zs.d2, MARGIN = 1, median)

# compute adjusted p-values from the combined z-scores
lambda = median(zs.median.d2^2) / 0.456
lambda # 1.04
adj.p.values = pchisq(zs.median.d2^2 / lambda, df = 1, lower = F)

#histogram of p-values
hist(adj.p.values, col = "red")

pdf(file = "Histogram_AdjustedPvalues_Tmax.pdf",width = 8,height = 4) 
hist(adj.p.values, col = "red")
dev.off()

## FDR control: Benjamini-Hochberg at level q
L = 128147 # number of loci
q = 0.1 # fdr level
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
length(candidates.bh) # 1120

write.table(candidates.bh,file=paste("WIFL_contemp_samples_221_imputed_gl_tmax_candidates.csv",sep=""),row.names=FALSE)

### VARIABLE 3 (tdmean) ###
#Record z-scores from the 5 runs in the zs matrix
zs.d3 = z.scores(project, K = 3, d = 3)
#Combine z-scores using the median
zs.median.d3 = apply(zs.d3, MARGIN = 1, median)

# compute adjusted p-values from the combined z-scores
lambda = median(zs.median.d3^2) / 0.456
lambda # 1.198078
adj.p.values = pchisq(zs.median.d3^2 / lambda, df = 1, lower = F)

#histogram of p-values
hist(adj.p.values, col = "red")

pdf(file = "Histogram_AdjustedPvalues_Tdmean.pdf",width = 8,height = 4) 
hist(adj.p.values, col = "red")
dev.off()

## FDR control: Benjamini-Hochberg at level q
L = 128147 # number of loci
q = 0.1 # fdr level
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
length(candidates.bh) # 745

write.table(candidates.bh,file=paste("WIFL_contemp_samples_221_imputed_gl_tdmean_candidates.csv",sep=""),row.names=FALSE)

