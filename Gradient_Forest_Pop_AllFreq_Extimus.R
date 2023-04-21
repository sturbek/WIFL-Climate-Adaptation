
library(gradientForest)
library(data.table)
library(raster)
library(gstat)
library(rgdal)
library(sf)
library(RColorBrewer)
library(rasterVis)
library(tidyverse)

Ewifl<-read_csv("WIFL_contemp_sites_extimus_environmental_data_1996_2015.csv")
dim(Ewifl)
# 11 sites, 9 variables

Ewifl2<-Ewifl[,-c(1,2)]
# 11 sites, no lat, lon

Gwifl<-read_csv("WIFL_contemp_samples_allele_frequencies_extimus.csv",col_names=T)
Gwifl <- Gwifl[-1,]
dim(Gwifl)
# 11 sites, 3,757,986 loci

Gwifl2 <- Gwifl
Gwifl2 <- Gwifl2[ ,colSums(is.na(Gwifl2)) == 0]
dim(Gwifl2)
# 11 sites, 2,274,783 loci

Gwifl2 <- data.frame(lapply(Gwifl2, function(x) as.numeric(as.character(x))),
           check.names=F, row.names = rownames(Gwifl2))
Gwifl2 <- tibble(Gwifl2)

# take random subset of 500,000 loci
Gwifl3 <- Gwifl2[, sample(ncol(Gwifl2), 500000) ]
dim(Gwifl3)

preds <- colnames(Ewifl2)
specs <- colnames(Gwifl3)

nSites <- dim(Gwifl3)[1]
nSpecs <- dim(Gwifl3)[2]

# set depth of conditional permutation
lev <- floor(log2(nSites*0.368/2))
lev

wiflforest=gradientForest(cbind(Ewifl2,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 1
Ewifl_R1 <- Ewifl2
Ewifl_R1$precip <- sample(Ewifl_R1$precip)
Ewifl_R1$tmax <- sample(Ewifl_R1$tmax)
Ewifl_R1$tmin <- sample(Ewifl_R1$tmin)
Ewifl_R1$tmean <- sample(Ewifl_R1$tmean)
Ewifl_R1$tdmean <- sample(Ewifl_R1$tdmean)
Ewifl_R1$vpdmax <- sample(Ewifl_R1$vpdmax)
Ewifl_R1$vpdmin <- sample(Ewifl_R1$vpdmin)
wiflforest_R1=gradientForest(cbind(Ewifl_R1,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 2
Ewifl_R2 <- Ewifl2
Ewifl_R2$precip <- sample(Ewifl_R2$precip)
Ewifl_R2$tmax <- sample(Ewifl_R2$tmax)
Ewifl_R2$tmin <- sample(Ewifl_R2$tmin)
Ewifl_R2$tmean <- sample(Ewifl_R2$tmean)
Ewifl_R2$tdmean <- sample(Ewifl_R2$tdmean)
Ewifl_R2$vpdmax <- sample(Ewifl_R2$vpdmax)
Ewifl_R2$vpdmin <- sample(Ewifl_R2$vpdmin)
wiflforest_R2=gradientForest(cbind(Ewifl_R2,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 3
Ewifl_R3 <- Ewifl2
Ewifl_R3$precip <- sample(Ewifl_R3$precip)
Ewifl_R3$tmax <- sample(Ewifl_R3$tmax)
Ewifl_R3$tmin <- sample(Ewifl_R3$tmin)
Ewifl_R3$tmean <- sample(Ewifl_R3$tmean)
Ewifl_R3$tdmean <- sample(Ewifl_R3$tdmean)
Ewifl_R3$vpdmax <- sample(Ewifl_R3$vpdmax)
Ewifl_R3$vpdmin <- sample(Ewifl_R3$vpdmin)
wiflforest_R3=gradientForest(cbind(Ewifl_R3,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 4
Ewifl_R4 <- Ewifl2
Ewifl_R4$precip <- sample(Ewifl_R4$precip)
Ewifl_R4$tmax <- sample(Ewifl_R4$tmax)
Ewifl_R4$tmin <- sample(Ewifl_R4$tmin)
Ewifl_R4$tmean <- sample(Ewifl_R4$tmean)
Ewifl_R4$tdmean <- sample(Ewifl_R4$tdmean)
Ewifl_R4$vpdmax <- sample(Ewifl_R4$vpdmax)
Ewifl_R4$vpdmin <- sample(Ewifl_R4$vpdmin)
wiflforest_R4=gradientForest(cbind(Ewifl_R4,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 5
Ewifl_R5 <- Ewifl2
Ewifl_R5$precip <- sample(Ewifl_R5$precip)
Ewifl_R5$tmax <- sample(Ewifl_R5$tmax)
Ewifl_R5$tmin <- sample(Ewifl_R5$tmin)
Ewifl_R5$tmean <- sample(Ewifl_R5$tmean)
Ewifl_R5$tdmean <- sample(Ewifl_R5$tdmean)
Ewifl_R5$vpdmax <- sample(Ewifl_R5$vpdmax)
Ewifl_R5$vpdmin <- sample(Ewifl_R5$vpdmin)
wiflforest_R5=gradientForest(cbind(Ewifl_R5,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 6
Ewifl_R6 <- Ewifl2
Ewifl_R6$precip <- sample(Ewifl_R6$precip)
Ewifl_R6$tmax <- sample(Ewifl_R6$tmax)
Ewifl_R6$tmin <- sample(Ewifl_R6$tmin)
Ewifl_R6$tmean <- sample(Ewifl_R6$tmean)
Ewifl_R6$tdmean <- sample(Ewifl_R6$tdmean)
Ewifl_R6$vpdmax <- sample(Ewifl_R6$vpdmax)
Ewifl_R6$vpdmin <- sample(Ewifl_R6$vpdmin)
wiflforest_R6 <- readRDS(file = "wiflforest.R6.extimus.ntree500.rds")

# random 7
Ewifl_R7 <- Ewifl2
Ewifl_R7$precip <- sample(Ewifl_R7$precip)
Ewifl_R7$tmax <- sample(Ewifl_R7$tmax)
Ewifl_R7$tmin <- sample(Ewifl_R7$tmin)
Ewifl_R7$tmean <- sample(Ewifl_R7$tmean)
Ewifl_R7$tdmean <- sample(Ewifl_R7$tdmean)
Ewifl_R7$vpdmax <- sample(Ewifl_R7$vpdmax)
Ewifl_R7$vpdmin <- sample(Ewifl_R7$vpdmin)
wiflforest_R7=gradientForest(cbind(Ewifl_R7,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 8
Ewifl_R8 <- Ewifl2
Ewifl_R8$precip <- sample(Ewifl_R8$precip)
Ewifl_R8$tmax <- sample(Ewifl_R8$tmax)
Ewifl_R8$tmin <- sample(Ewifl_R8$tmin)
Ewifl_R8$tmean <- sample(Ewifl_R8$tmean)
Ewifl_R8$tdmean <- sample(Ewifl_R8$tdmean)
Ewifl_R8$vpdmax <- sample(Ewifl_R8$vpdmax)
Ewifl_R8$vpdmin <- sample(Ewifl_R8$vpdmin)
wiflforest_R8=gradientForest(cbind(Ewifl_R8,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 9
Ewifl_R9 <- Ewifl2
Ewifl_R9$precip <- sample(Ewifl_R9$precip)
Ewifl_R9$tmax <- sample(Ewifl_R9$tmax)
Ewifl_R9$tmin <- sample(Ewifl_R9$tmin)
Ewifl_R9$tmean <- sample(Ewifl_R9$tmean)
Ewifl_R9$tdmean <- sample(Ewifl_R9$tdmean)
Ewifl_R9$vpdmax <- sample(Ewifl_R9$vpdmax)
Ewifl_R9$vpdmin <- sample(Ewifl_R9$vpdmin)
wiflforest_R9=gradientForest(cbind(Ewifl_R9,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# random 10
Ewifl_R10 <- Ewifl2
Ewifl_R10$precip <- sample(Ewifl_R10$precip)
Ewifl_R10$tmax <- sample(Ewifl_R10$tmax)
Ewifl_R10$tmin <- sample(Ewifl_R10$tmin)
Ewifl_R10$tmean <- sample(Ewifl_R10$tmean)
Ewifl_R10$tdmean <- sample(Ewifl_R10$tdmean)
Ewifl_R10$vpdmax <- sample(Ewifl_R10$vpdmax)
Ewifl_R10$vpdmin <- sample(Ewifl_R10$vpdmin)
wiflforest_R10=gradientForest(cbind(Ewifl_R10,Gwifl3), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=T,nbin=101, maxLevel=lev,trace=T)

# summmary stats for randomization
wiflforest$species.pos.rsq # 89,562
mean(wiflforest$res$rsq) # 0.1884983

wiflforest_R1$species.pos.rsq # 57,660
mean(wiflforest_R1$res$rsq) # 0.1236358

wiflforest_R2$species.pos.rsq # 55,141
mean(wiflforest_R2$res$rsq) # 0.1172123

wiflforest_R3$species.pos.rsq # 54,164
mean(wiflforest_R3$res$rsq) # 0.1089199

wiflforest_R4$species.pos.rsq # 49,974
mean(wiflforest_R4$res$rsq) # 0.1092889

wiflforest_R5$species.pos.rsq # 60,329
mean(wiflforest_R5$res$rsq) # 0.1272331

wiflforest_R6$species.pos.rsq # 41,154
mean(wiflforest_R6$res$rsq) # 0.1103168

wiflforest_R7$species.pos.rsq # 68,644
mean(wiflforest_R7$res$rsq) # 0.1188598

wiflforest_R8$species.pos.rsq # 47,693
mean(wiflforest_R8$res$rsq) # 0.09366692

wiflforest_R9$species.pos.rsq # 45,525
mean(wiflforest_R9$res$rsq) # 0.1027376

wiflforest_R10$species.pos.rsq # 56,231
mean(wiflforest_R10$res$rsq) # 0.1107559

plot(wiflforest,plot.type="Overall.Importance")
plot(wiflforest,plot.type="Cumulative.Importance")
