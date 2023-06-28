#FAST modelling - Maxent on DISMO

library(dismo)
library(raster)
library(rJava)
library(rgdal)
library(virtualspecies)
library(maptools)

#Load CSV (sp, lat, lon)

read.csv(proc)

#Remover duplicados

data(wrld_simpl)

uni <- unique(proc[, 2:3])
proc <- uni

dups <- duplicated(proc[, 1:2])
dups
sum(dups)

rm(dups)
rm(uni)

#Georreferenciar

coordinates(proc) <- ~longitude+latitude
crs(proc) <- crs(wrld_simpl)
class(proc)

#Criar shape de corte

e <- as(extent(-50, -34, -20, 0), "SpatialPolygons")
proc <- crop(proc, e)
class(proc)

#Localizar os dados ambientais e criar preditores, os arquivos bio devem estar dentro da pasta

setwd("source directory")
files <- list.files(path=".", pattern='tif', full.names=TRUE)
predictors <- stack(files)

#área de estudo

crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
predictors <- crop (predictors, e)
gc()

# Co-linearity between  predictors 

predictors.red <- removeCollinearity(predictors, multicollinearity.cutoff = 0.85,
                                     select.variables = FALSE, sample.points = TRUE, nb.points =1000,  plot = TRUE, method = "spearman")

#predictors.red 

gc()

# use this info in the next command line

rasters.selected <- subset(predictors, c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_18","wc2.1_30s_bio_1","wc2.1_30s_bio_19","wc2.1_30s_bio_12", 
                                           "wc2.1_30s_bio_15", "wc2.1_30s_bio_14", "wc2.1_30s_bio_6", "wc2.1_30s_bio_2", "wc2.1_30s_bio_8", "wc2.1_30s_bio_5"))


predictors <-rasters.selected
rm(rasters.selected)
rm(files)
gc()

# reduzir dos dados de proc (um ponto por célula de ocorrência)
cells <- cellFromXY(predictors[[1]], proc)
dups <- duplicated(cells)
procun <- proc[!dups, ]
cat(nrow(procun) - nrow(proc), "records are removed")

proc <- procun

rm(procun)
rm(files)
rm(cells)
rm(dups)
gc()

#Criar um poligono como area de distribuicao

r <- raster(proc)
res(r) <- 1

procsel <- gridSample(proc, r, n= .5)

#Background

bg <- sampleRandom(x=predictors,
                   size=20000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 



#Criar arquivos de treino e teste

fold <- kfold(procsel, k=2)
proctest <- procsel[fold == 1, ]
proctrain <- procsel[fold != 1, ]


#modelo maxent

procfitm <- maxent(predictors, proctrain)
modelm <- predict(procfitm, predictors, progress='text')
plot(modelm)
e1mcc = evaluate(procfitm, p=proctest, a=bg, x=predictors)
pvtest <- data.frame(extract(predictors, proctest))
avtest <- data.frame(extract(predictors, bg))
e2mcc = evaluate(procfitm, p=pvtest, a=avtest)
testp <- predict(procfitm, pvtest)
testa <- predict(procfitm, avtest)
e3mcc = evaluate(p=testp, a=testa)
e3mcc
writeRaster(modelm, filename="Pleroma_catingae.tif", overwrite=TRUE, bylayer=TRUE, format="GTiff")



#Removing bias

sb[,1] / sb[,2]

i <- pwdSample(proctest, bg, proctrain, n=1, tr=1)
pres_test_pwd <- proctest[!is.na(i[,1]), ]
back_test_pwd <- bg[na.omit(as.vector(i)), ]
sb2 <- dismo::ssb(pres_test_pwd, back_test_pwd, proctrain)
sb2[1]/ sb2[2]



e4mcc <- evaluate(procfitm, p=pres_test_pwd, a=back_test_pwd, x=predictors)

