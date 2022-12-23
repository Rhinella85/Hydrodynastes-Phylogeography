#Load necessary libraries
library(raster)
library(rworldmap)
library(rgdal)

#Input your spp locality data
DataSpecies=read.csv(file="bicinctus_data_nucDNA.csv")
#Plot the extent of your study region
plot(getMap(), xlim = c(-86, -33), ylim = c(-56, 14), asp=1)
#Plot your species localities
points(DataSpecies$long,DataSpecies$lat,pch=20,col="red")

#Load the WorldClim data
bio1=raster('bio_1.tif')
bio2=raster("bio_2.tif") 
bio3=raster("bio_3.tif") 
bio4=raster("bio_4.tif")
bio5=raster("bio_5.tif")
bio6=raster("bio_6.tif")
bio7=raster("bio_7.tif")
bio8=raster("bio_8.tif") 
bio9=raster("bio_9.tif")
bio10=raster("bio_10.tif")
bio11=raster("bio_11.tif")
bio12=raster("bio_12.tif")
bio13=raster("bio_13.tif")
bio14=raster("bio_14.tif")
bio15=raster("bio_15.tif")
bio16=raster("bio_16.tif")
bio17=raster("bio_17.tif")
bio18=raster("bio_18.tif")
bio19=raster("bio_19.tif")

#Below makes the bioclim rasters smaller in extent, perhaps easier to work but might also not matter...
#Define the extent of your study region, change this to wherever you're study region is.
ext_user=extent(c(-86, -33, -56, 14))
#Crop the climate layers
AnMeanTemp=crop(bio1,ext_user)
MeanDRange=crop(bio2,ext_user)
Iso=crop(bio3,ext_user)
TempSeaso=crop(bio4,ext_user)
MaxTempWarmMonth=crop(bio5,ext_user)
MintempColdMonth=crop(bio6,ext_user)
TempeAnRange=crop(bio7,ext_user)
MeanTempWetQuarter=crop(bio8,ext_user)
MeanTempDriQuarter=crop(bio9,ext_user)
MeanTempWarmQuarter=crop(bio10,ext_user)
MeanTempColdQuarter=crop(bio11,ext_user)
AnPreci=crop(bio12,ext_user)
PreciWetMonth=crop(bio13,ext_user)
PreciDriMonth=crop(bio14,ext_user)
PreciSeaso=crop(bio15,ext_user)
PreciWetQuarter=crop(bio16,ext_user)
PreciDriQuarter=crop(bio17,ext_user)
PreciWarmQuarter=crop(bio18,ext_user)
PreciColdQuarter=crop(bio19,ext_user)

#Make sure you did it right
plot(PreciSeaso)
#Add your taxon data if you want
points(DataSpecies$long,DataSpecies$lat,pch=20,col="red")

#Extract the clim data from your locality data
newbio1 <- subset((raster::extract(AnMeanTemp, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_1)
newbio2 <- subset((raster::extract(MeanDRange, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_2)
newbio3 <- subset((raster::extract(Iso, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_3)
newbio4 <- subset((raster::extract(TempSeaso, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_4)
newbio5 <- subset((raster::extract(MaxTempWarmMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_5)
newbio6 <- subset((raster::extract(MintempColdMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_6)
newbio7 <- subset((raster::extract(TempeAnRange, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_7)
newbio8 <- subset((raster::extract(MeanTempWetQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_8)
newbio9 <- subset((raster::extract(MeanTempDriQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_9)
newbio10 <- subset((raster::extract(MeanTempWarmQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_10)
newbio11 <- subset((raster::extract(MeanTempColdQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_11)
newbio12 <- subset((raster::extract(AnPreci, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_12)
newbio13 <- subset((raster::extract(PreciWetMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_13)
newbio14 <- subset((raster::extract(PreciDriMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_14)
newbio15 <- subset((raster::extract(PreciSeaso, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_15)
newbio16 <- subset((raster::extract(PreciWetQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_16)
newbio17 <- subset((raster::extract(PreciDriQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_17)
newbio18 <- subset((raster::extract(PreciWarmQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_18)
newbio19 <- subset((raster::extract(PreciColdQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio_19)
#View the first few rows of the vector you just made
head(newbio1)

#Remove correlated variables
#rm(newbio3, newbio7)

#Bind these together with your original locality file
gigas_clim_data_extracted<-cbind(DataSpecies, newbio1, newbio2, newbio4, newbio5, newbio6, newbio8, newbio9, newbio10, newbio11, newbio12, newbio13, newbio14, newbio15, newbio16, newbio17, newbio18, newbio19)
head(gigas_clim_data_extracted)

#Write your finished table to save these data!
write.table(gigas_clim_data_extracted, file= "gigas_clim_data_extracted.txt")



#Now we can do the RDA analysis
#First call 'ape', 'vegan', and get your genetic data
library(ape)
library(vegan)
#Load DNA
dna <- read.dna(file = "H_bicinctus_nuDNA.fasta", format = "fasta")
#Generate a genetic distance matrix, may want to use a different model of substitution
#"GG95" was used for mtDNA
#"JC69" for the nuDNA
gigasDist <- dist.dna(dna, model = "JC69", pairwise.deletion = TRUE)
#for nucDNA
gigasDist[is.na(gigasDist)] <- 0
#PCoA for gen dist matrix
as.dist(scale(gigasDist, F))->gigasMatrix
gigasMatrix[is.na(gigasMatrix)] <- 0
#pcoa((as.dist(scale(gigasDist, F))), correction = "none", rn= NULL)->PCOA_gigas
pcoa(gigasMatrix, correction = "none", rn= NULL)->PCOA_gigas

#Run Full Model with all variables
rdaFullModel <- rda(PCOA_gigas$vectors~long+lat+bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19+locale, data=gigas_clim_data_extracted)
#Summarize this run
#summary(rdaFullModel)
#Get Rsquared ajusted values
RsquareAdj(rdaFullModel)
#Conduct the rda anova
anova(rdaFullModel)
#Plot your results
#plot(rdaFullModel)


#Condition out Geography and Locale so that you are only looking at climate data
rdaClim <- rda(PCOA_gigas$vectors~bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19  + Condition(lat + long + locale), data=gigas_clim_data_extracted)

#Summarize this run
#summary(rdaClim)
#Get Rsquared ajusted values
RsquareAdj(rdaClim)
#Conduct the rda anova
anova(rdaClim)
#Plot your results
#plot(rdaClim)


#Condition out Climate Variables and Locale so you are only looking at the distance between samples
rdaDistance <- rda(PCOA_gigas$vectors~ lat + long + Condition(bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19+locale), data=gigas_clim_data_extracted)

#Summarize this run
#summary(rdaDistance)
#Get Rsquared ajusted values
RsquareAdj(rdaDistance)
#Conduct the rda anova
anova(rdaDistance)
#Plot your results
#plot(rdaDistance)

#Condition out Geography and Climate Variables so you are only looking at locality with respect to phylogeo barrier
rdaLocale <- rda(PCOA_gigas$vectors~ locale + Condition(lat+long+bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19), data=gigas_clim_data_extracted)

#Summarize this run
#summary(rdaLocale)
#Get Rsquared ajusted values
RsquareAdj(rdaLocale)
#Conduct the rda anova
anova(rdaLocale)
#Plot your results
#plot(rdaLocale)

#Condition out Climate so that you are only looking at both Geography and Locale
rdaGeoLoc <- rda(PCOA_gigas$vectors~lat + long + locale + Condition(bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19), data=gigas_clim_data_extracted)

#Summarize this run
#summary(rdaGeoLoc)
#Get Rsquared ajusted values
RsquareAdj(rdaGeoLoc)
#Conduct the rda anova
anova(rdaGeoLoc)
#Plot your results
#plot(rdaGeoLoc)

#Condition out Locale so that you are only looking at both Geography and Climate
rdaGeoClim <- rda(PCOA_gigas$vectors~lat + long + bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19 + Condition(locale), data=gigas_clim_data_extracted)

#Summarize this run
#summary(rdaGeoClim)
#Get Rsquared ajusted values
RsquareAdj(rdaGeoClim)
#Conduct the rda anova
anova(rdaGeoClim)
#Plot your results
#plot(rdaGeoClim)

#Condition out Geography so that you are only looking at both Locale and Climate
rdaLocClim <- rda(PCOA_gigas$vectors~bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19 + locale + Condition(lat + long), data=gigas_clim_data_extracted)

#Summarize this run
#summary(rdaLocClim)
#Get Rsquared ajusted values
RsquareAdj(rdaLocClim)
#Conduct the rda anova
anova(rdaLocClim)
#Plot your results
#plot(rdaLocClim)

#To get the % varaince explained
var <- varpart(PCOA_gigas$vectors,~bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19, ~lat + long, ~locale, data=gigas_clim_data_extracted)

pdf("gigas_varpartsprop.pdf",paper="a4r")
plot(var)
dev.off()