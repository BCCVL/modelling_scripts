#script to create biodiversity measures for BCCVL

library(raster)
library(SDMTools)

# define working directory
wd = "/home/jc140298/ibrahim"

# get a list of species with projection maps
species.list = c("ABT", "ANTADUS", "ANTFLAV")

#########################################	
#	
# species richness
#
#########################################

# set the threshold
thresholds = c(0.215, 0.056, 0.123)
# EMG These are the KAPPA cutoff's reported by the modelEvaluation.csv, other 
#	thresholds should be available for selection

# get the .tif file as raster for the first species to start
first.raster = raster(paste(wd, "/", species.list[1], "/output_bioclim/future.tif", sep=""))
# convert raster to asc to apply threshold
richness.asc = asc.from.raster(first.raster)
# apply threshold
richness.asc[which(richness.asc<thresholds[1])]=0
richness.asc[which(richness.asc>=thresholds[1])]=1

# add the rest of the species
for (sp in species.list[-1]) {

	# create an index to access threshold values
	index = which(species.list == sp)
	
	# get the .tif file as raster
	sp.raster = raster(paste(wd, "/", sp, "/output_bioclim/future.tif", sep=""))
	# convert raster to asc to apply threshold
	sp.asc = asc.from.raster(sp.raster)
	# apply the threshold
	sp.asc[which(sp.asc<thresholds[index])]=0
	sp.asc[which(sp.asc>=thresholds[index])]=1
	
	# add them together
	richness.asc = richness.asc + sp.asc
}

# convert richness asc back to raster
richness.raster = raster.from.asc(richness.asc)
# save richness raster as tif
writeRaster(richness.raster, paste(wd, "richness", sep="/"), format="GTiff", 
	options="COMPRESS=LZW", overwrite=TRUE)

#########################################	
#	
# calculate the gain/loss in species
#
#########################################

# get the baseline .tif as raster
current.raster = raster(paste(wd, "/", species.list[1], "/output_bioclim/current.tif", sep=""))
# convert raster to asc to apply threshold
current.asc = asc.from.raster(current.raster)
# apply the threshold
current.asc[which(current.asc<thresholds[1])]=0
current.asc[which(current.asc>=thresholds[1])]=1

# get the .tif as raster to be used for comparison
future.raster = raster(paste(wd, "/", species.list[1], "/output_bioclim/future.tif", sep=""))
# convert raster to asc to apply threshold
future.asc = asc.from.raster(future.raster)
# apply the threshold
future.asc[which(future.asc<thresholds[1])]=0
future.asc[which(future.asc>=thresholds[1])]=1
# EMG Use the same thresholdf for both?

# calculate the change in species number
change.asc = current.asc - future.asc
# convert asc to raster
change.raster = raster.from.asc(change.asc)

# save the change raster as tif
writeRaster(change.raster, paste(wd, "change", sep="/"), format="GTiff", 
	options="COMPRESS=LZW", overwrite=TRUE)