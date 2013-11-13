src/org/bccvl/compute/rscripts/domain.R# set CRAN mirror in case we need to download something
# TODO: this should be done on demand or on user basis...
r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# TODO: alse creating and populating add on package location is something that should not be done system wide

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("dismo","SDMTools", "gstat", "rgdal", "pROC", "R2HTML", "png") #list the libraries needed
installed = necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >=1) {
    install.packages(necessary[!installed], dep = T) #if library is not installed, install it
}
for (lib in necessary) {
    library(lib,character.only=T) #load the libraries
}

###read in the necessary observation, background and environmental data
#setwd(wd) #set the working directory
populate.data = FALSE #variable to define if there is a need to generate occur & background environmental info
if (file.exists(paste(wd, "/occur.RData", sep=""))
    && file.exists(paste(wd, "/bkgd.RData", sep=""))) {
    load(paste(wd, "/occur.RData", sep=""))
    load(paste(wd, "/bkgd.RData", sep="")) #if files already exist, load in the data
    if (!all(colnames(occur)==c('lon','lat',enviro.data.names))) {
        populate.data=TRUE #not the right data, we need to repopulate it
    }
} else {
    populate.data=TRUE # data does not exist, we need to generate it
}
if (populate.data) {
    occur = read.csv(occur.data) #read in the observation data lon/lat
    if (!is.null(bkgd.data)) {
        bkgd = read.csv(bkgd.data) #read in teh background position data lon.lat
    }
    for (ii in 1:length(enviro.data.current)) {
        cat(ii,'of',length(enviro.data.current),'\n') #cycle through each of the environmental datasets and append the data
        #tasc = read.asc(enviro.data.current[ii]) #read in the envirodata
        tasc = readGDAL(enviro.data.current[ii]) #read in the envirodata
        occur[,enviro.data.names[ii]] = extract.data(cbind(occur$lon,occur$lat),tasc) #extract envirodata for observations
        if (!is.null(bkgd.data)) bkgd[,enviro.data.names[ii]] = extract.data(cbind(bkgd$lon,bkgd$lat),tasc) #extract envirodata for background data
    }
    save(occur,file=paste(wd, "/occur.RData", sep="")) #write out the raw data for analysis
    if (!is.null(bkgd.data)) {
        save(bkgd,file=paste(wd, "/bkgd.RData", sep="")) #write out the raw data for analysis
    }
}

current.climate.scenario = stack(enviro.data.current)
if (is.null(enviro.data.future)) {
    project.domain=FALSE
} else {
    future.climate.scenario = stack(enviro.data.future)
}

# source helper functions (err.null, getModelObject, checkModelLayers, saveModelProject)
source(paste(function.path, "/my.Helper.Functions.R", sep=""))

###run the models and store models
#############################################################################################
#
# GEOGRAPHIC MODELS - use the geographic location of known occurrences
#
#############################################################################################

############### Spatial-only models for presence/background (or absence) data ###############

###############
#
# GEOIDW - inverse distance weighted interpolation
#
###############

# geoIDW(p, a, ...)
# p presence points; two column matrix, data.frame or SpatialPoints* object
# a absence points; must be of the same class as 'p'
# ... none implemented

if (model.geoIDW) {
	outdir = paste(wd,'/output_geoIDW',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	gidw = tryCatch(geoIDW(p=occur[,c('lon','lat')], a=bkgd[,c('lon','lat')]), error = err.null) #run the algorithm
	if (!is.null(gidw)) {	
		save(gidw,file=paste(outdir,"/model.object.RData",sep='')) #save out the model object
		bc.proj = predict(bc, current.climate.scenario, ext=opt.ext)   # predict for given climate scenario
		saveModelProjection(bc.proj, "bioclim", "current") # save output
		rm(list=c("bc", "bc.proj")); #clean up memory
	} else {
		write(paste("FAIL!", species, "Cannot create geoIDW model object", sep=": "), stdout())
	}
}


###############
#
# predict(object, x, ext=NULL, filename="", progress='text', ...)
#
# object A fitted model of class Bioclim, Domain, MaxEnt, ConvexHull, or Mahalanobis (classes that inherit from DistModel)
# x A Raster* object or a data.frame
# ext An extent object to limit the prediction to a sub-region of 'x'. Or an object that can be coerced to an Extent object by extent; such as a Raster* or Spatial* object
# filename Output filename for a new raster; if NA the result is not written to a file but returned with the RasterLayer object, in the data slot
# progress Character. Valid values are "" (no progress bar), "text" and "windows" (on that platform only)
# ... Additional model specific arguments. And additional arguments for file writing as for writeRaster
#
# For maxent models, there is an additional argument 'args' used to pass arguments (options) to the maxent software.
# For bioclim models, there is an additional argument 'tails' which you can use to ignore the left or right tail of the percentile distribution for a variable.
# For geoDist models, there is an additional argument fun that allows you to use your own (inverse) distance function, and argument scale=1 that allows you to scale
# the values (distances smaller than this value become one, and the others are divided by this value before computing the inverse distance).
# For spatial predictions with BRT, randomForest, etc., see 'predict' in the Raster package
#
###############

###project the models onto FUTURE climate and save raster files
if (project.geoIDW) {
	geoIDW.obj = getModelObject("geoIDW") # get the model object
	if (!is.null(geoIDW.obj)) {
		predictors = checkModelLayers(geoIDW.obj)
		geoIDW.proj = predict(geoIDW.obj, predictors, ext=opt.ext) # predict for given climate scenario
		saveModelProjection(geoIDW.proj, "geoIDW", "future") 	# save output
		rm(list=c("geoIDW.obj", "geoIDW.proj")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load geoIDW.obj from", wd, "output_geoIDW", sep=": "), stdout())
	}
}


###############
#
# evaluate(p, a, model, x, tr, ...)
#
# p presence points (x and y coordinate or SpatialPoints* object)
# Or, if x is missing, values at presence points (EMG: values returned by a predict())
# Or, a matrix with values to compute predictions for
# a absence points (x and y coordinate or SpatialPoints* object)
# Or, if x is missing, values at absence points (EMG: values returned by a predict())
# Or, a matrix with values to compute predictions for
# model any fitted model, including objects inheriting from 'DistModel'; not used when x is missing
# x Optional. Predictor values (object of class Raster*). If present, p and a are interpreted
# as (spatial) points (EMG: lon/lat)
# tr Optional. a vector of threshold values to use for computing the confusion matrices
# ... Additional arguments for the predict function (EMG: evaluate() calls predict())
#
# 'ModelEvaluation' output based on Fielding and Bell (1997) with attributes:
# presence - presence data used
# absence - absence data used
# np - number of presence points
# na - number of absence points
# auc - Area under the receiver operator (ROC) curve
# pauc - p-value for the AUC (for the Wilcoxon test W statistic
# cor - Correlation coefficient
# pcor - p-value for correlation coefficient
# t - vector of thresholds used to compute confusion matrices
# confusion - confusion matrices
# prevalence - Prevalence
# ODP - Overall diagnostic power
# CCR - Correct classification rate
# TPR - True positive rate
# TNR - True negative rate
# FPR - False positive rate
# FNR - False negative rate
# PPP - Positive predictive power
# NPP - Negative predictive power
# MCR - Misclassification rate
# OR - Odds-ratio
# kappa - Cohen's kappa
#
###############

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy = c(dismo.eval.method, biomod.models.eval.meth)

###evaluate the models and save the outputs
if (evaluate.geoIDW) {
	geoIDW.obj = getModelObject("geoIDW") # get the model object
	if (!is.null(geoIDW.obj)) {
		geoIDW.eval = evaluate(model=geoIDW.obj, p=occur[c("lon","lat")], 
			a=bkgd[c("lon","lat")]) # evaluate model using dismo's evaluate
		
		# need predictions and observed values to create confusion matrices for accuracy statistics
		geoIDW.fit = c(geoIDW.eval@presence, geoIDW.eval@absence)
		geoIDW.obs = c(rep(1, length(geoIDW.eval@presence)), rep(0, length(geoIDW.eval@absence)))

		# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
		geoIDW.combined.eval = sapply(model.accuracy, function(x){
			return(my.Find.Optim.Stat(Stat = x, Fit = geoIDW.fit, Obs = geoIDW.obs))
		})
		saveModelEvaluation(geoIDW.eval, geoIDW.combined.eval, "geoIDW")	# save output
						
		# create response curves
		createMarginalResponseCurves(geoIDW.obj, "geoIDW")

		# calculate variable importance (like biomod2, using correlations between predictions)
		calculateVariableImpt(geoIDW.obj, "geoIDW", 3)
		
		# calculate variable importance (like maxent, using decrease in AUC)
		calculatePermutationVarImpt(geoIDW.obj, geoIDW.eval, "geoIDW")
		
		# create HTML file with accuracy measures
		generateHTML(species, paste(wd, "/output_geoIDW", sep='')) 
		
		rm(list=c("geoIDW.obj", "geoIDW.eval", "geoIDW.combined.eval")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load geoIDW.obj from", wd, "output_geoIDW", sep=": "), stdout())
	}
}