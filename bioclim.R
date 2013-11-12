src/org/bccvl/compute/rscripts/bioclim.R# set CRAN mirror in case we need to download something
# TODO: this should be done on demand or on user basis...
r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# TODO: alse creating and populating add on package location is something that should not be done system wide

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("dismo","SDMTools", "rgdal", "pROC", "R2HTML", "png") #list the libraries needed
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
    project.bioclim=FALSE
} else {
    future.climate.scenario = stack(enviro.data.future)
}

# source helper functions (err.null, getModelObject, checkModelLayers, saveModelProject)
source(paste(function.path, "/my.Helper.Functions.R", sep=""))

###run the models and store models
#################################################################################
#
# PROFILE METHODS - only consider presence points: Bioclim, Domain, and Mahal
#
#################################################################################

###############
#
# BIOCLIM
#
###############

# bioclim(x, p, ...)
# x is a Raster* object or matrix
# p is a two column matrix or SpatialPoints* object
# if p is missing, x is a matrix of values of env vars at known locations of occurrence
# if p is present, it is the location of occurrence and used to extract values for env vars from x,
#       a Raster* object
# NOTE: env vars must be numerical

if (model.bioclim) {
    if (!all(enviro.data.type=="continuous")) {
        warning("bioclim not run because categorical data cannot be used")
    } else {
        outdir = paste(wd,'/output_bioclim',sep='')
        dir.create(outdir,recursive=TRUE) #create the output directory
        bc = tryCatch(bioclim(x=occur[,enviro.data.names]), error = err.null) #run bioclim with matrix of enviro data
        if (!is.null(bc)) {
            save(bc, file=paste(outdir,"/model.object.RData",sep='')) #save out the model object
            bc.proj = predict(bc, current.climate.scenario, tails=opt.tails, ext=opt.ext)   # predict for given climate scenario
            saveModelProjection(bc.proj, "bioclim", "current") # save output
			rm(list=c("bc", "bc.proj")); #clean up memory
        } else {
            write(paste("FAIL!", species, "Cannot create bioclim model object", sep=": "), stdout())
        }
    } # end if continuous
} # end if


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
if (project.bioclim) {
    bioclim.obj = getModelObject("bioclim")	# get the model object
    if (!is.null(bioclim.obj)) {
		predictors = checkModelLayers(bioclim.obj)
        bioclim.proj = predict(bioclim.obj, predictors, tails=opt.tails, ext=opt.ext)	# predict for given climate scenario
        saveModelProjection(bioclim.proj, "bioclim", "future") # save output
		rm(list=c("bioclim.obj", "bioclim.proj")) #clean up the memory
    } else {
        write(paste("FAIL!", species, "Cannot load bioclim.obj from", wd, "output_bioclim", sep=": "), stdout())
    }
} # end if bioclim


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
if (evaluate.bioclim) {
	bioclim.obj = getModelObject("bioclim")	# get the model object
	if (!is.null(bioclim.obj)) {
		bioclim.eval = evaluate(p=occur, a=bkgd, model=bioclim.obj)	# evaluate model using dismo's evaluate
		
		# need predictions and observed values to create confusion matrices for accuracy statistics
		bioclim.fit = c(bioclim.eval@presence, bioclim.eval@absence)
		bioclim.obs = c(rep(1, length(bioclim.eval@presence)), rep(0, length(bioclim.eval@absence)))

		# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
		bioclim.combined.eval = sapply(model.accuracy, function(x){
			return(my.Find.Optim.Stat(Stat = x, Fit = bioclim.fit, Obs = bioclim.obs))
		})
		saveModelEvaluation(bioclim.eval, bioclim.combined.eval, "bioclim")	# save output
		
		# create response curves
		createMarginalResponseCurves(bioclim.obj, "bioclim")

		# calculate variable importance (like biomod2, using correlations between predictions)
		calculateVariableImpt(bioclim.obj, "bioclim", 3)
		
		# calculate variable importance (like maxent, using decrease in AUC)
		calculatePermutationVarImpt(bioclim.obj, bioclim.eval, "bioclim")
		
		# create HTML file with accuracy measures
		generateHTML(species, paste(wd, "/output_bioclim", sep='')) 
		
		rm(list=c("bioclim.obj", "bioclim.eval", "bioclim.combined.eval")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load bioclim.obj from", wd, "output_bioclim", sep=": "), stdout())
	}
} # end if bioclim