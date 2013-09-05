#script to run to develop distribution models

# read in the arguments listed at the command line in shell script
args=(commandArgs(TRUE))  
# check to see if arguments are passed
if(length(args)==0){
    print("No arguments supplied.")
    # leave all args as default values
} else {
	for(i in 1:length(args)) { 
		eval(parse(text=args[[i]])) 
	}
	# expecting wd to be able to locate arguments file
}

# load arguments file 
source(paste(wd, "02.sprint.init.args.txt", sep=""))

###check if libraries are installed, install if necessary and then load them
necessary=c("dismo","SDMTools", "gbm") #list the libraries needed
installed = necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >=1) install.packages(necessary[!installed], dep = T) #if library is not installed, install it
for (lib in necessary) library(lib,character.only=T)#load the libraries

###read in the necessary observation, background and environmental data
populate.data = FALSE #variable to define if there is a need to generate occur & background environmental info
if (file.exists(paste(wd, "occur.RData", sep="")) && file.exists(paste(wd, "bkgd.RData", sep=""))) {
	load(paste(wd, "occur.RData", sep="")); load(paste(wd, "bkgd.RData", sep="")); #if files already exist, load in the data
	if (!all(colnames(occur)==c('lon','lat',enviro.data.names))) { populate.data=TRUE } #not the right data, we need to repopulate it
} else { populate.data=TRUE } # data does not exist, we need to generate it
if (populate.data) {
	occur = read.csv(occur.data) #read in the observation data lon/lat
	if (!is.null(bkgd.data)) bkgd = read.csv(bkgd.data) #read in teh background position data lon.lat
	for (ii in 1:length(current.enviro.data)) { cat(ii,'of',length(current.enviro.data),'\n') #cycle through each of the environmental datasets and append the data
		tasc = read.asc(current.enviro.data[ii]) #read in the envirodata
		occur[,enviro.data.names[ii]] = extract.data(cbind(occur$lon,occur$lat),tasc) #extract envirodata for observations
		if (!is.null(bkgd.data)) bkgd[,enviro.data.names[ii]] = extract.data(cbind(bkgd$lon,bkgd$lat),tasc) #extract envirodata for background data
	}
	save(occur,file=paste(wd, "occur.RData", sep="")); #write out the raw data for analysis
	if (!is.null(bkgd.data)) save(bkgd,file=paste(wd, "bkgd.RData", sep="")) #write out the raw data for analysis
}

# combine predictors into RasterStacks
current.climate.scenario = stack(current.enviro.data)
future.climate.scenario = stack(future.enviro.data)

### User-defined functions

## Needed for tryCatch'ing:
err.null <- function (e) return(NULL)

# function to save projection output raster
saveModelProjection = function(out.model, model.name, projectiontime) {
	model.dir = paste(wd, "output_", model.name, "/", sep="")
	writeRaster(out.model, paste(model.dir, projectiontime, sep="/"), format="GTiff")
}

# function to get model object
getModelObject = function(model.name) {
	model.dir = paste(wd, "output_", model.name, "/", sep=""); setwd(model.dir);
	model.obj = tryCatch(get(load(file=paste(model.dir, "model.object.RData", sep=""))), error = err.null)	
}

# function to save evaluate output
saveModelEvaluation = function(out.model, out.biomod.model) {
	model.dir = paste(getwd(), "/", sep="")
	save(out.model, file=paste(model.dir, "dismo.eval.object.RData", sep=''))	# save the 'dismo::ModelEvalution' object
	
	# save all the model accuracy statistics provided in both dismo and biomod2
	rownames(out.biomod.model) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
	write.csv(t(round(out.biomod.model, digits=3)), file=paste(getwd(), "/combined.modelEvaluation.csv", sep=""))
	# EMG no guarantee these value are correct
	
	# save AUROC curve
	png(file=paste(model.dir, "AUC.png", sep='')); plot(out.model, 'ROC'); dev.off()
}

# source my modified version of biomod2's Evaluate.models.R for consistent model accuracy statistics
source("/home/jc140298/sprint2/my.Evaluate.models.R")

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
#	a Raster* object
# NOTE: env vars must be numerical

if (model.bioclim) {
	if (!all(enviro.data.type=="continuous")) {
		warning("bioclim not run because categorical data cannot be used")
	} else {
		outdir = paste(wd,'output_bioclim/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
		bc = tryCatch(bioclim(x=occur[,enviro.data.names]), error = err.null) #run bioclim with matrix of enviro data
		if (!is.null(bc)) {		
			save(bc,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
			save(bc,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
			bioclim.proj = predict(bc, current.climate.scenario, tails=opt.tails)	# predict for CURRENT climate scenario
			saveModelProjection(bioclim.proj, "bioclim", "current") # save output
		} else {
			write(paste("FAIL!", species, "Cannot create bioclim model object", sep=": "), stdout())
		}			
	} # end if continuous
} # end if



#############################################################################################
#
# MACHINE LEARNING METHODS - use both presence and absence or background data: Maxent, BRT
#
#############################################################################################

###############
#
# BRT
#
###############

# gbm.step(data, gbm.x, gbm.y, offset = NULL, fold.vector = NULL, tree.complexity = 1, 
#	learning.rate = 0.01, bag.fraction = 0.75, site.weights = rep(1, nrow(data)), 
#	var.monotone = rep(0, length(gbm.x)), n.folds = 10, prev.stratify = TRUE, 
#	family = "bernoulli", n.trees = 50, step.size = n.trees, max.trees = 10000,
#	tolerance.method = "auto", tolerance = 0.001, keep.data = FALSE, plot.main = TRUE, 
#	plot.folds = FALSE, verbose = TRUE, silent = FALSE, keep.fold.models = FALSE, 
#	keep.fold.vector = FALSE, keep.fold.fit = FALSE, ...)
# data input data.frame
# gbm.x predictor variables
# gbm.y	response variable
# offset = NULL
# fold.vector = NULL	a fold vector to be read in for cross validation with offsets
# tree.complexity = 1	sets the complexity of individual trees
# learning.rate = 0.01	sets the weight applied to individual trees
# bag.fraction = 0.75	sets the proportion of observations used in selecting variables
# site.weights = rep(1, nrow(data))	allows varying weighting for sites
# var.monotone = rep(0, length(gbm.x))	restricts responses to individual predictors to monotone
# n.folds = 10	number of folds
# prev.stratify = TRUE	prevalence stratify the folds - only for presence/absence data
# family = "bernoulli"	family - bernoulli (=binomial), poisson, laplace or gaussian
# n.trees = 50	number of initial trees to fit
# step.size = n.trees	numbers of trees to add at each cycle
# max.trees = 10000	max number of trees to fit before stopping
# tolerance.method = "auto"	method to use in deciding to stop - "fixed" or "auto"
# tolerance = 0.001	tolerance value to use - if method == fixed is absolute, 
#	if auto is multiplier * total mean deviance
# keep.data = FALSE	Logical. keep raw data in final model
# plot.main = TRUE	Logical. plot hold-out deviance curve
# plot.folds = FALSE	Logical. plot the individual folds as well
# verbose = TRUE	Logical. control amount of screen reporting
# silent = FALSE	Logical. to allow running with no output for simplifying model)
# keep.fold.models = FALSE 	Logical. keep the fold models from cross valiation
# keep.fold.vector = FALSE	Logical. allows the vector defining fold membership to be kept
# keep.fold.fit = FALSE	Logical. allows the predicted values for observations from cross-validation 
#	to be kept

if (model.brt) {
	outdir = paste(wd,'output_brt/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	brt.data = rbind(occur,bkgd); brt.data$pa = c(rep(1,nrow(occur)),rep(0,nrow(bkgd))) #setup the data as needed
	brt = tryCatch(gbm.step(data=brt.data, gbm.x=which(names(brt.data) %in% enviro.data.names), 
		gbm.y=which(names(brt.data)=='pa'), 
		fold.vector = brt.fold.vector, 
		tree.complexity = brt.tree.complexity, 
		learning.rate = brt.learning.rate, 
		bag.fraction = brt.bag.fraction, 
		#site.weights = brt.site.weights, 
		#var.monotone = brt.var.monotone, 
		n.folds = brt.n.folds, 
		prev.stratify = brt.prev.stratify, 
		family = brt.family, 
		n.trees = brt.n.trees, 
		step.size = brt.step.size, 
		max.trees = brt.max.trees, 
		tolerance.method = brt.tolerance.method, 
		tolerance = brt.tolerance, 
		keep.data = brt.keep.data, 
		plot.main = brt.plot.main, 
		plot.folds = brt.plot.folds, 
		verbose = brt.verbose, 
		silent = brt.silent, 
		keep.fold.models = brt.keep.fold.models, 
		keep.fold.vector = brt.keep.fold.vector, 
		keep.fold.fit = brt.keep.fold.fit), error = err.null) #run the algorithm
	if (!is.null(brt)) {	
		save(brt,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
		save(brt,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
		# NOTE the order of arguments in the  predict function for brt; this is because
		#	the function is defined outside of the dismo package
		brt.proj = predict(current.climate.scenario, brt, n.trees=brt$gbm.call$best.trees) # predict for CURRENT climate scenario
		saveModelProjection(brt.proj, "brt", "current") 
	} else {
		write(paste("FAIL!", species, "Cannot create brt model object", sep=": "), stdout())
	}
}

###############
#
# predict(object, x, ext=NULL, filename="", progress='text', ...)
#
# object A fitted model of class Bioclim, Domain, MaxEnt, ConvexHull, or Mahalanobis (classes that inherit from DistModel)
# x  A Raster* object or a data.frame 
# ext  An extent object to limit the prediction to a sub-region of 'x'. Or an object that can be coerced to an Extent object by extent; such as a Raster* or Spatial* object  
# filename  Output filename for a new raster; if NA the result is not written to a file but returned with the RasterLayer object, in the data slot  
# progress  Character. Valid values are "" (no progress bar), "text" and "windows" (on that platform only)  
# ...  Additional model specific arguments. And additional arguments for file writing as for writeRaster  
#
# For maxent models, there is an additional argument 'args' used to pass arguments (options) to the maxent software.
# For bioclim models, there is an additional argument 'tails' which you can use to ignore the left or right tail of the percentile distribution for a variable.
# For geoDist models, there is an additional argument fun that allows you to use your own (inverse) distance function, and argument scale=1 that allows you to scale 
#	the values (distances smaller than this value become one, and the others are divided by this value before computing the inverse distance).
# For spatial predictions with BRT, randomForest, etc., see 'predict' in the Raster package
#
###############
 

###project the models onto FUTURE climate and save raster files
if (project.bioclim) {
	bioclim.obj = getModelObject("bioclim")	# get the model object
	if (!is.null(bioclim.obj)) {
		bioclim.proj = predict(bioclim.obj, future.climate.scenario, tails=opt.tails)	# predict for given climate scenario
		saveModelProjection(bioclim.proj, "bioclim", "future") # save output
	} else {
		write(paste("FAIL!", species, "Cannot load bioclim.obj from", wd, "output_bioclim", sep=": "), stdout())
	}
} # end if bioclim


if (project.brt) {
	brt.obj = getModelObject("brt") # get the model object
	if (!is.null(brt.obj)) {
		# NOTE the order of arguments in the  predict function for brt; this is because
		#	the function is defined outside of the dismo package
		brt.proj = predict(future.climate.scenario, brt.obj, n.trees=brt.obj$gbm.call$best.trees) # predict for given climate scenario
		saveModelProjection(brt.proj, "brt", "future") 	# save output
	} else {
		write(paste("FAIL!", species, "Cannot load brt.obj from", wd, "output_brt", sep=": "), stdout())
	}
}

###############
#
# evaluate(p, a, model, x, tr, ...)
#
# p presence points (x and y coordinate or SpatialPoints* object)
#	Or, if x is missing, values at presence points (EMG: values returned by a predict())
#	Or, a matrix with values to compute predictions for
# a absence points (x and y coordinate or SpatialPoints* object)
#	Or, if x is missing, values at absence points (EMG: values returned by a predict())
#	Or, a matrix with values to compute predictions for
# model any fitted model, including objects inheriting from 'DistModel'; not used when x is missing
# x 	Optional. Predictor values (object of class Raster*). If present, p and a are interpreted
#	as (spatial) points (EMG: lon/lat)
# tr 	Optional. a vector of threshold values to use for computing the confusion matrices
# ...	Additional arguments for the predict function (EMG: evaluate() calls predict())
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
		# need the predictions and observed values to create confusion matrices for accuracy statistics
		bioclim.fit = c(bioclim.eval@presence, bioclim.eval@absence)
		bioclim.obs = c(rep(1, nrow(occur)), rep(0, nrow(bkgd)))
		# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
		bioclim.combined.eval = sapply(model.accuracy,
								function(x){
                                return(my.Find.Optim.Stat(Stat = x,
														  Fit = bioclim.fit,
                                                          Obs = bioclim.obs))
                                })
		saveModelEvaluation(bioclim.eval, bioclim.combined.eval)	# save output
		rm(list=c("bioclim.obj", "bioclim.eval", "bioclim.combined.eval")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load bioclim.obj from", wd, "output_bioclim", sep=": "), stdout())
	}
} # end if bioclim

if (evaluate.brt) {
	brt.obj = getModelObject("brt") # get the model object
	if (!is.null(brt.obj)) {
		# NOTE the order of arguments in the  predict function for brt; this is because
		#	the function is defined outside of the dismo package
		brt.eval = evaluate(p=occur, a=bkgd, model=brt.obj, n.trees=brt.obj$gbm.call$best.trees) # evaluate model
		brt.fit = c(brt.eval@presence, brt.eval@absence)
		brt.obs = c(rep(1, nrow(occur)), rep(0, nrow(bkgd)))
		# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
		brt.combined.eval = sapply(model.accuracy,
								function(x){
                                return(my.Find.Optim.Stat(Stat = x,
														  Fit = brt.fit,
                                                          Obs = brt.obs))
                                })
		saveModelEvaluation(brt.eval, brt.combined.eval) 	# save output
		rm(list=c("brt.obj", "brt.eval", "brt.combined.eval")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load brt.obj from", wd, "output_brt", sep=": "), stdout())
	}
}