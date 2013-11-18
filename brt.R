src/org/bccvl/compute/rscripts/brt.R# set CRAN mirror in case we need to download something
# TODO: this should be done on demand or on user basis...
r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# TODO: alse creating and populating add on package location is something that should not be done system wide

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("dismo","SDMTools", "gbm", "rgdal", "pROC", "R2HTML", "png") #list the libraries needed
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
    project.brt=FALSE
} else {
    future.climate.scenario = stack(enviro.data.future)
}

# source helper functions (err.null, getModelObject, checkModelLayers, saveModelProject)
source(paste(function.path, "/my.Helper.Functions.R", sep=""))

###run the models and store models
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
# learning.rate = 0.01, bag.fraction = 0.75, site.weights = rep(1, nrow(data)),
# var.monotone = rep(0, length(gbm.x)), n.folds = 10, prev.stratify = TRUE,
# family = "bernoulli", n.trees = 50, step.size = n.trees, max.trees = 10000,
# tolerance.method = "auto", tolerance = 0.001, keep.data = FALSE, plot.main = TRUE,
# plot.folds = FALSE, verbose = TRUE, silent = FALSE, keep.fold.models = FALSE,
# keep.fold.vector = FALSE, keep.fold.fit = FALSE, ...)
# data input data.frame
# gbm.x predictor variables
# gbm.y response variable
# offset = NULL
# fold.vector = NULL a fold vector to be read in for cross validation with offsets
# tree.complexity = 1 sets the complexity of individual trees
# learning.rate = 0.01 sets the weight applied to individual trees
# bag.fraction = 0.75 sets the proportion of observations used in selecting variables
# site.weights = rep(1, nrow(data)) allows varying weighting for sites
# var.monotone = rep(0, length(gbm.x)) restricts responses to individual predictors to monotone
# n.folds = 10 number of folds
# prev.stratify = TRUE prevalence stratify the folds - only for presence/absence data
# family = "bernoulli" family - bernoulli (=binomial), poisson, laplace or gaussian
# n.trees = 50 number of initial trees to fit
# step.size = n.trees numbers of trees to add at each cycle
# max.trees = 10000 max number of trees to fit before stopping
# tolerance.method = "auto" method to use in deciding to stop - "fixed" or "auto"
# tolerance = 0.001 tolerance value to use - if method == fixed is absolute,
# if auto is multiplier * total mean deviance
# keep.data = FALSE Logical. keep raw data in final model
# plot.main = TRUE Logical. plot hold-out deviance curve
# plot.folds = FALSE Logical. plot the individual folds as well
# verbose = TRUE Logical. control amount of screen reporting
# silent = FALSE Logical. to allow running with no output for simplifying model)
# keep.fold.models = FALSE Logical. keep the fold models from cross valiation
# keep.fold.vector = FALSE Logical. allows the vector defining fold membership to be kept
# keep.fold.fit = FALSE Logical. allows the predicted values for observations from cross-validation
# to be kept

if (model.brt) {
    outdir = paste(wd,'/output_brt',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
    brt.data = rbind(occur[,c("lon", "lat", enviro.data.names)],bkgd[,c("lon", "lat", enviro.data.names)]) #setup the data as needed
	brt.data$pa = c(rep(1,nrow(occur)),rep(0,nrow(bkgd))) #setup the data as needed
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
        save(brt,file=paste(outdir,"/model.object.RData",sep='')) #save out the model object
        # NOTE the order of arguments in the predict function for brt; this is because
        # the function is defined outside of the dismo package
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

if (project.brt) {
    brt.obj = getModelObject("brt") # get the model object
    if (!is.null(brt.obj)) {
		predictors = checkModelLayers(brt.obj)
        # NOTE the order of arguments in the predict function for brt; this is because
        # the function is defined outside of the dismo package
        brt.proj = predict(future.climate.scenario, brt.obj, n.trees=brt.obj$gbm.call$best.trees) # predict for given climate scenario
        saveModelProjection(brt.proj, "brt", "future") # save output
		rm(list=c("brt.obj", "brt.proj")) #clean up the memory
    } else {
        write(paste("FAIL!", species, "Cannot load brt.obj from", wd, "output_brt", sep=": "), stdout())
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

if (evaluate.brt) {
    brt.obj = getModelObject("brt") # get the model object
    if (!is.null(brt.obj)) {
        # NOTE the order of arguments in the predict function for brt; this is because
        # the function is defined outside of the dismo package
        brt.eval = evaluate(p=occur, a=bkgd, model=brt.obj, n.trees=brt.obj$gbm.call$best.trees) # evaluate model
        brt.fit = c(brt.eval@presence, brt.eval@absence)
        brt.obs = c(rep(1, length(brt.eval@presence)), rep(0, length(brt.eval@absence)))
       # get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
        brt.combined.eval = sapply(model.accuracy,
            function(x){
                return(my.Find.Optim.Stat(Stat = x,
                                          Fit = brt.fit,
                                          Obs = brt.obs))
            })
        saveModelEvaluation(brt.eval, brt.combined.eval, "brt") # save output
		
		# create response curves
		createMarginalResponseCurves(brt.obj, "brt")

		# calculate variable importance (like biomod2, using correlations between predictions)
		calculateVariableImpt(brt.obj, "brt", 3)
		
		# calculate variable importance (like maxent, using decrease in AUC)
		calculatePermutationVarImpt(brt.obj, brt.eval, "brt")
				
		# create HTML file with accuracy measures
		generateHTML(paste(wd, "/output_brt", sep='')) 
		
        rm(list=c("brt.obj", "brt.eval", "brt.combined.eval")) #clean up the memory
    } else {
        write(paste("FAIL!", species, "Cannot load brt.obj from", wd, "output_brt", sep=": "), stdout())
	}
} # end if brt