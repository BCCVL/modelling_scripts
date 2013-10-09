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
	# expecting wd and species to be able to locate arguments file
}

# load arguments file
load(paste(wd, "01.init.args.model.current.", species, ".RData", sep=""))

###check if libraries are installed, install if necessary and then load them
necessary=c("dismo","SDMTools","gbm","gstat","deldir", "biomod2") #list the libraries needed
installed = necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >=1) install.packages(necessary[!installed], dep = T) #if library is not installed, install it
for (lib in necessary) library(lib,character.only=T)#load the libraries

###read in the necessary observation, background and environmental data
#setwd(wd) #set the working directory
populate.data = FALSE #variable to define if there is a need to generate occur & background environmental info
if (file.exists(paste(wd, "occur.RData", sep="")) && file.exists(paste(wd, "bkgd.RData", sep=""))) {
	load(paste(wd, "occur.RData", sep="")); load(paste(wd, "bkgd.RData", sep="")); #if files already exist, load in the data
	if (!all(colnames(occur)==c('lon','lat',enviro.data.names))) { populate.data=TRUE } #not the right data, we need to repopulate it
} else { populate.data=TRUE } # data does not exist, we need to generate it
if (populate.data) {
	occur = read.csv(occur.data) #read in the observation data lon/lat
	bkgd = read.csv(bkgd.data) #read in teh background position data lon.lat
	for (ii in 1:length(enviro.data)) { cat(ii,'of',length(enviro.data),'\n') #cycle through each of the environmental datasets and append the data
		tasc = read.asc(enviro.data[ii]) #read in the envirodata
		occur[,enviro.data.names[ii]] = extract.data(cbind(occur$lon,occur$lat),tasc) #extract envirodata for observations
		bkgd[,enviro.data.names[ii]] = extract.data(cbind(bkgd$lon,bkgd$lat),tasc) #extract envirodata for background data
	}
	save(occur,file=paste(wd, "occur.RData", sep="")); save(bkgd,file=paste(wd, "bkgd.RData", sep="")) #write out the raw data for analysis
}

## Needed for tryCatch'ing:
err.null <- function (e) return(NULL)

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
			rm(bc); #clean up memory
		} else {
			write(paste("FAIL!", species, "Cannot create bioclim model object", sep=": "), stdout())
		}			
	} # end if continuous
} # end if

###############
#
# DOMAIN
#
###############

# domain(x, p, ...)
# x is a Raster* object or matrix
# p is a two column matrix or SpatialPoints* object
# if p is missing, x is a matrix of values of env vars at known locations of occurrence
# if p is present, it is the location of occurrence and used to extract values for env vars from x,
#	a Raster* object
# NOTE: env vars must be numerical

if (model.domain) {
	if (!all(enviro.data.type=="continuous")) {
		warning("domain not run because categorical data cannot be used")
	} else {
		outdir = paste(wd,'output_domain/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
		dm = tryCatch(domain(x=occur[,enviro.data.names]), error = err.null) #run domain with matrix of enviro data
		if (!is.null(dm)) {	
			save(dm,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
			save(dm,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
			rm(dm); #clean up memory
		} else {
			write(paste("FAIL!", species, "Cannot create domain model object", sep=": "), stdout())
		}
	}
}

###############
#
# MAHALANOBIS
#
###############

# mahal(x, p, ...)
# x is a Raster* object or matrix
# p is a two column matrix or SpatialPoints* object
# if p is missing, x is a matrix of values of env vars at known locations of occurrence
# if p is present, it is the location of occurrence and used to extract values for env vars from x,
#	a Raster* object
# NOTE: env vars must be numerical

if (model.mahal) {
	if (!all(enviro.data.type=="continuous")) {
		warning("Mahal not run because categorical data cannot be used")
	} else {
		outdir = paste(wd,'output_mahal/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
		mm = tryCatch(mahal(x=occur[,enviro.data.names]), error = err.null) #run mahal with matrix of enviro data
		if (!is.null(mm)) {	
			save(mm,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
			save(mm,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
			rm(mm); #clean up memory
		} else {
			write(paste("FAIL!", species, "Cannot create mahal model object", sep=": "), stdout())
		}
	}
}

#############################################################################################
#
# GEOGRAPHIC MODELS - use the geographic location of known occurrences
#
#############################################################################################

############### Spatial-only models for presence data ###############

###############
#
# GEOGRAPHIC DISTANCE
#
###############

# geoDist(p, ...)
# p point locations (presence); two column matrix, data.frame or SpatialPoints* object
# ... you must supply a lonlat= argument(logical), unless p is a SpatialPoints* object and has a
#	valid CRS
# ... you can also supply an additional argument 'a' for absence points (currently ignored.); 
#	argument 'a' should be of the same class as argument 'p'

if (model.geodist) {
	outdir = paste(wd,'output_geodist/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	gd = tryCatch(geoDist(p=occur[,c('lon','lat')], lonlat=TRUE), error = err.null) #run geodist 
	if (!is.null(gd)) {	
		save(gd,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
		save(gd,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
		rm(gd); #clean up memory
	} else {
		write(paste("FAIL!", species, "Cannot create geodist model object", sep=": "), stdout())
	}
}

###############
#
# CONVEX HULLS
#
###############

# convHull(p, ...)
# p point locations (presence), two column matrix, data.frame or SpatialPoints* object
# ... you can supply an argument n (>=1) to get n convex hulls around subset of the points
# ... you can also set n=1:x, to get a set of overlapping polygons consisting of 1 to x parts; i.e.
#	the first polygon has 1 part, the second has 2 parts and x has x parts

if (model.convHull) {
	outdir = paste(wd,'output_convHull/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	ch = tryCatch(convHull(p=occur[,c('lon','lat')]), error = err.null) #run convex hull 
	if (!is.null(ch)) {		
		save(ch,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
		save(ch,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
		rm(ch); #clean up memory
	} else {
		write(paste("FAIL!", species, "Cannot create convHull model object", sep=": "), stdout())
	}
}

###############
#
# CIRCLES
#
###############

# circles(p, ...)
# p point locations (presence), two column matrix, data.frame or SpatialPoints* object
# d the radius of each circle in meters; a single number or a vector with elements corresponding to
#	rows in 'p'; if missing the diameter is computed from the inter-point distance
# n how many vertices in the circle? default is 360
# lonlat are these longitude/latitude data? default value is false
# r radius of the earth; only relevant for longitude/latitude data; default is 6378137 m

if (model.circles) {
	outdir = paste(wd,'output_circles/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	cc = tryCatch(circles(p=occur[,c('lon','lat')], lonlat=TRUE), error = err.null) #run circles 
	if (!is.null(cc)) {	
		save(cc,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
		save(cc,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
		rm(cc); #clean up memory
	} else {
		write(paste("FAIL!", species, "Cannot create circles model object", sep=": "), stdout())
	}
}

############### Spatial-only models for presence/background (or absence) data ###############

###############
#
# GEOIDW - inverse distance weighted interpolation
#
###############

# geoIDS(p, a, ...)
# p presence points; two column matrix, data.frame or SpatialPoints* object
# a absence points; must be of the same class as 'p'
# ... none implemented

if (model.geoIDW) {
	outdir = paste(wd,'output_geoIDW/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	gidw = tryCatch(geoIDW(p=occur[,c('lon','lat')], a=bkgd[,c('lon','lat')]), error = err.null) #run the algorithm
	if (!is.null(gidw)) {	
		save(gidw,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
		save(gidw,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
		rm(gidw); #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot create geoIDW model object", sep=": "), stdout())
	}
}

###############
#
# VORONOIHULL
#
###############

# voronoiHull(p, a, ...)
# p presence points; two column matrix, data.frame or SpatialPoints* object
# a absence points; must be of the same class as 'p'

if (model.voronoiHull) {
	outdir = paste(wd,'output_voronoiHull/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	vh = tryCatch(voronoiHull(p=occur[,c('lon','lat')], a=bkgd[,c('lon','lat')]), error = err.null) #run the algorithm
	if (!is.null(vh)) {	
		save(vh,file=paste(outdir,"model.object.RData",sep='')) #save out the model object
		save(vh,file=paste(outdir,"model.object.Rascii",sep=''), ascii=TRUE) #save out the model object as ascii for Daniel
		rm(vh); #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot create voronoiHull model object", sep=": "), stdout())
	}
}

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
		rm(brt); #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot create brt model object", sep=": "), stdout())
	}
}

###############
#
# MAXENT
#
###############

# maxent is being run as a system call with the same data. Arguments include:
# table below is "Flag -- Type -- Default -- Description"
# responsecurves -- boolean -- FALSE -- Create graphs showing how predicted relative probability of occurrence depends on the value of each environmental variable
# pictures -- boolean -- TRUE -- Create a .png image for each output grid
# jackknife -- boolean -- FALSE -- Measure importance of each environmental variable by training with each environmental variable first omitted, then used in isolation
# outputformat -- string -- logistic -- Representation of probabilities used in writing output grids. See Help for details
# outputfiletype -- string -- asc -- File format used for writing output grids
# outputdirectory -- directory --  -- Directory where outputs will be written. This should be different from the environmental layers directory.
# projectionlayers -- file/directory --  -- Location of an alternate set of environmental variables. Maxent models will be projected onto these variables.
# --  --  -- Can be a .csv file (in SWD format) or a directory containing one file per variable.
# --  --  -- Multiple projection files/directories can be separated by commas.
# samplesfile -- file --  -- Please enter the name of a file containing presence locations for one or more species.
# environmentallayers -- file/directory --  -- Environmental variables can be in a directory containing one file per variable,
# --  --  -- or all together in a .csv file in SWD format. Please enter a directory name or file name.
# randomseed -- boolean -- FALSE -- If selected, a different random seed will be used for each run, so a different random test/train partition
# --  --  -- will be made and a different random subset of the background will be used, if applicable.
# logscale -- boolean -- TRUE -- If selected, all pictures of models will use a logarithmic scale for color-coding.
# warnings -- boolean -- TRUE -- Pop up windows to warn about potential problems with input data.
# --  --  -- Regardless of this setting, warnings are always printed to the log file.
# tooltips -- boolean -- TRUE -- Show messages that explain various parts of the interface, like this message
# askoverwrite -- boolean -- TRUE -- If output files already exist for a species being modeled,
# --  --  -- pop up a window asking whether to overwrite or skip. Default is to overwrite.
# skipifexists -- boolean -- FALSE -- If output files already exist for a species being modeled,
# --  --  -- skip the species without remaking the model.
# removeduplicates -- boolean -- TRUE -- Remove duplicate presence records.
# --  --  -- If environmental data are in grids, duplicates are records in the same grid cell.
# --  --  -- Otherwise, duplicates are records with identical coordinates.
# writeclampgrid -- boolean -- TRUE -- Write a grid that shows the spatial distribution of clamping.
# --  --  -- At each point, the value is the absolute difference between prediction values with and without clamping.
# writemess -- boolean -- TRUE -- A multidimensional environmental similarity surface (MESS) shows where novel climate conditions exist in the projection layers.
# --  --  -- The analysis shows both the degree of novelness and the variable that is most out of range at each point.
# randomtestpoints -- integer -- 0 -- Percentage of presence localities to be randomly set aside as test points, used to compute AUC, omission etc.
# betamultiplier -- double -- 1 -- Multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.
# maximumbackground -- integer -- 10000 -- If the number of background points / grid cells is larger than this number, then this number of cells is chosen randomly for background points
# biasfile -- file --  -- Sampling is assumed to be biased according to the sampling distribution given in this grid file.
# --  --  -- Values in this file must not be zero or negative. MaxEnt will factor out the bias.
# --  --  -- Requires environmental data to be in grids, rather than a SWD format file
# testsamplesfile -- file --  -- Use the presence localities in this file to compute statistics (AUC, omission etc.)
# --  --  -- The file can contain different localities for different species.
# --  --  -- It takes precedence over the random test percentage.
# replicates -- integer -- 1 -- Number of replicate runs to do when cross-validating, bootstrapping or doing sampling with replacement runs
# replicatetype -- string -- crossvalidate -- If replicates > 1, do multiple runs of this type:
# --  --  -- Crossvalidate: samples divided into replicates folds; each fold in turn used for test data.
# --  --  -- Bootstrap: replicate sample sets chosen by sampling with replacement.
# --  --  -- Subsample: replicate sample sets chosen by removing random test percentage without replacement to be used for evaluation.
# perspeciesresults -- boolean -- FALSE -- Write separate maxentResults file for each species
# writebackgroundpredictions -- boolean -- FALSE -- Write .csv file with predictions at background points
# responsecurvesexponent -- boolean -- FALSE -- Instead of showing the logistic value for the y axis in response curves, show the exponent (a linear combination of features)
# linear -- boolean -- TRUE -- Allow linear features to be used
# quadratic -- boolean -- TRUE -- Allow quadratic features to be used
# product -- boolean -- TRUE -- Allow product features to be used
# threshold -- boolean -- TRUE -- Allow threshold features to be used
# hinge -- boolean -- TRUE -- Allow hinge features to be used
# addsamplestobackground -- boolean -- TRUE -- Add to the background any sample for which has a combination of environmental values that isn't already present in the background
# addallsamplestobackground -- boolean -- FALSE -- Add all samples to the background, even if they have combinations of environmental values that are already present in the background
# autorun -- boolean -- FALSE -- Start running as soon as the the program starts up
# writeplotdata -- boolean -- FALSE -- Write output files containing the data used to make response curves, for import into external plotting software
# fadebyclamping -- boolean -- FALSE -- Reduce prediction at each point in projections by the difference between
# --  --  -- clamped and non-clamped output at that point
# extrapolate -- boolean -- TRUE -- Predict to regions of environmental space outside the limits encountered during training
# visible -- boolean -- TRUE -- Make the Maxent user interface visible
# autofeature -- boolean -- TRUE -- Automatically select which feature classes to use, based on number of training samples
# doclamp -- boolean -- TRUE -- Apply clamping when projecting
# outputgrids -- boolean -- TRUE -- Write output grids. Turning this off when doing replicate runs causes only the summary grids (average, std deviation etc.) to be written, not those for the individual runs.
# plots -- boolean -- TRUE -- Write various plots for inclusion in .html output
# appendtoresultsfile -- boolean -- FALSE -- If false, maxentResults.csv file is reinitialized before each run
# maximumiterations -- integer -- 500 -- Stop training after this many iterations of the optimization algorithm
# convergencethreshold -- double -- 1.00E-05 -- Stop training when the drop in log loss per iteration drops below this number
# adjustsampleradius -- integer -- 0 -- Add this number of pixels to the radius of white/purple dots for samples on pictures of predictions.
# --  --  -- Negative values reduce size of dots.
# threads -- integer -- 1 -- Number of processor threads to use. Matching this number to the number of cores on your computer speeds up some operations, especially variable jackknifing.
# lq2lqptthreshold -- integer -- 80 -- Number of samples at which product and threshold features start being used
# l2lqthreshold -- integer -- 10 -- Number of samples at which quadratic features start being used
# hingethreshold -- integer -- 15 -- Number of samples at which hinge features start being used
# beta_threshold -- double -- -1 -- Regularization parameter to be applied to all threshold features; negative value enables automatic setting
# beta_categorical -- double -- -1 -- Regularization parameter to be applied to all categorical features; negative value enables automatic setting
# beta_lqp -- double -- -1 -- Regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
# beta_hinge -- double -- -1 -- Regularization parameter to be applied to all hinge features; negative value enables automatic setting
# logfile -- string -- maxent.log -- File name to be used for writing debugging information about a run in output directory
# cache -- boolean -- TRUE -- Make a .mxe cached version of ascii files, for faster access
# defaultprevalence -- double -- 0.5 -- Default prevalence of the species: probability of presence at ordinary occurrence points.
# --  --  -- See Elith et al., Diversity and Distributions, 2011 for details.
# applythresholdrule -- string --  -- Apply a threshold rule, generating a binary output grid in addition to the regular prediction grid. Use the full name of the threshold rule in Maxent's html output as the argument. For example, 'applyThresholdRule=Fixed cumulative value 1'.
# togglelayertype -- string --  -- Toggle continuous/categorical for environmental layers whose names begin with this prefix (default: all continuous)
# togglespeciesselected -- string --  -- Toggle selection of species whose names begin with this prefix (default: all selected)
# togglelayerselected -- string --  -- Toggle selection of environmental layers whose names begin with this prefix (default: all selected)
# verbose -- boolean -- FALSE -- Gived detailed diagnostics for debugging
# allowpartialdata -- boolean -- FALSE -- During model training, allow use of samples that have nodata values for one or more environmental variables.
# prefixes -- boolean -- TRUE -- When toggling samples or layers selected or layer types, allow toggle string to be a prefix rather than an exact match.
# nodata -- integer -- -9999 -- Value to be interpreted as nodata values in SWD sample data

if (model.maxent) {
	outdir = paste(wd,'output_maxent/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	write.csv(data.frame(species=species,occur),paste(outdir,"occur.csv",sep=''),row.names=FALSE)### create occur.csv for maxent
	write.csv(data.frame(species="bkgd",bkgd),paste(outdir,"bkgd.csv",sep=''),row.names=FALSE)### create bkgd.csv for maxent
	###not user modified section
	tstr = paste('java -mx2048m -jar ',maxent.jar,' ',sep='') #start the maxent string
	tstr = paste(tstr,'environmentallayers=',outdir,'bkgd.csv ',sep='')
	tstr = paste(tstr,'samplesfile=',outdir,'occur.csv ',sep='')
	tstr = paste(tstr,'outputdirectory=',outdir,' ',sep='')
	tstr = paste(tstr,'autorun=TRUE visible=FALSE warnings=FALSE tooltips=FALSE ',sep='')
	tstr = paste(tstr,'askoverwrite=FALSE skipifexists=FALSE prefixes=TRUE verbose=FALSE ',sep='')
	tstr = paste(tstr,'responsecurves=TRUE pictures=TRUE jackknife=TRUE writeclampgrid=TRUE ',sep='')
	tstr = paste(tstr,'writemess=TRUE writebackgroundpredictions=TRUE writeplotdata=FALSE outputgrids=TRUE ',sep='')
	tstr = paste(tstr,'plots=TRUE appendtoresultsfile=FALSE threads=1 adjustsampleradius=0 ',sep='')
	tstr = paste(tstr,'logfile=maxent.log cache=FALSE allowpartialdata=FALSE outputfiletype="asc" ',sep='')
	tstr = paste(tstr,'perspeciesresults=FALSE responsecurvesexponent=FALSE	dontcache nocache ',sep='')
	if (any(enviro.data.type!='continuous')){
		catvals = which(enviro.data.type!='continuous')
		for (ii in catvals) {
			tstr = paste(tstr,'togglelayertype=',enviro.data.names[ii],' ',sep='') #toggle the layer type
		}
	}
	### based on user modified
	tstr = paste(tstr,'outputformat=',outputformat,' ',sep='')
	tstr = paste(tstr,'randomseed=',randomseed,' ',sep='')
	tstr = paste(tstr,'logscale=',logscale,' ',sep='')
	tstr = paste(tstr,'removeduplicates=',removeduplicates,' ',sep='')
	tstr = paste(tstr,'randomtestpoints=',randomtestpoints,' ',sep='')
	tstr = paste(tstr,'betamultiplier=',betamultiplier,' ',sep='')
	tstr = paste(tstr,'maximumbackground=',maximumbackground,' ',sep='')
	tstr = paste(tstr,'biasfile=',biasfile,' ',sep='')
	tstr = paste(tstr,'testsamplesfile=',testsamplesfile,' ',sep='')
	tstr = paste(tstr,'replicates=',replicates,' ',sep='')
	tstr = paste(tstr,'replicatetype=',replicatetype,' ',sep='')
	tstr = paste(tstr,'linear=',linear,' ',sep='')
	tstr = paste(tstr,'quadratic=',quadratic,' ',sep='')
	tstr = paste(tstr,'product=',product,' ',sep='')
	tstr = paste(tstr,'threshold=',threshold,' ',sep='')
	tstr = paste(tstr,'hinge=',hinge,' ',sep='')
	tstr = paste(tstr,'addsamplestobackground=',addsamplestobackground,' ',sep='')
	tstr = paste(tstr,'addallsamplestobackground=',addallsamplestobackground,' ',sep='')
	tstr = paste(tstr,'fadebyclamping=',fadebyclamping,' ',sep='')
	tstr = paste(tstr,'extrapolate=',extrapolate,' ',sep='')
	tstr = paste(tstr,'autofeature=',autofeature,' ',sep='')
	tstr = paste(tstr,'doclamp=',doclamp,' ',sep='')
	tstr = paste(tstr,'maximumiterations=',maximumiterations,' ',sep='')
	tstr = paste(tstr,'convergencethreshold=',convergencethreshold,' ',sep='')
	tstr = paste(tstr,'lq2lqptthreshold=',lq2lqptthreshold,' ',sep='')
	tstr = paste(tstr,'l2lqthreshold=',l2lqthreshold,' ',sep='')
	tstr = paste(tstr,'hingethreshold=',hingethreshold,' ',sep='')
	tstr = paste(tstr,'beta_threshold=',beta_threshold,' ',sep='')
	tstr = paste(tstr,'beta_categorical=',beta_categorical,' ',sep='')
	tstr = paste(tstr,'beta_lqp=',beta_lqp,' ',sep='')
	tstr = paste(tstr,'beta_hinge=',beta_hinge,' ',sep='')
	tstr = paste(tstr,'defaultprevalence=',defaultprevalence,' ',sep='')
	tstr = paste(tstr,'nodata=',nodata,' ',sep='')
	system(tstr)	
}


############### BIOMOD2 Models ###############
# 1. Format the data
# 2. Define the model options
# 3. Compute the model
# NOTE: Model evaluation is included as part of model creation

# BIOMOD_FormatingData(resp.var, expl.var, resp.xy = NULL, resp.name = NULL, eval.resp.var = NULL, 
#	eval.expl.var = NULL, eval.resp.xy = NULL, PA.nb.rep = 0, PA.nb.absences = 1000, PA.strategy = 'random',
#	PA.dist.min = 0, PA.dist.max = NULL, PA.sre.quant = 0.025, PA.table = NULL, na.rm = TRUE)
#
# resp.var a vector, SpatialPointsDataFrame (or SpatialPoints if you work with ‘only presences’ data) containing species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminated ) that will be used to build the species distribution models.
# expl.var a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to build your models.
# resp.xy optional 2 columns matrix containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be used to build your models.
# eval.resp.var	a vector, SpatialPointsDataFrame your species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminated ) that will be used to evaluate the models with independant data (or past data for instance).
# eval.expl.var	a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to evaluate the models with independant data (or past data for instance).
# eval.resp.xy opional 2 columns matrix containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be used to evaluate the modelswith independant data (or past data for instance).
# resp.name	response variable name (character). The species name.
# PA.nb.rep	number of required Pseudo Absences selection (if needed). 0 by Default.
# PA.nb.absences number of pseudo-absence selected for each repetition (when PA.nb.rep > 0) of the selection (true absences included)
# PA.strategy strategy for selecting the Pseudo Absences (must be ‘random’, ‘sre’, ‘disk’ or ‘user.defined’)
# PA.dist.min minimal distance to presences for ‘disk’ Pseudo Absences selection (in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )
# PA.dist.max maximal distance to presences for ‘disk’ Pseudo Absences selection(in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )
# PA.sre.quant quantile used for ‘sre’ Pseudo Absences selection
# PA.table a matrix (or a data.frame) having as many rows than resp.var values. Each column correspund to a Pseudo-absences selection. It contains TRUE or FALSE indicating which values of resp.var will be considered to build models. It must be used with ‘user.defined’ PA.strategy.
# na.rm	logical, if TRUE, all points having one or several missing value for environmental data will be removed from analysis

# format the data as required by the biomod package
formatBiomodData = function() {
	biomod.data = rbind(occur[,c("lon","lat")],bkgd[,c("lon","lat")])
	biomod.data.pa = c(rep(1,nrow(occur)),rep(0,nrow(bkgd)))
	myBiomodData <- BIOMOD_FormatingData(resp.var = biomod.data.pa, expl.var = stack(enviro.data),	
		resp.xy = biomod.data, resp.name = species)
	return(myBiomodData)
}

	
# BIOMOD_Modeling(data, models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT'), models.options = NULL, 
#	NbRunEval=1, DataSplit=100, Yweights=NULL, Prevalence=NULL, VarImport=0, models.eval.meth = c('KAPPA','TSS','ROC'), 
#	SaveObj = TRUE, rescal.all.models = TRUE, do.full.models = TRUE, modeling.id = as.character(format(Sys.time(), '%s')),
#	...)
#
# data	BIOMOD.formated.data object returned by BIOMOD_FormatingData
# models vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF' and 'MAXENT'
# models.options BIOMOD.models.options object returned by BIOMOD_ModelingOptions
# NbRunEval	Number of Evaluation run
# DataSplit	% of data used to calibrate the models, the remaining part will be used for testing
# Yweights response points weights
# Prevalence either NULL (default) or a 0-1 numeric used to build 'weighted response weights'
# VarImport	Number of permutation to estimate variable importance
# models.eval.meth vector of names of evaluation metric among 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI' and 'ETS'
# SaveObj keep all results and outputs on hard drive or not (NOTE: strongly recommended)
# rescal.all.models	if true, all model prediction will be scaled with a binomial GLM
# do.full.models if true, models calibrated and evaluated with the whole dataset are done
# modeling.id character, the ID (=name) of modeling procedure. A random number by default.
# ... further arguments :
# DataSplitTable : a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a 'RUN'. If filled, args NbRunEval, DataSplit and do.full.models will be ignored.


		
###############
#
# GLM - generalized linear model (glm)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic', interaction.level = 0, myFormula = NULL, 
#	test = 'BIC', family = 'binomial', control = glm.control(epsilon = 1e-08, maxit = 1000, trace = FALSE)))				  
# myFormula	: a typical formula object (see example). If not NULL, type and interaction.level args are switched off
#	You can choose to either: 
#	1) generate automatically the GLM formula by using the type and interaction.level arguments 
#		type : formula given to the model ('simple', 'quadratic' or 'polynomial')
#		interaction.level : integer corresponding to the interaction level between variables considered. Consider that 
#			interactions quickly enlarge the number of effective variables used into the GLM
#	2) or construct specific formula
# test : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC'). 'none' is also a supported value which implies to concider only the full model (no stepwise selection). This can lead to convergence issu and strange results.
# family : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See family for details of family functions.) BIOMOD only runs on presence-absence data so far, so binomial family by default.
# control : a list of parameters for controlling the fitting process. For glm.fit this is passed to glm.control
#	glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
#		epsilon	- positive convergence tolerance e; the iterations converge when |dev - dev_{old}|/(|dev| + 0.1) < e
#		maxit - integer giving the maximal number of IWLS iterations
#		trace - logical indicating if output should be produced for each iteration

if (model.glm) {
	outdir = paste(wd,'output_glm/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(GLM = glm.myBiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.glm <- BIOMOD_Modeling(data = myBiomodData, models=c('GLM'), models.options= myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth=biomod.models.eval.meth, SaveObj=TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	# model output saved as part of BIOMOD_Modeling() # EMG not sure how to retrieve
	if (!is.null(myBiomodModelOut.glm)) {		
			save(myBiomodModelOut.glm, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}
	
###############
#
# GAM - generalized additive model (gam::gam, mgcv::gam or mgcv::bam)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(GAM = list(algo = 'GAM_mgcv', type = 's_smoother', k = NULL, 
#	interaction.level = 0, myFormula = NULL, family = 'binomial', 
#	control = gam.control(epsilon = 1e-06, trace = FALSE, maxit = 100)))
# algo : either "GAM_mgcv" (default), "GAM_gam" or "BAM_mgcv" defining the chosen GAM function (see gam, gam resp. bam for more details)
# myFormula : a typical formula object (see example). If not NULL, type and interaction.level args are switched off. 
#	You can choose to either:
#	1) generate automatically the GAM formula by using the type and interaction.level arguments 
#		type : the smother used to generate the formula. Only "s_smoother" available at time. 
#		interaction.level : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GAM. Interaction are not considered if you choosed "GAM_gam" algo
#	2) or construct specific formula
# k : a smooth term in a formula argument to gam (see gam s or mgcv s)
# family : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See family for details of family functions.) BIOMOD only runs on presence-absence data so far, so binomial family by default.
# control : see gam::gam.control or mgcv::gam.control
# 	gam::gam.control(epsilon=1e-07, bf.epsilon = 1e-07, maxit=30, bf.maxit = 30, trace=FALSE,...)
# 	mgcv::gam.control(irls.reg=0.0,epsilon = 1e-06, maxit = 100, mgcv.tol=1e-7,mgcv.half=15, trace = FALSE, rank.tol=.Machine$double.eps^0.5, nlm=list(), optim=list(), newton=list(), outerPIsteps=0, idLinksBases=TRUE, scalePenalty=TRUE, keepData=FALSE) 

if (model.gam) {
	outdir = paste(wd,'output_gam/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(GAM = gam.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.gam <- BIOMOD_Modeling(data = myBiomodData, models = c('GAM'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.gam)) {		
			save(myBiomodModelOut.gam, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}
	
###############
#
# GBM - generalized boosting model (usually called boosted regression trees) (gbm)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(GBM = list(distribution = 'bernoulli', interaction.depth = 7, 
#	shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, n.trees = 500, cv.folds = 5))
# n.trees : the total number of trees to fit. This is equivalent to the number of iterations and the number of basis functions in the additive expansion.
# cv.folds : Number of cross-validation folds to perform. If cv.folds>1 then gbm, in addition to the usual fit, will perform a cross-validation, calculate an estimate of generalization error returned in cv.error.
# distribution : either a character string specifying the name of the distribution to use or a list with a component name specifying the distribution and any additional parameters needed. 
#	If not specified, gbm will try to guess: if the response has only 2 unique values, bernoulli is assumed; otherwise, if the response is a factor, multinomial is assumed; otherwise, if the response has class "Surv", coxph is assumed; otherwise, gaussian is assumed.
# 	Currently available options are "gaussian" (squared error), "laplace" (absolute loss), "tdist" (t-distribution loss), "bernoulli" (logistic regression for 0-1 outcomes), "huberized" (huberized hinge loss for 0-1 outcomes), "multinomial" (classification when there are more than 2 classes), "adaboost" (the AdaBoost exponential loss for 0-1 outcomes), "poisson" (count outcomes), "coxph" (right censored observations), "quantile", or "pairwise" (ranking measure using the LambdaMart algorithm).
# 	If quantile regression is specified, distribution must be a list of the form list(name="quantile",alpha=0.25) where alpha is the quantile to estimate. The current version's quantile regression method does not handle non-constant weights and will stop.
# 	If "tdist" is specified, the default degrees of freedom is 4 and this can be controlled by specifying distribution=list(name="tdist", df=DF) where DF is your chosen degrees of freedom.
# 	If "pairwise" regression is specified, distribution must be a list of the form list(name="pairwise",group=...,metric=...,max.rank=...) (metric and max.rank are optional, see below). 
#		group is a character vector with the column names of data that jointly indicate the group an instance belongs to (typically a query in Information Retrieval applications). For training, only pairs of instances from the same group and with different target labels can be considered. 
#		metric is the IR measure to use, one of:
#			conc - Fraction of concordant pairs; for binary labels, this is equivalent to the Area under the ROC Curve
#			mrr - Mean reciprocal rank of the highest-ranked positive instance
#			map - Mean average precision, a generalization of mrr to multiple positive instances
# 			ndcg - Normalized discounted cumulative gain. The score is the weighted sum (DCG) of the user-supplied target values, weighted by log(rank+1), and normalized to the maximum achievable value. This is the default if the user did not specify a metric.
#				ndcg and conc allow arbitrary target values, while binary targets {0,1} are expected for conc and ndcg. 
#		For ndcg and mrr, a cut-off can be chosen using a positive integer parameter max.rank. If left unspecified, all ranks are taken into account.
#	Note that splitting of instances into training and validation sets follows group boundaries and therefore only approximates the specified train.fraction ratio (the same applies to cross-validation folds). Internally, queries are randomly shuffled before training, to avoid bias.
#	Weights can be used in conjunction with pairwise metrics, however it is assumed that they are constant for instances from the same group.
#	For details and background on the algorithm, see e.g. Burges (2010).
# interaction.depth : The maximum depth of variable interactions. 1 implies an additive model, 2 implies a model with up to 2-way interactions, etc.
# shrinkage : a shrinkage parameter applied to each tree in the expansion. Also known as the learning rate or step-size reduction.
# bag.fraction : the fraction of the training set observations randomly selected to propose the next tree in the expansion. This introduces randomnesses into the model fit. 
#	If bag.fraction<1 then running the same model twice will result in similar but different fits. gbm uses the R random number generator so set.seed can ensure that the model can be reconstructed. 
#	Preferably, the user can save the returned gbm.object using save.
# train.fraction : The first train.fraction * nrows(data) observations are used to fit the gbm and the remainder are used for computing out-of-sample estimates of the loss function.

if (model.gbm) {
	outdir = paste(wd,'output_gbm/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(GBM = gbm.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.gbm <- BIOMOD_Modeling(data = myBiomodData, models = c('GBM'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth=biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models= biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.gbm)) {		
			save(myBiomodModelOut.gbm, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
#
# CTA - classification tree analysis (rpart)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(CTA = list(method = 'class', parms = 'default', cost = NULL, 
#		control = rpart.control(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25))
# method : one of "anova", "poisson", "class" or "exp". If method is missing then the routine tries to make an intelligent guess. 
#	If y is a survival object, then method = "exp" is assumed, if y has 2 columns then method = "poisson" is assumed, if y is a factor then method = "class" is assumed, otherwise method = "anova" is assumed. 
#	It is wisest to specify the method directly, especially as more criteria may added to the function in future. 
# parms : optional parameters for the splitting function.
#	Anova splitting has no parameters.
# 	Poisson splitting has a single parameter, the coefficient of variation of the prior distribution on the rates. The default value is 1.
#	Exponential splitting has the same parameter as Poisson.
#	For classification splitting, the list can contain any of: the vector of prior probabilities (component prior), the loss matrix (component loss) or the splitting index (component split). 
#		The priors must be positive and sum to 1. The loss matrix must have zeros on the diagonal and positive off-diagonal elements. The splitting index can be gini or information. 
#		The default priors are proportional to the data counts, the losses default to 1, and the split defaults to gini.
# cost : a vector of non-negative costs, one for each variable in the model. Defaults to one for all variables. 
#	These are scalings to be applied when considering splits, so the improvement on splitting on a variable is divided by its cost in deciding which split to choose.
# control : a list of options that control details of the rpart algorithm. See rpart.control.
#	rpart.control(minsplit = 20, minbucket = round(minsplit/3), cp = 0.01, maxcompete = 4, maxsurrogate = 5, 
#		usesurrogate = 2, xval = 10, surrogatestyle = 0, maxdepth = 30, ...)
# NOTE: for method and parms, you can give a 'real' value as described in the rpart help file or 'default' that implies default rpart values.

if (model.cta) {
	outdir = paste(wd,'output_cta/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(CTA = cta.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.cta <- BIOMOD_Modeling(data = myBiomodData, models = c('CTA'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.cta)) {		
			save(myBiomodModelOut.cta, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
# 
# ANN - artificial neural network (nnet)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(ANN = list(NbCV = 5, rang = 0.1, maxit = 200))
# NbCV : nb of cross validation to find best size and decay parameters
# rang : Initial random weights on [-rang, rang]
# maxit : maximum number of iterations. Default 100

if (model.ann) {
	outdir = paste(wd,'output_ann/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(ANN = ann.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.ann <- BIOMOD_Modeling(data = myBiomodData, models = c('ANN'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.ann)) {		
			save(myBiomodModelOut.ann, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
#
# SRE - surface range envelop (usually called BIOCLIM)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(SRE = list(quant = 0.025))
# quant : quantile of 'extreme environmental variable' removed for selection of species envelops

if (model.sre) {
	outdir = paste(wd,'output_sre/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(SRE = sre.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.sre <- BIOMOD_Modeling(data = myBiomodData, models = c('SRE'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.sre)) {		
			save(myBiomodModelOut.sre, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
#
# FDA - flexible discriminant analysis (fda)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(FDA = list(method = 'mars'))
# method : regression method used in optimal scaling. 
#	Default is linear regression via the function polyreg, resulting in linear discriminant analysis. 
#	Other possibilities are mars and bruto. For Penalized Discriminant analysis gen.ridge is appropriate.

if (model.fda) {
	outdir = paste(wd,'output_fda/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(FDA = fda.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.fda <- BIOMOD_Modeling(data = myBiomodData, models = c('FDA'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.fda)) {		
			save(myBiomodModelOut.fda, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
#
# MARS - multiple adaptive regression splines (mars)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(MARS = list(degree = 2, penalty = 2, thresh = 0.001, prune = TRUE))
# degree : an optional integer specifying maximum interaction degree (default is 1)
# penalty : an optional value specifying the cost per degree of freedom charge (default is 2)
# thresh : an optional value specifying forward stepwise stopping threshold (default is 0.001)
# prune : an optional logical value specifying whether the model should be pruned in a backward stepwise fashion (default is TRUE)

if (model.mars) {
	outdir = paste(wd,'output_mars/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(MARS = mars.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.mars <- BIOMOD_Modeling(data = myBiomodData, models = c('MARS'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.mars)) {		
			save(myBiomodModelOut.mars, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
#
# RF - random forest (randomForest)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(RF = list(do.classif = TRUE, ntree = 50, mtry = 'default'))
# do.classif : if TRUE classification random.forest computed else regression random.forest will be done
# ntree : Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
# mtry : Number of variables randomly sampled as candidates at each split. 
#	Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3)
# NOTE: for mtry, you can give a 'real' value as described in randomForest help file or 'default' that implies default randomForest values

if (model.rf) {
	outdir = paste(wd,'output_rf/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(RF = rf.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.rf <- BIOMOD_Modeling(data = myBiomodData, models = c('RF'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.rf)) {		
			save(myBiomodModelOut.rf, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}

###############
#
# MAXENT - maxent
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(MAXENT = list(path_to_maxent.jar = 'C:/userdata', maximumiterations = 200, 
#	visible = FALSE, linear = TRUE, quadratic = TRUE, product = TRUE, threshold = TRUE, hinge = TRUE, 
#	lq2lqptthreshold = 80, l2lqthreshold = 10, hingethreshold = 15, beta_threshold = -1, beta_categorical = -1, 
#	beta_lqp = -1, beta_hinge = -1, defaultprevalence = 0.5))
# path_to_maxent.jar : character, the link to maxent.jar file (the working directory by default)
# maximumiterations : integer (default 200), maximum iteration done
# visible : logical (default FALSE), make the Maxent user interface visible
# linear : logical (default TRUE), allow linear features to be used
# quadratic : logical (default TRUE), allow quadratic features to be used
# product : logical (default TRUE), allow product features to be used
# threshold : logical (default TRUE), allow threshold features to be used
# hinge : logical (default TRUE), allow hinge features to be used
# lq2lqptthreshold : integer (default 80), number of samples at which product and threshold features start being used
# l2lqthreshold : integer (default 10), number of samples at which quadratic features start being used
# hingethreshold : integer (default 15), number of samples at which hinge features start being used
# beta_threshold : numeric (default -1.0), regularization parameter to be applied to all threshold features; negative value enables automatic setting
# beta_categorical : numeric (default -1.0), regularization parameter to be applied to all categorical features; negative value enables automatic setting
# beta_lqp : numeric (default -1.0), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
# beta_hinge : numeric (default -1.0), regularization parameter to be applied to all hinge features; negative value enables automatic setting
# defaultprevalence : numeric (default 0.5), default prevalence of the species: probability of presence at ordinary occurrence points

if (model.biomod.maxent) {
	outdir = paste(wd,'output_biomod.maxent/',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(MAXENT = biomod.maxent.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.biomod.maxent <- BIOMOD_Modeling(data = myBiomodData, models = c('MAXENT'),	models.options = myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth = biomod.models.eval.meth, SaveObj = TRUE,
		rescal.all.models = biomod.rescal.all.models, do.full.models = biomod.do.full.models, 
		modeling.id = biomod.modeling.id)
	if (!is.null(myBiomodModelOut.biomod.maxent)) {		
			save(myBiomodModelOut.biomod.maxent, file=paste(outdir,"model.object.RData",sep='')) #save out the model object
	}
}
