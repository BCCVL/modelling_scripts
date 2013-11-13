src/org/bccvl/compute/rscripts/brt.R# set CRAN mirror in case we need to download something
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
    project.brt=FALSE
} else {
    future.climate.scenario = stack(enviro.data.future)
}

###run the models and store models
#############################################################################################
#
# MACHINE LEARNING METHODS - use both presence and absence or background data: Maxent, BRT
#
#############################################################################################

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
	tstr = paste(tstr,'outputformat=',maxent.outputformat,' ',sep='')
	tstr = paste(tstr,'randomseed=',maxent.randomseed,' ',sep='')
	tstr = paste(tstr,'logscale=',maxent.logscale,' ',sep='')
	tstr = paste(tstr,'removeduplicates=',maxent.removeduplicates,' ',sep='')
	tstr = paste(tstr,'randomtestpoints=',maxent.randomtestpoints,' ',sep='')
	tstr = paste(tstr,'betamultiplier=',maxent.betamultiplier,' ',sep='')
	tstr = paste(tstr,'maximumbackground=',maxent.maximumbackground,' ',sep='')
	tstr = paste(tstr,'biasfile=',maxent.biasfile,' ',sep='')
	tstr = paste(tstr,'testsamplesfile=',maxent.testsamplesfile,' ',sep='')
	tstr = paste(tstr,'replicates=',maxent.replicates,' ',sep='')
	tstr = paste(tstr,'replicatetype=',maxent.replicatetype,' ',sep='')
	tstr = paste(tstr,'linear=',maxent.linear,' ',sep='')
	tstr = paste(tstr,'quadratic=',maxent.quadratic,' ',sep='')
	tstr = paste(tstr,'product=',maxent.product,' ',sep='')
	tstr = paste(tstr,'threshold=',maxent.threshold,' ',sep='')
	tstr = paste(tstr,'hinge=',maxent.hinge,' ',sep='')
	tstr = paste(tstr,'addsamplestobackground=',maxent.addsamplestobackground,' ',sep='')
	tstr = paste(tstr,'addallsamplestobackground=',maxent.addallsamplestobackground,' ',sep='')
	tstr = paste(tstr,'fadebyclamping=',maxent.fadebyclamping,' ',sep='')
	tstr = paste(tstr,'extrapolate=',maxent.extrapolate,' ',sep='')
	tstr = paste(tstr,'autofeature=',maxent.autofeature,' ',sep='')
	tstr = paste(tstr,'doclamp=',maxent.doclamp,' ',sep='')
	tstr = paste(tstr,'maximumiterations=',maxent.maximumiterations,' ',sep='')
	tstr = paste(tstr,'convergencethreshold=',maxent.convergencethreshold,' ',sep='')
	tstr = paste(tstr,'lq2lqptthreshold=',maxent.lq2lqptthreshold,' ',sep='')
	tstr = paste(tstr,'l2lqthreshold=',maxent.l2lqthreshold,' ',sep='')
	tstr = paste(tstr,'hingethreshold=',maxent.hingethreshold,' ',sep='')
	tstr = paste(tstr,'beta_threshold=',maxent.beta_threshold,' ',sep='')
	tstr = paste(tstr,'beta_categorical=',maxent.beta_categorical,' ',sep='')
	tstr = paste(tstr,'beta_lqp=',maxent.beta_lqp,' ',sep='')
	tstr = paste(tstr,'beta_hinge=',maxent.beta_hinge,' ',sep='')
	tstr = paste(tstr,'defaultprevalence=',maxent.defaultprevalence,' ',sep='')
	tstr = paste(tstr,'nodata=',maxent.nodata,' ',sep='')
	system(tstr)	
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

if (project.maxent) {

	#*************** UNDER CONSTRUCTION ***************
	# maxent model creation was run as a system call outside of R, need to do the same for projection
	# EMG check to see if argument defaults / modifiables are the same as during creation
	
	# create output directory
	model.dir = paste(wd, "output_maxent/", sep="")
	
	### not user modified section
	tstr = paste("java -cp ", maxent.jar, " density.Project ", model.dir, species, ".lambdas ", sep="")
	# where to find the climate scenarios
	tstr = paste(tstr, dirname(enviro.data[1]), " ", sep="")
	# where to put, what to name the output
	tstr = paste(tstr, model.dir, es.name, ".asc", sep="")
	# optional arguments
	tstr = paste(tstr, " nowriteclampgrid nowritemess fadebyclamping dontcache", sep="")
	system(tstr)
	
	# EMG cache=FALSE nocache and dontcache all manage to be ignored and a maxent.cache is created
	# EMG 'outputfiletype' = asc, mxe, grd, bil only NOT geotiff; can create *.png ('pictures=TRUE')
}


######################################################################################
# model accuracy helpers
######################################################################################

# function to save evaluate output
saveModelEvaluation = function(out.model, out.biomod.model, model.name) {
    model.dir = paste(wd, "/output_", model.name, "/", sep="")
    save(out.model, file=paste(model.dir, "dismo.eval.object.RData", sep=''))	# save the 'dismo::ModelEvalution' object

    # save all the model accuracy statistics provided in both dismo and biomod2
    rownames(out.biomod.model) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
    write.csv(t(round(out.biomod.model, digits=3)), file=paste(model.dir, "combined.modelEvaluation.csv", sep=""))
    # EMG no guarantee these value are correct

    # save AUROC curve

    png(file=paste(model.dir, "AUC.png", sep=''));
    plot(out.model, 'ROC');
    dev.off()
}

my.Find.Optim.Stat <- function(Stat='TSS',Fit,Obs,Precision = 5, Fixed.thresh = NULL){
    if(length(unique(Obs)) == 1 | length(unique(Fit)) == 1){
        # warning("\nObserved or fited data contains only a value.. Evaluation Methods switched off\n",immediate.=T)
        # best.stat <- cutoff <- true.pos <- sensibility <- true.neg <- specificity <- NA
        warning("\nObserved or fited data contains a unique value.. Be carefull with this models predictions\n",immediate.=T)
        #best.stat <- cutoff <- true.pos <- sensibility <- true.neg <- specificity <- NA
    } #else {
    if(Stat != 'ROC'){
        StatOptimum <- my.getStatOptimValue(Stat)
        if(is.null(Fixed.thresh)){ # test a range of threshold to get the one giving the best score
            if(length(unique(Fit)) == 1){
                valToTest <- unique(Fit)
                valToTest <- round(c(mean(c(0,valToTest)), mean(c(1000,valToTest))))
            } else{
                mini <- max(min(quantile(Fit,0.05, na.rm=T), na.rm=T),0)
                maxi <- min(max(quantile(Fit,0.95, na.rm=T), na.rm=T),1000)
                # valToTest <- unique( round(c(seq(mini,maxi,length.out=100), mini, maxi)) )
                # EMG no idea why the round() is here, it makes vals between 0 and 1 (ie bioclim) all 0
                valToTest <- unique( c(seq(mini,maxi,length.out=100)))
                # deal with unique value to test case
                if(length(valToTest)<3){
                    valToTest <- round(c(mean(0,mini), valToTest, mean(1000,maxi)))
                }
            }
            # valToTest <- unique( c(seq(mini,maxi,by=Precision), mini, maxi) )
        } else{
            valToTest <- Fixed.thresh
        }

        calcStat <- sapply(lapply(valToTest, function(x){return(table(Fit>x,Obs))} ), my.calculate.stat, stat=Stat)

        # scal on 0-1 ladder.. 1 is the best
        calcStat <- 1 - abs(StatOptimum - calcStat)

        best.stat <- max(calcStat, na.rm=T)

        cutoff <- median(valToTest[which(calcStat==best.stat)]) # if several values are selected

        misc <- table(Fit >= cutoff, Obs)
        misc <- .contagency.table.check(misc)
        true.pos <- misc['TRUE','1']
        true.neg <- misc['FALSE','0']
        specificity <- (true.neg * 100)/sum(misc[,'0'])
        sensibility <- (true.pos * 100)/sum(misc[,'1'])
    } else{
        roc1 <- roc(Obs, Fit, percent=T)
        roc1.out <- coords(roc1, "best", ret=c("threshold", "sens", "spec"))
        best.stat <- as.numeric(auc(roc1))/100
        cutoff <- as.numeric(roc1.out["threshold"])
        sensibility <- as.numeric(roc1.out["sensitivity"])
        specificity <- as.numeric(roc1.out["specificity"])
    }
  #}
    return(cbind(best.stat,cutoff,sensibility,specificity))
}

my.getStatOptimValue <- function(stat){
    if(stat == 'TSS') return(1)
    if(stat == 'KAPPA') return(1)
    if(stat == 'ACCURACY') return(1)
    if(stat == 'BIAS') return(1)
    if(stat == 'POD') return(1)
    if(stat == 'FAR') return(0)
    if(stat == 'POFD') return(0)
    if(stat == 'SR') return(1)
    if(stat == 'CSI') return(1)
    if(stat == 'ETS') return(1)
    if(stat == 'HK') return(1)
    if(stat == 'HSS') return(1)
    if(stat == 'OR') return(1000000)
    if(stat == 'ORSS') return(1)

    #dismo
    if(stat == 'ODP') return(1)
    # if(stat == 'CCR') return(1) # same as ACCURACY
    # if(stat == 'TPR') return(1) # same as POD
    if(stat == 'TNR') return(1)
    if(stat == 'FPR') return(0)
    if(stat == 'FNR') return(0)
    # if(stat == 'PPP') return(1) # same as SR
    if(stat == 'NPP') return(1)
    if(stat == 'MCR') return(0)
    if(stat == 'OR') return(1000000)
    # if(stat == 'kappa') return(1) # same as KAPPA
}

my.calculate.stat <-
    function(Misc, stat='TSS') {
        # Contagency table checking
        Misc <- .contagency.table.check(Misc)

        # Defining Classification index
        hits <- Misc['TRUE','1']
        misses <- Misc['FALSE','1']
        false_alarms <- Misc['TRUE','0']
        correct_negatives <- Misc['FALSE','0']

        total <- sum(Misc)
        forecast_1 <- sum(Misc['TRUE',])
        forecast_0 <- sum(Misc['FALSE',])
        observed_1 <- sum(Misc[,'1'])
        observed_0 <- sum(Misc[,'0'])

        # Calculating choosen evaluating metric
        if(stat=='TSS'){
            return( (hits/(hits+misses)) + (correct_negatives/(false_alarms+correct_negatives)) -1 )
        }

        if(stat=='KAPPA'){
            Po <- (1/total) * (hits + correct_negatives)
            Pe <- ((1/total)^2) * ((forecast_1 * observed_1) + (forecast_0 * observed_0))
            return( (Po - Pe) / (1-Pe) )
        }

        if(stat=='ACCURACY'){
            return( (hits + correct_negatives) / total)
        }

        if(stat=='BIAS'){
            return( (hits + false_alarms) / (hits + misses))
        }

        if(stat=='POD'){
            return( hits / (hits + misses))
        }

        if(stat=='FAR'){
            return(false_alarms/(hits+false_alarms))
        }

        if(stat=='POFD'){
            return(false_alarms / (correct_negatives + false_alarms))
        }

        if(stat=='SR'){
            return(hits / (hits + false_alarms))
        }

        if(stat=='CSI'){
            return(hits/(hits+misses+false_alarms))
        }

        if(stat=='ETS'){
            hits_rand <- ((hits+misses)*(hits+false_alarms)) / total
            return( (hits-hits_rand) / (hits+misses+false_alarms-hits_rand))
        }

        # if(stat=='HK'){
        # return((hits/(hits+misses)) - (false_alarms/(false_alarms + correct_negatives)))
        # }

        # if(stat=='HSS'){
        # expected_correct_rand <- (1/total) * ( ((hits+misses)*(hits+false_alarms)) +
        # ((correct_negatives + misses)*(correct_negatives+false_alarms)) )
        # return((hits+correct_negatives-expected_correct_rand) / (total - expected_correct_rand))
        # }

        # if(stat=='OR'){
        # return((hits*correct_negatives)/(misses*false_alarms))
        # }

        # if(stat=='ORSS'){
        # return((hits*correct_negatives - misses*false_alarms) / (hits*correct_negatives + misses*false_alarms))
        # }

        # if(stat=="BOYCE"){
        #
        # }

        #dismo
        if(stat=='ODP'){
            return((false_alarms + correct_negatives) / total)
        }

        # if(stat=='CCR'){
        # return((hits + correct_negatives) / total)
        # }

        # if(stat=='TPR'){
        # return(hits / (hits + misses))
        # }

        if(stat=='TNR'){
            return(correct_negatives / (false_alarms + correct_negatives))
        }

        if(stat=='FPR'){
            return(false_alarms / (false_alarms + correct_negatives))
        }

        if(stat=='FNR'){
            return(misses / (hits + misses))
        }

        # if(stat=='PPP'){
        # return(hits / (hits + false_alarms))
        # }

        if(stat=='NPP'){
            return(correct_negatives / (misses + correct_negatives))
        }

        if(stat=='MCR'){
            return((false_alarms + misses) / total)
        }

        if(stat=='OR'){
            return((hits * correct_negatives) / (misses * false_alarms))
        }

        # if(stat=='kappa'){
        # return(((hits + correct_negatives) - (((hits + misses)*(hits + false_alarms) + (false_alarms + correct_negatives)*(misses + correct_negatives)) / total)) /
        # (total -(((hits + misses)*(hits + false_alarms) + (false_alarms + correct_negatives)*(misses + correct_negatives)) / total)))
        # }
    }

.contagency.table.check <- function(Misc){
    # Contagency table checking
    if(dim(Misc)[1]==1){
        if(row.names(Misc)[1]=="FALSE"){
            Misc <- rbind(Misc, c(0,0))
            rownames(Misc) <- c('FALSE','TRUE')
        } else{
            a <- Misc
            Misc <- c(0,0)
            Misc <- rbind(Misc, a)
            rownames(Misc) <- c('FALSE','TRUE')
        }
    }

    if(ncol(Misc) != 2 | nrow(Misc) !=2 ){
        Misc = matrix(0, ncol=2, nrow=2, dimnames=list(c('FALSE','TRUE'), c('0','1')))
    }

    if((sum(colnames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) | (sum(rownames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) ){
        stop("Unavailable contagency table given")
    }

    if('0' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='0')] <- 'FALSE'
    if('1' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='1')] <- 'TRUE'

    return(Misc)
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
if (evaluate.maxent) {
	# read in the Maxent predictions at the presence and background points, and 
	#	extract the columns we need
	model.dir <- paste(wd, "output_maxent", sep=""); setwd(model.dir);
	presence <- read.csv(paste(model.dir, "/", species, "_samplePredictions.csv", sep=""))
	background <- read.csv(paste(model.dir, "/", species, "_backgroundPredictions.csv", sep=""))
	log.presence <- presence$Logistic.prediction
	log.absence <- background$logistic
	maxent.eval.obj = dismoModelEvaluation(log.presence, log.absence) # use predictions to generate dismo-like model evaluation object
		
	# need predictions and observed values to create confusion matrices for accuracy statistics
	maxent.fit = c(log.presence, log.absence)
	maxent.obs = c(rep(1, length(log.presence)), rep(0, length(log.absence)))

	# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
	maxent.combined.eval = sapply(model.accuracy, function(x){
			return(my.Find.Optim.Stat(Stat = x, Fit = maxent.fit, Obs = maxent.obs))
		})
	saveModelEvaluation(maxent.eval.obj, maxent.combined.eval)	# save output
	rm(list=c("maxent.eval.obj", "maxent.combined.eval")); #clean up the memory
}