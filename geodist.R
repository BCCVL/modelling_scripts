src/org/bccvl/compute/rscripts/geodist.R# set CRAN mirror in case we need to download something
# TODO: this should be done on demand or on user basis...
r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# TODO: alse creating and populating add on package location is something that should not be done system wide

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("dismo","SDMTools", "gbm", "rgdal", "pROC") #list the libraries needed
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
    && file.exists(paste(wd, "/lbkgd.RData", sep=""))) {
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

## Needed for tryCatch'ing:
err.null <- function (e) return(NULL)

# function to save projection output raster
saveModelProjection = function(out.model, model.name, projectiontime) {
    model.dir = paste(wd, "/output_", model.name, "/", sep="")
    filename = paste(projectiontime, '.tif')
    writeRaster(out.model, paste(model.dir, projectiontime, sep="/"), format="GTiff")
}

# function to get model object
getModelObject = function(model.name) {
    model.dir = paste(wd, "/output_", model.name, "/", sep="")
    model.obj = tryCatch(get(load(file=paste(model.dir, "model.object.RData", sep=""))), error = err.null)
}

###run the models and store models
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
		gd.proj = predict(gd, current.climate.scenario, fun=opt.fun, scale=opt.scale, ext=opt.ext)   # predict for given climate scenario
		saveModelProjection(gd.proj, "geodist", "current") # save output
		rm(list=c("gd", "gd.proj")); #clean up memory
	} else {
		write(paste("FAIL!", species, "Cannot create geodist model object", sep=": "), stdout())
	} # end if null
} # end if

# functions for assitance with projections
# function to check that the environmental layers used to project the  model are the same as the ones used
# 	to create the model object 
checkModelLayers = function(model.obj) {

	message("Checking environmental layers used for projection")
	# get the names of the environmental layers from the original model
	if (inherits(model.obj, "DistModel")) { # dismo package
		model.layers = colnames(model.obj@presence)
	} else if (inherits(model.obj, "gbm")) { # brt package
		model.layers = summary(brt.obj)$var
	} else if (inherits(model.obj, "BIOMOD.models.out")) { # biomod package
		model.layers = model.obj@expl.var.names
	}
	
	# get the names of the climater scenario's env layers
	pred.layers = names(future.climate.scenario)
	
	# check if the env layers were in the original model
    if(sum(!(pred.layers %in% model.layers)) > 0 ){
		message("Dropping environmental layers not used in the original model creation...")
		# create a new list of env predictors by dropping layers not in the original model
		new.predictors = climate.scenario
		for (pl in pred.layers) {
			if (!(pl %in% model.layers)) {
				new.predictors = dropLayer(new.predictors, pl)
			}	
		}
		return(new.predictors)
	} else {
		return(future.climate.scenario)
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
if (project.geodist) {
	geodist.obj = getModelObject("geodist") # get the model object
	if (!is.null(geodist.obj)) {
		predictors = checkModelLayers(geodist.obj)
		geodist.proj = predict(geodist.obj, predictors, scale=opt.scale, ext=opt.ext) # predict for given climate scenario
		saveModelProjection(geodist.proj, "geodist") 	# save output
		rm(list=c("geodist.obj", "geodist.proj")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load geodist.obj from", wd, "output_geodist", sep=": "), stdout())
	}
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

# function to generate marginal (mean) response curves for dismo models
# i.e., hold all but one predictor variable to its mean value and recalculate model predictions
createMarginalResponseCurves = function(out.model, model.name) {
	model.dir = paste(getwd(), "/", sep="")

	# get the enviromental variables and values used to create the model
	if (model.name == "brt") {
		model.values = matrix(out.model$data$x, ncol=length(out.model$var.names)); env.vars = out.model$var.names;
	} else if (model.name %in% c("geoIDW", "voronoiHull")) {
		model.values = rbind(out.model@presence, out.model@absence); env.vars = colnames(model.values);
	} else {
		model.values = out.model@presence; env.vars = colnames(model.values);
	}

	if (!(length(model.values)==0)) {

		# create a matrix to hold average values for each environmental variable
		mean.values = matrix(data = NA, nrow = 100, ncol = length(env.vars)); colnames(mean.values) = env.vars;
		# for each variable, populate the column with the mean value
		for (i in 1:ncol(mean.values)) {
			mean.values[,i] = rep(mean(model.values[,i], na.rm=TRUE), 100)
		}
	
		# allow each environmental variable to vary, keeping other variable values at average, and predict suitability
		for (j in 1:ncol(mean.values)) {
			range.values = seq(min(model.values[,j]), max(model.values[,j]), length.out=100)
			temp.data = mean.values
			temp.data[,j] = range.values
			if (model.name == "brt") {
				colnames(temp.data) = env.vars
				new.predictions = predict(out.model, as.data.frame(temp.data), n.trees = out.model$gbm.call$best.trees, type = "response")
			} else {
				new.predictions = predict(out.model, temp.data)
			}
			
			# create separate file for each response curve
			save.name = env.vars[j]
			png(file=paste(model.dir, save.name, "_response.png", sep=""))
				plot(range.values, new.predictions, ylim=c(0,1), xlab="", ylab="", main=save.name, type="l")
				rug(model.values[,j])
			dev.off()
		}
	} else {
		write(paste(species, ": Cannot create response curves from", model.name, "object", sep=" "), stdout())
	}
}

# function to calculate variable importance values for dismo models based on biomod2's correlation between predictions
# i.e., hold all but one predictor variable to its actual values, resample that one predictor and recalculate model predictions
calculateVariableImpt = function(out.model, model.name, num_samples) {
# EMG num_samples should be same as biomod.VarImport arg set in 01.init.args.model.current.R 
	model.dir = paste(getwd(), "/", sep="")
	
	# get the enviromental variables and values used to create the model
	# EMG this is duplicated from above, should be able to combine
	if (model.name == "brt") {
		model.values = matrix(out.model$data$x, ncol=length(out.model$var.names)); env.vars = out.model$var.names;
		colnames(model.values) = env.vars
	} else if (model.name %in% c("geoIDW", "voronoiHull")) {
		model.values = rbind(out.model@presence, out.model@absence); env.vars = colnames(model.values);
	} else {
		model.values = out.model@presence; env.vars = colnames(model.values);
	}

	if (!(length(model.values)==0)) {
	
		# predict using actual values
		if (model.name == "brt") {
			actual.predictions = predict(out.model, as.data.frame(model.values), n.trees = out.model$gbm.call$best.trees, type = "response")
		} else {
			actual.predictions = predict(out.model, model.values)
		} 

		# create a table to hold the output
		varimpt.out = matrix(NA, nrow=length(env.vars), ncol=num_samples+2)
		dimnames(varimpt.out) = list(env.vars, c(paste("sample_", c(1:num_samples, "mean")), "percent"))
		
		# create a copy of the env data matrix
		sample.data = model.values
		
		# for each predictor variable 
		for (p in 1:ncol(sample.data)) {

			# for each num_sample
			for (s in 1:num_samples) {
					
				# resample from that variables' values, keeping other variable values the same, and predict suitability
				sample.data[,p] = sample(x=sample.data[,p], replace=FALSE)

				# predict using sampled values
				if (model.name == "brt") {
					new.predictions = predict(out.model, as.data.frame(sample.data), n.trees = out.model$gbm.call$best.trees, type = "response")
				} else {
					new.predictions = predict(out.model, sample.data)
				}
			
				# calculate correlation between original predictions and new predictions
				varimpt.out[p,s] = 1-max(round(cor(x=actual.predictions, y=new.predictions, use="pairwise.complete.obs",
					method="pearson"), digits=3),0)
			}		
		}
		
		# calculate mean variable importance, normalize to percentages, and write results
		varimpt.out[,num_samples+1] = round(rowMeans(varimpt.out, na.rm=TRUE), digits=3)
		varimpt.out[,num_samples+2] = round((varimpt.out[,num_samples+1]/sum(varimpt.out[,num_samples+1]))*100, digits=0)
		write.csv(varimpt.out, file=paste(getwd(), "/biomod2_like_VariableImportance.csv", sep=""))
	} else {
		write(paste(species, ": Cannot calculate variable importance for ", model.name, "object", sep=" "), stdout())
	}
}

# function to calculate variable importance values for dismo models based on Maxent's decrease in AUC 
# i.e., hold all but one predictor variable to its original values, resample that one predictor and recalculate model AUC
calculatePermutationVarImpt = function(out.model, model.eval, model.name) {

	model.dir = paste(getwd(), "/", sep="")
	
	# get the enviromental variables and values used to create the model
	# EMG this is duplicated from above, should be able to combine or find an easier way to determine
	if (model.name == "brt") {
		model.values = matrix(out.model$data$x, ncol=length(out.model$var.names)); env.vars = out.model$var.names;
		colnames(model.values) = env.vars
	} else if (model.name %in% c("geoIDW", "voronoiHull")) {
		model.values = rbind(out.model@presence, out.model@absence); env.vars = colnames(model.values);
	} else {
		model.values = out.model@presence; env.vars = colnames(model.values);
	}

	if (!(length(model.values)==0)) {
	
		# get the occurrence and background environmental data used to evaluate the model
		p.swd=occur
		a.swd=bkgd
			
		# get the AUC from the original model evaluation
		init.auc = round(model.eval@auc, digits=3)
		
		# create a table to hold the output
		permvarimpt.out = matrix(NA, nrow=length(env.vars), ncol=4)
		dimnames(permvarimpt.out) = list(env.vars, c("init.auc", "sample.auc", "change.auc", "percent"))
		permvarimpt.out[,"init.auc"] = rep(init.auc, length(env.vars))
		
		# create a copy of the occurrence and background environmental data
		sample.p = p.swd[,env.vars]
		sample.a = a.swd[,env.vars]
			
		# for each predictor variable 
		for (v in 1:length(env.vars)) {
					
			# resample from that variables' values, keeping other variable values the same, and 
			sample.p[,v] = sample(x=sample.p[,v], replace=FALSE)
			sample.a[,v] = sample(x=sample.a[,v], replace=FALSE)

			# re-evaluate model with sampled env values
			if (model.name == "brt") {
				sample.eval = evaluate(p=sample.p, a=sample.a, model=brt.obj, n.trees=brt.obj$gbm.call$best.trees)
			} else {
				sample.eval = evaluate(p=sample.p, a=sample.a, model=out.model)
			}
			# get the new auc
			permvarimpt.out[v,"sample.auc"] = round(sample.eval@auc, digits=3)
		}
		
		# calculate the difference in auc, normalize to percentages, and write results
		permvarimpt.out[,"change.auc"] = permvarimpt.out[,"init.auc"] - permvarimpt.out[,"sample.auc"]
		for (r in 1:nrow(permvarimpt.out)) {
			if (permvarimpt.out[r,"change.auc"] < 0) {  # EMG what if AUC increases?
				permvarimpt.out[r,"change.auc"] = 0
			}
		}
		permvarimpt.out[,"percent"] = round((permvarimpt.out[,"change.auc"]/sum(permvarimpt.out[,"change.auc"]))*100, digits=0)
		write.csv(permvarimpt.out, file=paste(getwd(), "/maxent_like_VariableImportance.csv", sep=""))
	} else {
		write(paste(species, ": Cannot calculate maxent-like variable importance for ", model.name, "object", sep=" "), stdout())
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
if (evaluate.geodist) {
	geodist.obj = getModelObject("geodist") # get the model object
	if (!is.null(geodist.obj)) {
		geodist.eval = evaluate(model=geodist.obj, p=occur, a=bkgd) # evaluate model using dismo's evaluate
		#EMG NOTE: no error for p,a if columns not specified c.f. convHull
		
		# need predictions and observed values to create confusion matrices for accuracy statistics
		geodist.fit = c(geodist.eval@presence, geodist.eval@absence)
		geodist.obs = c(rep(1, length(geodist.eval@presence)), rep(0, length(geodist.eval@absence)))

		# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
		geodist.combined.eval = sapply(model.accuracy, function(x){
			return(my.Find.Optim.Stat(Stat = x, Fit = geodist.fit, Obs = geodist.obs))
		})
		saveModelEvaluation(geodist.eval, geodist.combined.eval)	# save output
				
		# create response curves
		createMarginalResponseCurves(geodist.obj, "geodist")

		# calculate variable importance (like biomod2, using correlations between predictions)
		calculateVariableImpt(geodist.obj, "geodist", 3)
		
		# calculate variable importance (like maxent, using decrease in AUC)
		calculatePermutationVarImpt(geodist.obj, geodist.eval, "geodist")
		
		rm(list=c("geodist.obj", "geodist.eval", "geodist.combined.eval")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load geodist.obj from", wd, "output_geodist", sep=": "), stdout())
	}
}