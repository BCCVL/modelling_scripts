#helper functions for BCCVL model, project, and evaulate scripts

## Needed for tryCatch'ing:
err.null <- function (e) return(NULL)

# function to get model object
getModelObject = function(model.name) {
	model.dir = paste(wd, "/output_", model.name, sep="")
	model.obj = tryCatch(get(load(file=paste(model.dir, "/model.object.RData", sep=""))), error = err.null)	
}

# function to save projection output raster
saveModelProjection = function(out.model, model.name, projectiontime) {
    model.dir = paste(wd, "/output_", model.name, "/", sep="")
    filename = paste(projectiontime, '.tif')
    writeRaster(out.model, paste(model.dir, projectiontime, sep="/"), format="GTiff", options="COMPRESS=LZW")
}

######################################################################################
# model projection
######################################################################################

# function to check that the environmental layers used to project the  model are the same as the ones used
# 	to create the model object 
checkModelLayers = function(model.obj) {

	message("Checking environmental layers used for projection")
	# get the names of the environmental layers from the original model
	if (inherits(model.obj, "DistModel")) { # dismo package
		model.layers = colnames(model.obj@presence)
	} else if (inherits(model.obj, "gbm")) { # brt package
		model.layers = summary(model.obj)$var
	} else if (inherits(model.obj, "BIOMOD.models.out")) { # biomod package
		model.layers = model.obj@expl.var.names
	}
	
	# get the names of the climate scenario's env layers
	pred.layers = names(future.climate.scenario)
	
	# check if the env layers were in the original model
    if(sum(!(pred.layers %in% model.layers)) > 0 ){
		message("Dropping environmental layers not used in the original model creation...")
		# create a new list of env predictors by dropping layers not in the original model
		new.predictors = future.climate.scenario
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

######################################################################################
# model accuracy
######################################################################################

# function to save evaluate output
saveModelEvaluation = function(out.model, out.biomod.model, model.name) {
    model.dir = paste(wd, "/output_", model.name, sep="")
    save(out.model, file=paste(model.dir, "/dismo.eval.object.RData", sep=''))	# save the 'dismo::ModelEvalution' object

    # save all the model accuracy statistics provided in both dismo and biomod2
    rownames(out.biomod.model) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
    write.csv(t(round(out.biomod.model, digits=3)), file=paste(model.dir, "/combined.modelEvaluation.csv", sep=""))
    # EMG no guarantee these value are correct

    # save AUROC curve
    png(file=paste(model.dir, "/AUC.png", sep=''));
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
    model.dir = paste(wd, "/output_", model.name, sep="")

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
			range.values = seq(min(model.values[,j], na.rm=TRUE), max(model.values[,j], na.rm=TRUE), length.out=100)
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
			png(file=paste(model.dir, "/", save.name, "_response.png", sep=""))
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
    model.dir = paste(wd, "/output_", model.name, sep="")
	
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
		write.csv(varimpt.out, file=paste(model.dir, "/biomod2_like_VariableImportance.csv", sep=""))
	} else {
		write(paste(species, ": Cannot calculate variable importance for ", model.name, "object", sep=" "), stdout())
	}
}

# function to calculate variable importance values for dismo models based on Maxent's decrease in AUC 
# i.e., hold all but one predictor variable to its original values, resample that one predictor and recalculate model AUC
calculatePermutationVarImpt = function(out.model, model.eval, model.name) {
    model.dir = paste(wd, "/output_", model.name, sep="")
	
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

		# check for and remove any NA's present in the data
		no.na.sample.p = na.omit(sample.p); no.na.sample.a = na.omit(sample.a)
		if (nrow(no.na.sample.p) != nrow(sample.p)) {
			write(paste("calculatePermutationVarImpt(): NA's were removed from presence data!"), stdout())
		}
		if (nrow(no.na.sample.a) != nrow(sample.a)) {
			write(paste("calculatePermutationVarImpt(): NA's were removed from absence data!"), stdout())
		}
		
		# for each predictor variable 
		for (v in 1:length(env.vars)) {
		
			# resample from that variables' values, keeping other variable values the same 
			no.na.sample.p[,v] = sample(x=no.na.sample.p[,v], replace=FALSE)
			no.na.sample.a[,v] = sample(x=no.na.sample.a[,v], replace=FALSE)
			
			# re-evaluate model with sampled env values and NA omitted p/a
			if (model.name == "brt") {
				sample.eval = evaluate(p=no.na.sample.p, a=no.na.sample.a, model=out.model, n.trees=out.model$gbm.call$best.trees)				
			} else {
				sample.eval = evaluate(p=no.na.sample.p, a=no.na.sample.a, model=out.model)
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
		write.csv(permvarimpt.out, file=paste(model.dir, "/maxent_like_VariableImportance.csv", sep=""))
	} else {
		write(paste(species, ": Cannot calculate maxent-like variable importance for ", model.name, "object", sep=" "), stdout())
	}
}

# function to create HTML file with accuracy measures
generateHTML = function(sp.name, outputdir) {

	setwd(outputdir)

	# read in model outputs
	auccurve = readPNG(paste(outputdir, "/AUC.png", sep=""))
	accuracystats = read.csv(paste(outputdir, "/combined.modelEvaluation.csv", sep=""),	row.names=c(1))

	# create the output file 
	target <- HTMLInitFile(outdir=outputdir, filename=paste(sp.name,"_output", sep=""), BackGroundColor="#CCCCCC")

	# add content
	HTML(paste("<center><br><H1>Model Output for ", sp.name, sep=""), file=target)

	HTML("<br><H2>AUC:ROC curve", file=target)
	HTMLInsertGraph("AUC.png", file=target)

	HTML("<br><H2>Accuracy measures",file=target)
	HTML(accuracystats, file=target)

	# close the file
	HTMLEndFile()
}


# function to save evaluate output for BIOMOD2 models
saveBIOMODModelEvaluation = function(loaded.name, biomod.model) {
	# get and save the model evaluation statistics
	# EMG these must specified during model creation with the arg "models.eval.meth"
	evaluation = getModelsEvaluations(biomod.model)
	write.csv(evaluation, file=paste(getwd(), "/biomod2.modelEvaluation.txt", sep=""))

	# get the model predictions and observed values
	predictions = getModelsPrediction(biomod.model); obs = getModelsInputData(biomod.model, "resp.var");

	# get the model accuracy statistics using a modified version of biomod2's Evaluate.models.R
	combined.eval = sapply(model.accuracy, function(x){
		return(my.Find.Optim.Stat(Stat = x, Fit = predictions, Obs = obs))
	})
	# save all the model accuracy statistics provided in both dismo and biomod2
	rownames(combined.eval) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
	write.csv(t(round(combined.eval, digits=3)), file=paste(getwd(), "/combined.modelEvaluation.csv", sep=""))
		
	# save AUC curve
	require(pROC, quietly=T)
    roc1 <- roc(as.numeric(obs), as.numeric(predictions), percent=T)
	png(file=paste(getwd(), "/pROC.png", sep=''))
	plot(roc1, main=paste("AUC=",round(auc(roc1)/100,3),sep=""), legacy.axes=TRUE)
	dev.off()
	
	# get and save the variable importance estimates
	variableImpt = getModelsVarImport(biomod.model)
	if (!is.na(variableImpt)) {
	#EMG Note this will throw a warning message if variables (array) are returned	
		write.csv(variableImpt, file=paste(getwd(), "/variableImportance.txt", sep=""))
	} else {
		message("VarImport argument not specified during model creation!")
		#EMG must create the model with the arg "VarImport" != 0
	}

	# save response curves (Elith et al 2005)
	png(file=paste(getwd(), "/mean_response_curves.png", sep=''))
		test=response.plot2(models = loaded.name, Data = getModelsInputData(biomod.model,"expl.var"),
			show.variables = getModelsInputData(biomod.model,"expl.var.names"), fixed.var.metric = "mean") 
			#, data_species = getModelsInputData(biomod.model,"resp.var"))
			# EMG need to investigate why you would want to use this option - uses presence data only
	dev.off()
}