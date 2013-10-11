src/org/bccvl/compute/rscripts/ann.R# set CRAN mirror in case we need to download something
# TODO: this should be done on demand or on user basis...
r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# TODO: alse creating and populating add on package location is something that should not be done system wide

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("biomod2","SDMTools", "rgdal", "pROC") #list the libraries needed
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
    project.bioclim=FALSE
} else {
    future.climate.scenario = stack(enviro.data.future)
}
# get the name of the scenario for use during projection
es.name=basename(enviro.data.current)

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
	myBiomodData <- BIOMOD_FormatingData(resp.var = biomod.data.pa, expl.var = stack(current.climate.scenario),	
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
		# predict for current climate scenario
		ann.proj.c = BIOMOD_Projection(modeling.output=myBiomodModelOut.ann, new.env=current.climate.scenario, proj.name=es.name, 
			xy.new.env = biomod.xy.new.env,	selected.models = biomod.selected.models, binary.meth = biomod.binary.meth, 
			filtered.meth = biomod.filtered.meth, compress = biomod.compress, 
			build.clamping.mask = biomod.build.clamping.mask, silent = opt.biomod.silent, do.stack = opt.biomod.do.stack, 
			keep.in.memory = opt.biomod.keep.in.memory,	output.format = opt.biomod.output.format)
		# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
	}
}

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

############### BIOMOD2 Models ###############
###############
#
# BIOMOD_Projection(modeling.output, new.env, proj.name, xy.new.env = NULL, selected.models = 'all', 
#	binary.meth = NULL, filtered.meth = NULL, compress = 'xz', build.clamping.mask = TRUE, ...)
#
# modeling.output "BIOMOD.models.out" object produced by a BIOMOD_Modeling run
# new.env A set of explanatory variables onto which models will be projected. It could be a data.frame, a matrix, 
#	or a rasterStack object. Make sure the column names (data.frame or matrix) or layerNames (rasterStack) perfectly 
#	match with the names of variables used to build the models in the previous steps.
# proj.name	a character defining the projection name (a new folder will be created with this name)
# xy.new.env optional coordinates of new.env data. Ignored if new.env is a rasterStack
# selected.models 'all' when all models have to be used to render projections or a subset vector of modeling.output 
#	models computed (accessing with the slot @models.computed of your "BIOMOD.models.out" object)
#	EMG: eg, myBiomodModelOut@models.computed[grep("_RF", getModelsBuiltModels(myBiomodModelOut))]
# binary.meth a vector of a subset of models evaluation method computed before (see BIOMOD_Modeling). If NULL then 
#	no binary transformation computed, else the given binary techniques will be used to transform the projection 
#	into 0/1 data.
#	EMG: dimnames(getModelsEvaluations(myBiomodModelOut)[[1]] to see which methods were used
# filtered.meth	a vector of a subset of models evaluation method computed before (see BIOMOD_Modeling). If NULL then 
#	no filtering transformation computed, else the given binary techniques will be used to transform the projection by 
#	settting to 0 the probability values below a specific threshold.
# compress compression format of objects stored on your hard drive. May be one of ‘xz’, ‘gzip’ or NULL
# build.clamping.mask if TRUE, a clamping mask that identifies locations where predictions are uncertain because the
#	values of the variables are outside the range used for calibrating the models will be saved
#	EMG: It appears that setting this to FALSE will still generate a mask file
# ... Additional arguments:
# silent logical, if TRUE, console outputs are turned off
# do.stack logical, if TRUE, attempt to save all projections in a unique object i.e RasterStack. If FALSE or if objects are
#	too heavy to be load all together in memory, projections will be stored into separated files
# keep.in.memory logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object. 
#	That can be usefull to prevent memory issues.
# output.format whether ‘.Rdata’, ‘.grd’ or ‘.img’ defining projections saving format (on hard drive). If new.env argument 
#	is under table format (data.frame or matrix), the only choice you have is ‘.Rdata’
#	EMG: If .grd is selected (meta-data), also get a .gri file (values)
#
###############

if (project.ann) {	
	ann.obj = getModelObject("ann") # get the model object
	outdir = paste(wd,'output_ann/',sep=''); setwd(outdir) #set the working directory
	if (!is.null(ann.obj)) {
		predictors = checkModelLayers(ann.obj)
		ann.proj = BIOMOD_Projection(modeling.output=ann.obj, new.env=predictors, proj.name=es.name, 
			xy.new.env = biomod.xy.new.env,	selected.models = biomod.selected.models, binary.meth = biomod.binary.meth, 
			filtered.meth = biomod.filtered.meth, compress = biomod.compress, 
			build.clamping.mask = biomod.build.clamping.mask, silent = opt.biomod.silent, do.stack = opt.biomod.do.stack, 
			keep.in.memory = opt.biomod.keep.in.memory,	output.format = opt.biomod.output.format)
		# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
		rm(list=c("ann.obj", "ann.proj")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load ann.obj from", wd, "output_ann", sep=": "), stdout())
	}
}

######################################################################################
# model accuracy helpers
######################################################################################

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy = c(dismo.eval.method, biomod.models.eval.meth)

# function to save evaluate output
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

if (evaluate.ann) {	
	ann.obj = getModelObject("ann") # get the model object
	if (!is.null(ann.obj)) {
		ann.loaded.model = BIOMOD_LoadModels(ann.obj, models="ANN") # load model
		saveBIOMODModelEvaluation(ann.loaded.model, ann.obj) 	# save output
		rm("ann.obj") #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load ann.obj from", getwd(), sep=": "), stdout())
	}
}