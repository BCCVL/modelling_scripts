src/org/bccvl/compute/rscripts/gbm.R# set CRAN mirror in case we need to download something
# TODO: this should be done on demand or on user basis...
r <- getOption("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos=r)
# TODO: alse creating and populating add on package location is something that should not be done system wide

#script to run to develop distribution models
###check if libraries are installed, install if necessary and then load them
necessary=c("rgdal", "SDMTools", "biomod2", "R2HTML", "png") #list the libraries needed
installed = necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >=1) {
    install.packages(necessary[!installed], dep = T) #if library is not installed, install it
}
for (lib in necessary) {
    library(lib,character.only=T) #load the libraries
}

###read in the necessary observation, background and environmental data
populate.data = TRUE #variable to define if there is a need to generate occur & background environmental info
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
}

current.climate.scenario = stack(enviro.data.current)
if (is.null(enviro.data.future)) {
    project.bioclim=FALSE
} else {
    future.climate.scenario = stack(enviro.data.future)
}

# source helper functions (err.null, getModelObject, checkModelLayers, saveModelProject)
source(paste(function.path, "/my.Helper.Functions.R", sep=""))

# source modified projeciton function to create .tif
source(paste(function.path, "/my.BIOMOD_Projection.R", sep=""))

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
	outdir = paste(wd,'/output_gbm',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(GBM = gbm.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.gbm <- BIOMOD_Modeling(data=myBiomodData, models=c('GBM'),	models.options=myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth=biomod.models.eval.meth, SaveObj=TRUE,
		rescal.all.models=biomod.rescal.all.models, do.full.models=biomod.do.full.models, 
		modeling.id=biomod.modeling.id)
	if (!is.null(myBiomodModelOut.gbm)) {		
		save(myBiomodModelOut.gbm, file=paste(outdir,"/model.object.RData",sep='')) #save out the model object
		# predict for current climate scenario
		gbm.proj.c = my.BIOMOD_Projection(modeling.output=myBiomodModelOut.gbm, new.env=current.climate.scenario, proj.name="current", 
			xy.new.env=biomod.xy.new.env, selected.models=biomod.selected.models, binary.meth=biomod.binary.meth, 
			filtered.meth= biomod.filtered.meth, compress=biomod.compress, 
			build.clamping.mask=biomod.build.clamping.mask, silent=opt.biomod.silent, do.stack=opt.biomod.do.stack, 
			keep.in.memory=opt.biomod.keep.in.memory, output.format=opt.biomod.output.format)
		# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
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

if (project.gbm) {
	gbm.obj = getModelObject("gbm") # get the model object
	outdir = paste(wd,'/output_gbm',sep=''); setwd(outdir) #set the working directory
	if (!is.null(gbm.obj)) {
		predictors = checkModelLayers(gbm.obj)
		gbm.proj = my.BIOMOD_Projection(modeling.output=gbm.obj, new.env=predictors, proj.name="future", 
			xy.new.env=biomod.xy.new.env, selected.models=biomod.selected.models, binary.meth=biomod.binary.meth, 
			filtered.meth=biomod.filtered.meth, compress=biomod.compress, 
			build.clamping.mask=biomod.build.clamping.mask, silent=opt.biomod.silent, do.stack=opt.biomod.do.stack, 
			keep.in.memory=opt.biomod.keep.in.memory, output.format=opt.biomod.output.format)
		# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
		rm(list=c("gbm.obj", "gbm.proj")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load gbm.obj from", wd, "output_gbm", sep=": "), stdout())
	}
}


######################################################################################
# evaluate
######################################################################################

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy = c(dismo.eval.method, biomod.models.eval.meth)

if (evaluate.gbm) {
	gbm.obj = getModelObject("gbm") # get the model object
	if (!is.null(gbm.obj)) {
		gbm.loaded.model = BIOMOD_LoadModels(gbm.obj, models="GBM") # load model
		saveBIOMODModelEvaluation(gbm.loaded.model, gbm.obj, "gbm") 	# save output
		rm("gbm.obj") #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load gbm.obj from", getwd(), sep=": "), stdout())
	}
}