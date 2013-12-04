src/org/bccvl/compute/rscripts/gam.R# set CRAN mirror in case we need to download something
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
	outdir = paste(wd,'/output_gam',sep=''); dir.create(outdir,recursive=TRUE); #create the output directory
	setwd(outdir) # set the working directory (where model results will be stored)
	myBiomodData = formatBiomodData() # 1. Format the data
	myBiomodOptions <- BIOMOD_ModelingOptions(GAM = gam.BiomodOptions) # 2. Define the model options
	# 3. Compute the model
	myBiomodModelOut.gam <- BIOMOD_Modeling(data=myBiomodData, models=c('GAM'),	models.options=myBiomodOptions,
		NbRunEval=biomod.NbRunEval,	DataSplit=biomod.DataSplit,	Yweights=biomod.Yweights, Prevalence=biomod.Prevalence,
		VarImport=biomod.VarImport,	models.eval.meth=biomod.models.eval.meth, SaveObj=TRUE,
		rescal.all.models=biomod.rescal.all.models, do.full.models=biomod.do.full.models, 
		modeling.id=biomod.modeling.id)
	if (!is.null(myBiomodModelOut.gam)) {		
		save(myBiomodModelOut.gam, file=paste(outdir,"/model.object.RData",sep='')) #save out the model object
		# predict for current climate scenario
		gam.proj.c = my.BIOMOD_Projection(modeling.output=myBiomodModelOut.gam, new.env=current.climate.scenario, proj.name="current", 
			xy.new.env=biomod.xy.new.env, selected.models=biomod.selected.models, binary.meth=biomod.binary.meth, 
			filtered.meth=biomod.filtered.meth, compress=biomod.compress, 
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

if (project.gam) {
	gam.obj = getModelObject("gam") # get the model object
	outdir = paste(wd,'/output_gam',sep=''); setwd(outdir) #set the working directory
	if (!is.null(gam.obj)) {
		predictors = checkModelLayers(gam.obj)
		gam.proj = my.BIOMOD_Projection(modeling.output=gam.obj, new.env=predictors, proj.name="future", 
			xy.new.env=biomod.xy.new.env, selected.models=biomod.selected.models, binary.meth=biomod.binary.meth, 
			filtered.meth=biomod.filtered.meth, compress=biomod.compress, 
			build.clamping.mask=biomod.build.clamping.mask, silent=opt.biomod.silent, do.stack=opt.biomod.do.stack, 
			keep.in.memory=opt.biomod.keep.in.memory, output.format=opt.biomod.output.format)
		# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
		rm(list=c("gam.obj", "gam.proj")) #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load gam.obj from", wd, "output_gam", sep=": "), stdout())
	}
}


######################################################################################
# evaluate
######################################################################################

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy = c(dismo.eval.method, biomod.models.eval.meth)

if (evaluate.gam) {
	gam.obj = getModelObject("gam") # get the model object
	if (!is.null(gam.obj)) {
		gam.loaded.model = BIOMOD_LoadModels(gam.obj, models="GAM") # load model
		saveBIOMODModelEvaluation(gam.loaded.model, gam.obj, "gam") 	# save output
		rm("gam.obj") #clean up the memory
	} else {
		write(paste("FAIL!", species, "Cannot load gam.obj from", getwd(), sep=": "), stdout())
	}
}