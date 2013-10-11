#initial arguments used to define the inputs, models, outputs, etc

# read in the arguments listed at the command line
args=(commandArgs(TRUE))  
# check to see if arguments are passed
if(length(args)==0){
    print("No arguments supplied.")
    # leave all args as default values

wd = "/home/jc140298/bccvl/ABT/" #define the working directory - where to put the outputs
species = "ABT"	#define the species of interest
occur.data = "/home/jc140298/bccvl/ABT/occur.csv" #define the lon/lat of the observation records -- 2 column matrix of longitude and latitude
bkgd.data = "/home/jc140298/bccvl/ABT/bkgd.csv" #define the the lon/lat of the background / psuedo absence points to use -- 2 column matrix of longitude and latitude

} else {
	for(i in 1:length(args)) { 
		eval(parse(text=args[[i]])) 
	}
	# expecting wd, species, occur.data, and bkgd.data
}
# EMG need to expand this to include all other args or come up with a way to parse this properly

enviro.data = c("/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_01.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_04.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_05.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_06.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_12.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_15.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_16.asc",
"/home/jc165798/working/BCCVL/envirodata/climate_1990/bioclim_17.asc") #define the enviro data to use -- assumed location of data files in ascii grid format
enviro.data.names = c("bioclim_01","bioclim_04","bioclim_05","bioclim_06",
"bioclim_12","bioclim_15","bioclim_16","bioclim_17") #define the names of the enviro data
enviro.data.type = c("continuous","continuous","continuous","continuous",
"continuous","continuous","continuous","continuous") #type in terms of continuous or categorical

### define the models to be used
model.bioclim = FALSE #boolean to run BIOCLIM algorithm -- all envirodata must be continuous
model.domain = FALSE #boolean to run DOMAIN algorithm -- all envirodata must be continuous
model.mahal = FALSE #boolean to run MAHALANOBIS algorithm -- all envirodata must be continuous
model.geodist = FALSE #boolean to run geographic distances algorithm -- only requires lon/lat of observations
model.convHull = FALSE #boolean to run convex hulls algorithm -- only requires lon/lat of observations
model.circles = FALSE #boolean to run "circles" algorithm -- only requires lon/lat of observations
model.geoIDW = FALSE #boolean to run inverse distance weighted algorithm -- only requires lon/lat of observations and pseudo-absences
model.voronoiHull = FALSE #boolean to run Voronoi Hulls algorithm -- only requires lon/lat of observations and pseudo-absences

# need to determine the number of sites for default brt.site.weights arg
model.brt = TRUE #boolean to run Boosted regression tree algorithm
if (model.brt) { #additional parameters to set
	brt.fold.vector = NULL #a fold vector to be read in for cross validation with offsets
	brt.tree.complexity = 7 #sets the complexity of individual trees
	brt.learning.rate = 0.001 #sets the weight applied to individual trees
	brt.bag.fraction = 0.5 #sets the proportion of observations used in selecting variables
	brt.site.weights = NULL #allows varying weighting for sites; rep(1, nrow(data))
	brt.var.monotone = NULL #restricts responses to individual predictors to monotone; rep(0, length(enviro.data))
	brt.n.folds = 5 #number of folds
	brt.prev.stratify = TRUE #prevalence stratify the folds - only for presence/absence data
	brt.family = "bernoulli" #family - bernoulli (=binomial), poisson, laplace or gaussian
	brt.n.trees = 2000 #number of initial trees to fit
	brt.step.size = brt.n.trees #numbers of trees to add at each cycle
	brt.max.trees = 10000 #max number of trees to fit before stopping
	brt.tolerance.method = "auto" #method to use in deciding to stop - "fixed" or "auto"
	brt.tolerance = 0.001 #tolerance value to use - if method == fixed is absolute, if auto is multiplier * total mean deviance
	brt.keep.data = FALSE #Logical. keep raw data in final model
	brt.plot.main = FALSE #Logical. plot hold-out deviance curve
	brt.plot.folds = FALSE #Logical. plot the individual folds as well
	brt.verbose = FALSE #Logical. control amount of screen reporting
	brt.silent = FALSE #Logical. to allow running with no output for simplifying model)
	brt.keep.fold.models = FALSE  #Logical. keep the fold models from cross valiation
	brt.keep.fold.vector = FALSE #Logical. allows the vector defining fold membership to be kept
	brt.keep.fold.fit = FALSE #Logical. allows the predicted values for observations from cross-validation to be kept
}

model.maxent = FALSE #boolean to run maxent algorithm
if (model.maxent) {
	maxent.jar = "/home/jc165798/working/BCCVL/maxent.jar" #define location of maxent.jar file
	outputformat="logistic" #options include logistic, cumulative, raw
	randomseed=FALSE
	logscale=TRUE
	removeduplicates=TRUE
	randomtestpoints=0
	betamultiplier=1
	maximumbackground=10000
	biasfile=NULL
	testsamplesfile=NULL
	replicates=10
	replicatetype="crossvalidate" #options include crossvalidate, bootstrap, subsample
	linear=TRUE
	quadratic=TRUE
	product=TRUE
	threshold=TRUE
	hinge=TRUE
	addsamplestobackground=TRUE
	addallsamplestobackground=FALSE
	fadebyclamping=FALSE
	extrapolate=TRUE
	autofeature=TRUE
	doclamp=TRUE
	maximumiterations=500
	convergencethreshold=1.00E-05
	lq2lqptthreshold=80
	l2lqthreshold=10
	hingethreshold=15
	beta_threshold=-1
	beta_categorical=-1
	beta_lqp=-1
	beta_hinge=-1
	defaultprevalence=0.5
	nodata=-9999
}

############### BIOMOD2 Models ###############
#
# general arguments to perform any biomod modelling
#
biomod.NbRunEval = 10  # n-fold cross-validation; ignored if DataSplitTable is filled
biomod.DataSplit = 100 ## % for calibrating/training, remainder for testing; ignored if DataSplitTable is filled
biomod.Yweights = NULL #response points weights
biomod.Prevalence = NULL #either NULL (default) or a 0-1 numeric used to build "weighted response weights"
biomod.VarImport = 3 #number of resampling of each explanatory variable to measure the relative importance of each variable for each selected model
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS") #vector of evaluation metrics 
biomod.rescal.all.models = TRUE #if true, all model prediction will be scaled with a binomial GLM
biomod.do.full.models = TRUE #if true, models calibrated and evaluated with the whole dataset are done; ignored if DataSplitTable is filled
biomod.modeling.id = species #character, the ID (=name) of modeling procedure. A random number by default
# biomod.DataSplitTable = NULL #a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a "RUN". If filled, args NbRunEval, DataSplit and do.full.models will be ignored
# EMG Need to test whether a NULL values counts as an argument

#
# model-specific arguments to create a biomod model
#

model.glm = FALSE #boolean to run generalized linear model algorithm
if (model.glm) {
	glm.myBiomodOptions <- list(
		type = "quadratic",	#"simple", "quadratic" or "polynomial"; switched off if myFormula is not NULL
		interaction.level = 0, #integer corresponding to the interaction level between variables considered; switched off if myFormula is not NULL
		myFormula = NULL, #specific formula; if not NULL, type and interaction.level are args are switched off
		test = "AIC", #"AIC", "BIC" or "none"
		family = "binomial", #"binomial", "gaussian", "gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"
		mustart = 0.5, #starting values for the vector of means
		control = list(
			epsilon = 1e-08, #positive convergence tolerance e
			maxit = 50, #integer giving the maximal number of IWLS iterations
			trace = FALSE #logical indicating if output should be produced for each iteration
		)
	)
}

model.gbm = FALSE #boolean to run generalized boosting model algorithm
if (model.gbm) {
	my.distribution = "bernoulli" #"bernoulli", "gaussian", "laplace", "tdist", "huberized", "multinomial", "adaboost", "poisson", "coxph", "quantile", or "pairwise"
	if (my.distribution == "quantile") {
		my.distribution = list(name="quantile", alpha=0.25) #alpha is the quantile to estimate
	} else if (my.distribution == "tdist") {
		my.distribution = list(name="tdist", df=4) #df is degrees of freedom
	} else if (my.distribution == "pairwise") {
		my.distribution = list(name="pairwise", 
			group="", #group is character vector with the column names of data that jointly indicate the group an instance belongs to
			metric="", #optional; "conc", "mrr", "map", or "ndcg"
			max.rank="" #optional; for ndcg and mrr, a cut-off can be chosen using a positive integer parameter
		)
		# EMG Need to test whether an empty string "" breaks things
	}		
	gbm.BiomodOptions <- list(distribution = my.distribution,
		interaction.depth = 7, #maximum depth of variable interactions
		shrinkage = 0.001, #a shrinkage parameter applied to each tree in the expansion
		bag.fraction = 0.5, #the fraction of the training set observations randomly selected to propose the next tree in the expansion
		train.fraction = 1, #The first train.fraction * nrows(data) observations are used to fit the gbm
		n.trees = 500, #the total number of trees to fit
		cv.folds = 5 #number of cross-validation folds to perform
	)
}

model.gam = FALSE #boolean to run generalized additive model algorithm
if (model.gam) {
	my.algo = "GAM_mgcv" #"GAM_mgcv", "GAM_gam" or "BAM_mgcv"
	if (my.algo == "GAM_gam") {	
		my.control = list(
			epsilon=1e-07, #convergence threshold for local scoring iterations
			bf.epsilon = 1e-07, #convergence threshold for backfitting iterations
			maxit=30, #maximum number of local scoring iterations
			bf.maxit = 30, #maximum number of backfitting iterations
			trace=FALSE #should iteration details be printed while gam is fitting the model
		) 
	} else if (my.algo == "GAM_mgcv") {
		my.control = list(
			irls.reg=0.0, #for most models this should be 0; the size of the ridge regression penalty to the model to impose identifiability
			epsilon = 1e-06, #this is used for judging conversion of the GLM IRLS loop
			maxit = 100, #maximum number of IRLS iterations to perform
			mgcv.tol=1e-7, #the convergence tolerance parameter to use in GCV/UBRE optimization
			mgcv.half=15, #if a step of the GCV/UBRE optimization method leads to a worse GCV/UBRE score, then the step length is halved; this is the number of halvings to try before giving up
			trace = FALSE, #set this to TRUE to turn on diagnostic output
			rank.tol=.Machine$double.eps^0.5, #the tolerance used to estimate the rank of the fitting problem
			nlm=list(), #list of control parameters to pass to nlm if this is used for outer estimation of smoothing parameters (not default)
			optim=list(), #list of control parameters to pass to optim if this is used for outer estimation of smoothing parameters (not default)
			newton=list(), #list of control parameters to pass to default Newton optimizer used for outer estimation of log smoothing parameters
			outerPIsteps=0, #the number of performance interation steps used to initialize outer iteration
			idLinksBases=TRUE, #if smooth terms have their smoothing parameters linked via the id mechanism (see s), should they also have the same bases. Set this to FALSE only if you are sure you know what you are doing
			scalePenalty=TRUE, #this option rescales the penalty matrices to accomodate this problem. Probably should be set to FALSE if you are linking smoothing parameters but have set idLinkBases to FALSE
			keepData=FALSE #should a copy of the original data argument be kept in the gam object
		)
	} else {
		my.control = list(epsilon = 1e-06, trace = FALSE, maxit = 100)
	}
	gam.BiomodOptions <- list(algo=my.algo, 
		type = "s_smoother", #the smoother used to generate the formula; only "s_smoother" available at time; switched off if myFormula is not NULL
		k = NULL, #a smooth term in a formula argument to gam (see gam s or mgcv s)
		interaction.level = 0, #integer corresponding to the interaction level between variables considered; switched off if myFormula is not NULL
		myFormula = NULL, #specific formula; if not NULL, type and interaction.level are args are switched off
		family = "binomial", #"bernoulli", "gaussian", "laplace", "tdist", "huberized", "multinomial", "adaboost", "poisson", "coxph", "quantile", or "pairwise"
		control=my.control 		
	)
}

model.cta = FALSE #boolean to run classification tree analysis algorithm
if (model.cta) {
	cta.BiomodOptions <- list(
		method = "class", #"anova", "poisson", "class" or "exp"
		parms = "default", #optional parameters for the splitting function
		cost = NULL, #a vector of non-negative costs, one for each variable in the model. Defaults to one for all variables
		control = list(
			xval = 5, #number of cross-validations
			minbucket = 5, #the minimum number of observations in any terminal <leaf> node
			minsplit = 5, #the minimum number of observations that must exist in a node in order for a split to be attempted
			cp = 0.001, #complexity parameter
			maxdepth = 25 #Set the maximum depth of any node of the final tree, with the root node counted as depth 0. Values greater than 30 rpart will give nonsense results on 32-bit machines
		)
	)
}

model.ann = FALSE #boolean to run artificial neural network algorithm
if (model.ann) {
	ann.BiomodOptions <- list(
		NbCV = 5, #nb of cross validation to find best size and decay parameters
		rang = 0.1, #Initial random weights on [-rang, rang]
		maxit = 200 #maximum number of iterations. Default 100
	)
}

model.sre = FALSE #boolean to run surface range envelop algorithm			
if (model.sre) {
	sre.BiomodOptions <- list(
		quant = 0.025 #quantile of "extreme environmental variable" removed for selection of species envelops
	)
}

model.fda = FALSE #boolean to run flexible discriminant analysis algorithm
if (model.fda) {
	fda.BiomodOptions <- list(
		method = "mars" #regression method used in optimal scaling; "polyreg", "mars", "bruto" or "gen.ridge"
	)
}

model.mars = FALSE #boolean to run multiple adaptive regression splines algorithm
if (model.mars) {
	mars.BiomodOptions <- list(
		degree = 2, #an optional integer specifying maximum interaction degree (default is 1)
		penalty = 2, #an optional value specifying the cost per degree of freedom charge (default is 2)
		thresh = 0.001, #an optional value specifying forward stepwise stopping threshold (default is 0.001)
		prune = TRUE #an optional logical value specifying whether the model should be pruned in a backward stepwise fashion (default is TRUE)
	)
}

model.rf = FALSE #boolean to run random forest algorithm
 if (model.rf) {
 	rf.BiomodOptions <- list(
		do.classif = TRUE, #if TRUE classification random.forest computed else regression random.forest will be done
		ntree = 50, #number of trees to grow
		mtry = "default" #number of variables randomly sampled as candidates at each split
	)
}

model.biomod.maxent = FALSE #boolean to run {biomod} maxent algorithm
if (model.biomod.maxent) {
	biomod.maxent.BiomodOptions <- list(
		path_to_maxent.jar = "/home/jc165798/working/BCCVL/maxent.jar", #character, the link to maxent.jar file (the working directory by default)
		maximumiterations = 200, #integer (default 200), maximum iteration done
		visible = FALSE, #logical (default FALSE), make the Maxent user interface visible
		linear = TRUE, #logical (default TRUE), allow linear features to be used
		quadratic = TRUE, #logical (default TRUE), allow quadratic features to be used
		product = TRUE, #logical (default TRUE), allow product features to be used
		threshold = TRUE, #logical (default TRUE), allow threshold features to be used
		hinge = TRUE, #logical (default TRUE), allow hinge features to be used
		lq2lqptthreshold = 80, #integer (default 80), number of samples at which product and threshold features start being used
		l2lqthreshold = 10, #integer (default 10), number of samples at which quadratic features start being used
		hingethreshold = 15, #integer (default 15), number of samples at which hinge features start being used
		beta_threshold = -1, #numeric (default -1.0), regularization parameter to be applied to all threshold features; negative value enables automatic setting
		beta_categorical = -1, #numeric (default -1.0), regularization parameter to be applied to all categorical features; negative value enables automatic setting
		beta_lqp = -1, #numeric (default -1.0), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
		beta_hinge = -1, #numeric (default -1.0), regularization parameter to be applied to all hinge features; negative value enables automatic setting
		defaultprevalence = 0.5 #numeric (default 0.5), default prevalence of the species: probability of presence at ordinary occurrence points
	)
}

# save workspace to set arguments used by 01.model.current.R
save.image(paste(wd, "/01.init.args.model.current.", species, ".RData", sep=""))
save.image(paste(wd, "/01.init.args.model.current.", species, ".Rascii", sep=""), ascii=TRUE) # for Daniel