# This is python string template which is currently tightly tied to
# ../biomod.maxent.py
{% macro strvector(value) -%}
  {%- if value -%}
    c(
    {%- for item in value -%}
      {{ '"%s"'|format(item) -}}
      {{ "," if not loop.last }}
    {%- endfor -%}
    )
  {%- else -%}
    NULL
  {%- endif -%}
{%- endmacro %}

.libPaths("{{ rlibdir }}")
wd = "{{ workdir }}" #define the working directory
species = "{{ species }}" #define the species of interest
occur.data = "{{ occurrence }}" #define the lon/lat of the observation records -- 2 column matrix of longitude and latitude
bkgd.data = {{ '"%s"' % background if background else "NULL" }} 
#define the the lon/lat of the background / psuedo absence points to use -- 2 column matrix of longitude and latitude
enviro.data.names = {{ strvector(enviro['names']) }} #define the names of the enviro data
enviro.data.current = {{ strvector(enviro['data']) }} #define the current enviro data to use
enviro.data.type = {{ strvector(enviro['type']) }} #type in terms of continuous or categorical
enviro.data.future = {{ strvector(future['data']) }} #define the future enviro data to use

model.biomod.maxent = TRUE #boolean to run {biomod} maxent algorithm
project.biomod.maxent = TRUE #boolean to project {biomod} maxent algorithm
evaluate.biomod.maxent = TRUE #boolean to evaluate {biomod} maxent algorithm

############### BIOMOD2 Models ###############
#
# general parameters to perform any biomod modelling
#
biomod.NbRunEval = {{ Nb_Run_Eval }}  # default 10; n-fold cross-validation; ignored if DataSplitTable is filled
biomod.DataSplit = {{ Data_Split }} # default 100; % for calibrating/training, remainder for testing; ignored if DataSplitTable is filled
biomod.Yweights = NULL #response points weights
biomod.Prevalence = NULL #either NULL (default) or a 0-1 numeric used to build "weighted response weights"
biomod.VarImport = {{ Var_Import }} # default 0; number of resampling of each explanatory variable to measure the relative importance of each variable for each selected model
#EMG this parameter needs to be specified in order to get VariableImportance metrics during model evaluation
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS") #vector of evaluation metrics 
biomod.rescal.all.models = TRUE #if true, all model prediction will be scaled with a binomial GLM
biomod.do.full.models = TRUE #if true, models calibrated and evaluated with the whole dataset are done; ignored if DataSplitTable is filled
biomod.modeling.id = "{{ species }}" #character, the ID (=name) of modeling procedure. A random number by default
# biomod.DataSplitTable = NULL #a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a "RUN". If filled, args NbRunEval, DataSplit and do.full.models will be ignored
# EMG Need to test whether a NULL values counts as an argument

# model-specific arguments to create a biomod model
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

############### BIOMOD2 Models ###############
#
# general parameters to project any biomod modelling
#
#modeling.output #"BIOMOD.models.out" object produced by a BIOMOD_Modeling run
#new.env #a set of explanatory variables onto which models will be projected; must match variable names used to build the models
#proj.name #a character defining the projection name (a new folder will be created with this name)
biomod.xy.new.env = NULL #optional coordinates of new.env data. Ignored if new.env is a rasterStack
#biomod.selected.models = "{{ species }}"'all' #'all' when all models have to be used to render projections or a subset vector of modeling.output models computed (eg, = grep(’_RF’, getModelsBuiltModels(myBiomodModelOut)))
# EMG If running one model at a time, this parameter becomes irrevelant
biomod.binary.meth = NULL #a vector of a subset of models evaluation method computed in model creation 
biomod.filtered.meth = NULL #a vector of a subset of models evaluation method computed in model creation
biomod.compress = "{{ compress }}" # default 'gzip'; compression format of objects stored on your hard drive. May be one of ‘xz’, ‘gzip’ or NULL
biomod.build.clamping.mask = TRUE #if TRUE, a clamping mask will be saved on hard drive
opt.biomod.silent = FALSE #logical, if TRUE, console outputs are turned off
opt.biomod.do.stack = TRUE #logical, if TRUE, attempt to save all projections in a unique object i.e RasterStack
opt.biomod.keep.in.memory = TRUE #logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object
opt.biomod.output.format = NULL #'.Rdata', '.grd' or '.img'; if NULL, and new.env is not a Raster class, output is .RData defining projections saving format (on hard drive)


# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS")