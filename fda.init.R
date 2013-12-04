# This is python string template which is currently tightly tied to
# ../fda.py
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

# define path to helper functions
function.path = "/home/jc140298/modelling_scripts"

model.fda = TRUE #boolean to run flexible discriminant analysis algorithm
project.fda = TRUE #boolean to project flexible discriminant analysis algorithm
evaluate.fda = TRUE #boolean to evaluate flexible discriminant analysis algorithm

############### BIOMOD2 Models ###############
#
# general parameters to perform any biomod modelling
#
biomod.NbRunEval = {{ Nb_Run_Eval }}  # default 10; n-fold cross-validation; ignored if DataSplitTable is filled
biomod.DataSplit = {{ Data_Split }} # default 100; % for calibrating/training, remainder for testing; ignored if DataSplitTable is filled
biomod.Yweights = NULL #response points weights
biomod.Prevalence = NULL #either NULL (default) or a 0-1 numeric used to build "weighted response weights"
# EMG Not currently implemented through the BCCVL
biomod.VarImport = {{ Var_Import }} # default 0; number of resampling of each explanatory variable to measure the relative importance of each variable for each selected model
# EMG this parameter needs to be specified in order to get VariableImportance metrics during model evaluation
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS") #vector of evaluation metrics 
biomod.rescal.all.models = FALSE #if true, all model prediction will be scaled with a binomial GLM
# EMG need to have the default 0-1000 scaling for visualizing the output maps
biomod.do.full.models = TRUE #if true, models calibrated and evaluated with the whole dataset are done; ignored if DataSplitTable is filled
biomod.modeling.id = "{{ species }}" #character, the ID (=name) of modeling procedure. A random number by default
# biomod.DataSplitTable = NULL #a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a "RUN". If filled, args NbRunEval, DataSplit and do.full.models will be ignored
# EMG Need to test whether a NULL values counts as an argument

# model-specific arguments to create a biomod model
fda.BiomodOptions <- list(
	method = "mars" #regression method used in optimal scaling; "polyreg", "mars", "bruto" or "gen.ridge"
)	

############### BIOMOD2 Models ###############
#
# general parameters to project any biomod modelling
#
#modeling.output #"BIOMOD.models.out" object produced by a BIOMOD_Modeling run
#new.env #a set of explanatory variables onto which models will be projected; must match variable names used to build the models
#proj.name #a character defining the projection name (a new folder will be created with this name)
biomod.xy.new.env = NULL #optional coordinates of new.env data. Ignored if new.env is a rasterStack
# EMG always using a rasterStack
biomod.selected.models = "all" #'all' when all models have to be used to render projections or a subset vector of modeling.output models computed (eg, = grep(’_RF’, getModelsBuiltModels(myBiomodModelOut)))
# EMG If running one model at a time, this parameter becomes irrevelant
biomod.binary.meth = NULL #a vector of a subset of models evaluation method computed in model creation 
biomod.filtered.meth = NULL #a vector of a subset of models evaluation method computed in model creation
#biomod.compress = "{{ compress }}" # default 'xz'; compression format of objects stored on your hard drive. May be one of ‘xz’, ‘gzip’ or NULL
# EMG compress = NULL fails; currently commented out of glm.R
biomod.build.clamping.mask = TRUE #if TRUE, a clamping mask will be saved on hard drive
opt.biomod.silent = FALSE #logical, if TRUE, console outputs are turned off
opt.biomod.do.stack = TRUE #logical, if TRUE, attempt to save all projections in a unique object i.e RasterStack
opt.biomod.keep.in.memory = TRUE #logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object
opt.biomod.output.format = ".tif" #'.Rdata', '.grd' or '.img'; if NULL, and new.env is not a Raster class, output is .RData defining projections saving format (on hard drive)
# EMG 'tif' added for BCCVL, requires my.BIOMOD_Projection.R

# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS")