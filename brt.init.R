# This is a Python string template file and is currently tightly tied
# to ../brt.py.
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

model.brt = TRUE #boolean to run Boosted regression tree algorithm
project.brt = TRUE #boolean to project Boosted regression tree algorithm
evaluate.brt = TRUE #boolean to evaluate Boosted regression tree algorithm

#additional parameters to set for modelling brt
brt.fold.vector = NULL #a fold vector to be read in for cross validation with offsets
brt.tree.complexity = {{ tree_complexity }} #sets the complexity of individual trees
brt.learning.rate = {{ learning_rate }} #sets the weight applied to individual trees
brt.bag.fraction = {{ bag_fraction }} #sets the proportion of observations used in selecting variables
brt.site.weights = NULL #allows varying weighting for sites; rep(1, nrow(data))
brt.var.monotone = NULL #restricts responses to individual predictors to monotone; rep(0, length(gbm.x))
brt.n.folds = {{ n_folds }} #number of folds
brt.prev.stratify = {{ prev_stratify }} #prevalence stratify the folds - only for presence/absence data
brt.family = "{{ family }}" #family - bernoulli (=binomial), poisson, laplace or gaussian
brt.n.trees = {{ n_trees }} #number of initial trees to fit
brt.step.size = brt.n.trees #numbers of trees to add at each cycle
brt.max.trees = {{ max_trees }} #max number of trees to fit before stopping
brt.tolerance.method = "{{ tolerance_method }}" #method to use in deciding to stop - "fixed" or "auto"
brt.tolerance = {{ tolerance_value }} #tolerance value to use - if method == fixed is absolute, if auto is multiplier * total mean deviance
brt.keep.data = FALSE #Logical. keep raw data in final model
brt.plot.main = FALSE #Logical. plot hold-out deviance curve
brt.plot.folds = FALSE #Logical. plot the individual folds as well
brt.verbose = FALSE #Logical. control amount of screen reporting
brt.silent = FALSE #Logical. to allow running with no output for simplifying model)
brt.keep.fold.models = FALSE #Logical. keep the fold models from cross valiation
brt.keep.fold.vector = FALSE #Logical. allows the vector defining fold membership to be kept
brt.keep.fold.fit = FALSE #Logical. allows the predicted values for observations from cross-validation to be kept

# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS")