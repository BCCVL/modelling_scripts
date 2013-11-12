# This is a Python string template file and is currently tightly tied
# to ../maxent.py.
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

model.maxent = TRUE #boolean to run MAXENT algorithm
project.maxent = TRUE #boolean to project MAXENT algorithm
evaluate.maxent = TRUE #boolean to evaluate MAXENT algorithm

#additional parameters to set for modelling maxent
maxent.jar = "/home/jc165798/working/BCCVL/maxent.jar" #define location of maxent.jar file
maxent.outputformat = "logistic" #options include logistic, cumulative, raw
maxent.randomseed =FALSE
maxent.logscale =TRUE
maxent.removeduplicates =TRUE
maxent.randomtestpoints = {{ random_test_points }} #default 0
maxent.betamultiplier = {{ beta_multiplier }} #default 1
maxent.maximumbackground = {{ maximum_background }} #default 10000
maxent.biasfile = NULL
maxent.testsamplesfile = NULL
maxent.replicates = {{ replicates }} #default 10
maxent.replicatetype = "{{ replicate_type }}" #default "crossvalidate" #options include crossvalidate, bootstrap, subsample
maxent.linear = TRUE
maxent.quadratic = TRUE
maxent.product = TRUE
maxent.threshold = TRUE
maxent.hinge = TRUE
maxent.addsamplestobackground = TRUE
maxent.addallsamplestobackground = FALSE
maxent.fadebyclamping = FALSE
maxent.extrapolate = TRUE
maxent.autofeature = TRUE
maxent.doclamp = TRUE
maxent.maximumiterations = {{ maximum_iterations }} #default 500
maxent.convergencethreshold = {{ convergence_threshold }} #default 1.00E-05
maxent.lq2lqptthreshold = {{ lq2lqpt_threshold }} #default 80
maxent.l2lqthreshold = {{ l2lq_threshold }} #default 10
maxent.hingethreshold = {{ hinge_threshold }} #default 15
maxent.beta_threshold = {{ beta_threshold }} #default -1
maxent.beta_categorical = {{ beta_categorical }} #default -1
maxent.beta_lqp = {{ beta_lqp }} #default -1
maxent.beta_hinge = {{ beta_hinge }} #default -1
maxent.defaultprevalence = {{ default_prevalence }} #default 0.5
maxent.nodata = {{ no_data }} #default -9999

# EMG there is also an additional argument 'args' used to pass arguments (options) to the maxent software for projection

# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth = c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS")