# this script recreates Keatley and Hudson analysis for flowing plants
# Keatley, M.R., and Hudson, I.L. (2012). Detecting change in an Australian 
#	flowering record: Comparisons of linear regression and cumulative sum 
	analysis change point analysis. Austral Ecology 37, 825-835.
# Keatley, M.R., and Hudson, I.L. (2007). Shift In Flowering Dates 
#	Of Australian Plants Related To Climate: 1983-2006.

# set the working directory
wd = "c:/userdata/FRM/lynda/"
setwd(wd)

# read in the species occurrence and first flowering day data
species.info = read.csv(paste(wd, "flowers_FFD.csv", sep=""))

#                    species     lon    lat yearday      FFD
#1      Comesperma volubile 145.415 37.996    1983       NA
#2      Comesperma volubile 145.415 37.996    1984 279.2265
#3      Comesperma volubile 145.415 37.996    1985 283.0939
#4      Comesperma volubile 145.415 37.996    1986 262.9282
#5      Comesperma volubile 145.415 37.996    1987 283.0939

# get a list of flower names
species = unique(species.info$species)

# load the libraries required
library(MASS)	# Stepwise Regression

# read in the env data created by lynda_data.R (downloaded from BOM and
#	manipulated accordingly)
env.data = read.csv(paste(wd, "env_vars.csv", sep=""), header=TRUE, 
	row.names=1)

#year annual.rain Jan_rain Feb_rain March_rain April_rain May_rain June_rain
#1983       830.2     54.2     13.8       55.8       49.2     87.2      69.0
#1984       898.0     67.2     44.2       74.8       72.2     30.8      53.0
#1985       859.2     20.8      7.6       37.0       93.6     83.1      76.2
#1986       756.4     36.8     25.0       11.0       86.8    109.2      52.8
#1987       820.3     60.4     55.0       76.1       37.6    110.9      80.0
#1988       822.3     60.6     22.4       35.4       28.4     85.7      94.8


# get the names of the columns to use as variable names
predictors = colnames(env.data)


######################################################
##
##	Is there a change in FFD through time?
## 	Linear Regression
##
######################################################

# for each species, plot the flowering through time regressions

# create a list to hold the model output
sp.lm.models = list()

# create a matrix to hold the model parameter estimates and impt info
# see Table 1 (Keatley and Hudson 2012)
lm.output = matrix(NA, nrow=length(species), ncol=7)
colnames(lm.output) = c("Mean date", "SD", "n", "R2", "P-value", "Shift_tot", 	
	"Shift_dpy")

# fit regression for each species
for (sp in species) {

	# get the index of the species
	i = species[which(species == sp)]

	# get that species' data
	sp.data = species.info[species.info$species == sp,]

	# fit the species model
	sp.lm.models[[i]] = lm(sp.data$FFD~sp.data$yearday)

	# plot the relationship
	# need to remove the NA's
	notNA <- !is.na(sp.data$FFD) 
	windows()
	plot(sp.data$yearday[notNA], sp.data$FFD[notNA], main=sp, xlab="", 
		ylab="", pch=19, xlim=c(1980, 2010), 
		ylim=c(min(sp.data$FFD[notNA])-10, max(sp.data$FFD[notNA])+10))
	# add the best fit line to the plot
	# EMG might user want to see these figures?
	abline(sp.lm.models[[i]])

	# calculate the statistics in Table 1 (Keatley and Hudson 2012)
	meandate = mean(sp.data$FFD, na.rm=TRUE)
	sddate = sd(sp.data$FFD, na.rm=TRUE)
	n = sum(!is.na(sp.data$FFD))
	r2 = summary(sp.lm.models[[i]])$r.squared
	pvalue = pf(summary(sp.lm.models[[i]])$fstatistic[1], 
		summary(sp.lm.models[[i]])$fstatistic[2], 
		summary(sp.lm.models[[i]])$fstatistic[3],
	      lower.tail = FALSE)

	# use best fit line to estimate shift in flowering
	modelpredictions = predict(sp.lm.models[[i]])
	shift = modelpredictions[1] - modelpredictions[length(modelpredictions)]
	peryear = shift/abs(min(sp.data$yearday)-max(sp.data$yearday))
	
	# save the species output
	lm.output[i,] = c(meandate, sddate, n, r2, pvalue, shift, peryear)
}

# combine species output with species names
df.lm.output = data.frame(species, lm.output)

# save 
write.csv(df.lm.output, file=paste(wd, "Response_output.csv", sep=""), 
	row.names=FALSE)


######################################################
##
##	Is there a change in env variables through time?
##	Linear Regression
##
######################################################

# for each variable, plot the flowering through time regressions

# create a list to hold the model output
env.lm.models = list()

# create a matrix to hold the model parameter estimates and impt info
env.lm.output = matrix(NA, nrow=length(predictors), ncol=4)
colnames(env.lm.output) = c("n", "P-value", "Shift", "R2") 	

years = as.numeric(rownames(env.data))

# fit regression for each predictor (enviromental) variable
for (env in predictors) {

	# get the index of the variable
	i = which(predictors == env)

	# fit the model
	env.lm.models[[i]] = lm(env.data[,i] ~ years)

	# plot the model
	windows()
	notNA <- !is.na(env.data[,i]) 
	plot(years[notNA], (env.data[,i])[notNA], main=env, xlab="", ylab="", 
		pch=19, xlim=c(1980, 2010), ylim=c(min((env.data[,i][notNA]))-10, 
		max((env.data[,i][notNA]))+10))
	abline(env.lm.models[[i]])
	
	n = sum(!is.na(env.data[,i]))
	pvalue = pf(summary(env.lm.models[[i]])$fstatistic[1], 
		summary(env.lm.models[[i]])$fstatistic[2], 
		summary(env.lm.models[[i]])$fstatistic[3],
	      lower.tail = FALSE)
	shift = as.numeric(coef(env.lm.models[[i]])[2])
	r2 = summary(env.lm.models[[i]])$r.squared

	
	env.lm.output[i,] = c(n, pvalue, shift, r2)
}

df.env.lm.output = data.frame(predictors, env.lm.output)
write.csv(df.env.lm.output, file=paste(wd, "Environment_output.csv", sep=""),
	row.names=FALSE)

######################################################
##
##	Which environmental variables are most important?  
##	Multiple linear regression
##
######################################################

# first make a list of important variables from Envrionment_output.csv
# get the output
env.out = read.csv(paste(wd, "Environment_output.csv", sep=""), row.names=1)

# search for pvalues less than 0.5
imptvars = rownames(env.out[env.out$P.value < 0.05,])

# create a list to hold output of multiple regression model
mlm.output.fit = list()
mlm.output.step = list()

# for each species, determine which variables are most important
for (sp in species) {

	i = species[which(species == sp)]

	# get single species data
	sp.data = species.info[species.info$species == sp,]

	# need a single data frame to fit the model
	mlm.data = cbind(sp.data$FFD, env.data[,imptvars])

	# create module formula with impt variables only
	mlm.formula <- as.formula(paste("sp.data$FFD ~", paste(colnames(mlm.data[-1]), 
		collapse="+")))

	# fit multiple regression model
	mlm.output.fit[[i]] <- lm(mlm.formula, data=mlm.data)
summary(mlm.output.fit)
	mlm.output.step[[i]] <- stepAIC(mlm.output.fit[[i]], direction="forward")
summary(mlm.output.step)

#	step$anova # display results
}
# EMG only one significant variable for one species
