# this script recreates Simpkins et al 2013 analysis for tadpoles
# Simpkins C, Shuker J.D., Lollback G.W., Castley J.G., Hero J. (2013).
#	 Environmental variables associated with the distribution and occupancy 
#	of habitat specialist tadpoles in naturally acidic, oligotrophic 
#	waterbodies. Austral Ecology (in press).

# set the working directory
wd = "c:/userdata/FRM/shuker/"
setwd(wd)

species = c("Litoria_olongburensis", "Crinia_tinnula")

# read in the species and env data
all.data = read.csv(paste(wd, "shuker.csv", sep=""), header=TRUE)

# get the names of the columns to use as variable names
predictors = colnames(all.data)[7:11]

# create a list of formulas 
glm.explanatories = 
	c("salinity", 
	"salinity + turbidity", 
	"salinity + turbidity + depth",
	"salinity + turbidity + depth + percent_cover",
	"salinity + turbidity + depth + percent_cover + predatory_fish",
	"turbidity", 
	"turbidity + depth",
	"turbidity + depth + percent_cover",
	"turbidity + depth + percent_cover + predatory_fish",
	"depth",
	"depth + percent_cover",
	"depth + percent_cover + predatory_fish",
	"percent_cover",
	"percent_cover + predatory_fish",
	"predatory_fish")


######################################################
##
##	1) Assess the importance of environmental variables 
##		on the relative abundance of tadpoles
##
######################################################

# create a list to hold the model output
sp.abund.models = list()

# fit regression to abundance data for each species and model formula
for (sp in species) {

	# get the index of the species
	i = which(species == sp)

	# create a list for outputs of each glm
	sp.abund.models[[i]] = list()

	# for each set of explanatories
	for (j in 1:length(glm.explanatories)) {

		# generate the model formula
		my.resp = paste(sp, "_A", sep="")
		my.formula = paste(my.resp, "~", glm.explanatories[j])

		# fit the model
		sp.abund.models[[i]][[j]] = glm(formula=my.formula, data=all.data,
			family=poisson)
	}

}


######################################################
##
##	2) Assess the importance of environmental variables 
#		on the occupancy of tadpoles
##
######################################################

# create a list to hold the model output
sp.occur.models = list()

# fit regression to occurrence data for each species
for (sp in species) {

	# get the index of the species
	k = which(species == sp)

	# create a list for outputs of each glm
	sp.occur.models[[k]] = list()

	# for each set of explanatories
	for (l in 1:length(glm.explanatories)) {

		# generate the model formula
		my.occur.resp = paste(sp, "_O", sep="")
		my.occur.formula = paste(my.occur.resp, "~", glm.explanatories[l])

		# fit the species model
		sp.occur.models[[k]][[l]] = glm(formula=my.occur.formula, data=all.data,
			family=binomial)
	}
}


######################################################
##
##	Model selection
##
######################################################


# rank the models using second order AICc
library(AICcmodavg)

# create model selection table for each species
for (s in 1:length(species)) {

	A_comparison = aictab(sp.abund.models[[s]], modnames=glm.explanatories)
	write.csv(A_comparison, file=paste(species[s], "_abundance_comparison.csv", sep=""))

	O_comparison = aictab(sp.occur.models[[s]], modnames=glm.explanatories)
	write.csv(O_comparison, file=paste(species[s], "_occurrence_comparison.csv", sep=""))
}


######################################################
##
##	Parameter estimates
##
######################################################

best.abund.model = which(glm.explanatories==A_comparison[[1]][1])
summary(sp.abund.models[[1]][best.abund.model][[1]])

best.occur.model = which(glm.explanatories==O_comparison[[1]][1])
summary(sp.occur.models[[1]][best.occur.model][[1]])