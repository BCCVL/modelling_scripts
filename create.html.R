# this script will combine's R's output files into a single html doc

install.packages(c("R2HTML", "png"))
library(R2HTML) 
library(png)

function generateHTML(sp, outputdir) {

# define the species
#species = "ABT"

# define where the model outputs are stored
#wd = "/home/jc140298/sprint2/ABT/output_bioclim/"

setwd(outputdir)

# read in model outputs
auccurve = readPNG(paste(outputdir, "AUC.png", sep=""))
accuracystats = read.csv(paste(outputdir, "combined.modelEvaluation.csv", sep=""),	row.names=c(1))

# create the output file 
target <- HTMLInitFile(outdir=outputdir, filename=paste(sp,"_output", sep=""), BackGroundColor="#CCCCCC")

# add content
HTML(paste("<center><br><H1>Model Output for ", sp, sep=""), file=target)

HTML("<br><H2>AUC:ROC curve", file=target)
HTMLInsertGraph("AUC.png", file=target)

HTML("<br><H2>Accuracy measures",file=target)
HTML(accuracystats, file=target)

# close the file
HTMLEndFile()

}