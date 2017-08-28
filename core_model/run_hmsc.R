# Set WD
setwd("~/Bureau/Documents/Workshops/metaco2.0/metaco_toy_model")

# Run simulations
source("run_scenarios.R")

# Run the analysis
XY <- read.table("XY.txt")
E <- read.table("E.txt")

# Wrapper to run the analysis
run_hmsc <- function(Y, X, XY, D, title) {

	library(HMSC)
	library(scatterplot3d)

	# Prepare the data
	siteID = as.factor(c(1:nrow(X)))
	coord <- data.frame(siteID = siteID, coordX = XY[,1], coorY  = XY[,2])

	# Fit the model with space as an autocorrelation and site as a random effect

	# With environment
	data_full <- as.HMSCdata(Y = Y, X = X, Random = siteID, Auto = coord)
	model_full <- hmsc(data_full, family = "probit", niter = 10000, nburn = 1000, thin = 10)

	# Without environment
	data_E0 <- as.HMSCdata(Y = Y, Random = siteID, Auto = coord)
	model_E0 <- hmsc(data_E0, family = "probit", niter = 10000, nburn = 1000, thin = 10)

	# Test for convergence
	mixingParamX <- as.mcmc(model_full, parameters = "paramX")
	tmp<-geweke.diag(mixingParamX)
	pvalue2sided<-2*pnorm(-abs(tmp$z))

	# Variance partitioning
	VP <- variPart(model_full, groupX = c("Intercept",rep("E", D)))

	# Plot the results of variation partitioning
	dev.new(width = 10, height = 10)
	library(vcd)
	ternaryplot(
		VP[,c(2,4,3)], 
		pch = 19, 
		bg = "lightgray",
		grid_color = "white",
		col = "darkblue",
		main = title,
		dimnames = c("Environment", "Spatial autocorrelation", "Co-distribution")
	)

	# Plot the heat maps
	library(corrplot)
	library(viridis)
	pal <- viridis(200)

	dev.new(width = 10, height = 10)
	codistr_full <- apply(corRandomEff(model_full)$Random, 1:2, mean)
	corrplot(codistr_full, method = "color", col = pal, type = "lower",
	diag = FALSE, tl.srt = 45)
	title(paste(title, "Controlling for E"), cex = 2)

	dev.new(width = 10, height = 10)
	codistr_E0 <- apply(corRandomEff(model_E0)$Random, 1:2, mean)
	corrplot(codistr_E0, method = "color", col = pal, type = "lower",
	diag = FALSE, tl.srt = 45)
	title(paste(title, "Not controlling for E"), cex = 2)

	# Return the resultsw
	return(list(VP, pvalue2sided, model_full))
}


res1 <- run_hmsc(Y = Y1, X = E, XY = XY, D = 1, title = "Metapopulation")

res2 <- run_hmsc(Y = Y2, X = E, XY = XY, D = 1, title = "Sorting, no competition")

Y3sum <- apply(Y3,2,sum)
Y3 <- Y3[,Y3sum!=0]
res3 <- run_hmsc(Y = Y3, X = E, XY = XY, D = 1, title = "Sorting with competition")

Y4sum <- apply(Y4,2,sum)
Y4 <- Y4[,Y4sum!=0]
res4 <- run_hmsc(Y = Y4, X = E, XY = XY, D = 1, title = "Neutral")
