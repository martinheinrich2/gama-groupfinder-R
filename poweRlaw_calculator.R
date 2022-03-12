# Power-law fit following the instructions in the poweRlaw package

# calculate x_min and alpha for all frequency files
# xmin is the lower cut-off x_min
# pars is a vector of parameter values (alpha)
# get ntail value and use it to calculate tail ratio of distributions

library (poweRlaw)
library (tools)
library (tictoc)
library (dplyr)
library (ptsuite)
library (qdapTools)
library (ggplot2)
tic()

# function to split filename into parts, remove "candidates"
SplitName <- function(name)  {
  name <- (file_path_sans_ext(basename(name)))
  name1 <- strsplit(name, "_")
  part1    <- unlist(name1)[2*(1:length(name))-1]
  part2    <- unlist(name1)[2*(1:length(name))  ]
  part3    <- unlist(name1)[2*(1:length(name))+1]
  SplitName <- paste0(part2, " km/s and ", part3, " Mpc")
  return(SplitName)
}

# function to load data
loadData <- function(filename) {
  filename <- paste(resultfolder, filename, sep="")
  richness <- read.csv(filename, header=TRUE, sep=",")
  return(richness)
}

# set path to data
resultfolder <- "./data/z_cmb/work_data/richness/"
outputfolder <- "./data/z_cmb/work_data/results/"

# create dataframes
pl_estimates <- data.frame()
tmp_est <- data.frame(los=numeric(), trans=numeric(), xmin=integer(), alpha=numeric(), maxrichness=numeric(), tfrac=numeric())

# get list of candidate files
temp = list.files(resultfolder, pattern = ".csv")

for (i in 1:length(temp)) {
  print(temp[i])
  filename <- temp[i]
  # split filename
  name <- (file_path_sans_ext(basename(filename)))
  name1 <- strsplit(filename, "_")
  param1    <- unlist(name1)[2*(1:length(filename))-1]
  param2    <- as.numeric(unlist(name1)[2*(1:length(filename))  ])
  param3    <- unlist(name1)[2*(1:length(filename))+1]
  param4    <- as.numeric(strsplit(param3, ".csv", fixed = TRUE))
  # extract parameters from filename
  param <- SplitName(filename)
  
  # load dataframe (i)
  myData <- loadData(filename)
  names(myData)[1] <- "richness"
  names(myData)[2] <- "count"
  newname <- paste0(param2, "_", sprintf("%.1f",param4))
  
  # calculate the 20/80 fraction of counts 
  # WARNING: might be better to calculate the real number of galaxies for each group (richness * counts)
  # scount is the sum of all groups, fcount is the sum of 80% of all groups
  scount <- sum(myData$count)
  eightycount <- round((4/5)*scount)
  twentycount <- round((1/5)*scount)
  fraccount <- (twentycount/eightycount)
  # calculate cumulative sums of counts
  myData$csum <- cumsum(myData$count)
  # get maximum richness in Data and use only Data with richness > 10
  maxrichness <- max(myData$richness)
  #x <- l
  if (maxrichness <= 11) {
    next()
  }
  
  # create subsample of counts 10 or less, use counts (N, number of members of each group)
  # set sample size to 10
  tailsampleN <- myData %>% filter(count <= 10)
  subsampleN <- myData %>% filter(count >= 10)
  sN <- sum(subsampleN$count)
  tN <- sum(tailsampleN$count)
  Ntail <- tN/sN
  iNtail <- sN/tN
  # use richness R for tail ratio
  lsample <- length(subsampleN$richness)
  tsample <- length(tailsampleN$richness)
  # number of R with X <= 10 divided by number of R with x > 10
  ltratio <- tsample/lsample
  
  # create subsample of richness 10 or less, use richness (N, number of members of each group)
  # set sample size to 10
  subsample <- myData %>% filter(richness <= 10)
  tailsample <- myData %>% filter(richness >= 11)
  subsample1 <- myData %>% filter(csum < eightycount)
  subsample1$csum1 <- cumsum(subsample1$count)
  tailsample1 <- myData %>% filter(csum > eightycount)
  tailsample1$csum1 <- cumsum(tailsample1$count)
  
  # generate a displ (discrete power law) object with empty model parameters "pars" and "no_pars"
  # from the vector "count" of the complete sample
  data_pl <- displ$new(myData$count)
  
  #estimate the x_min (lower cut-off) and exponent alpha of the power law and assign them to the power law object
  # x_min the lower theshold, for discrete distributions, x_min >0
  # estimate the lower bound, with respect to count
  # xmins and xmax are by default set to 1e+05 as the maximum value to explore
  # the default measure to calculate the distances is Kolomogorov Smirnoff statistic (KS)
  est <- estimate_xmin(data_pl)
  #est <- estimate_xmin(data_pl, pars = seq(1.5, 2.5, 0.01))
  # update the distribution object to the estimated lower bound Xmin
  data_pl$setXmin(est)
  # For a given value xmin, the scaling parameter is estimated by numerically optimising the log- likelihood.
  # The optimiser is initialised using the analytical MLE as stated in b_powerlaw_examples.pdf
  # plot the data (from x_min)
  #plot(data_pl)
  # Add the fitted distribution
  #lines(data_pl, col=2)
  
  # also possible to return the data invisibly for later use in ggplot2
  # dd = plot(data_pl)
  # head(dd, 3)
  
  # Clauset et al. reccomend a bootstrap procedure to get a handle on parameter uncertainty. 
  # check number of core with  parallel::detectCores()
  #bs = bootstrap(data_pl, no_of_sims=1000, threads=4)
  # The results of the bootstrap procedure can be investigated with histograms
  #hist(bs$bootstraps[,2], breaks="fd")
  #hist(bs$bootstraps[,3], breaks="fd")
  # and a bivariate scatter plot
  #plot(jitter(bs$bootstraps[,2], factor=1.2), bs$bootstraps[,3])
  
  # To determine whether the underlying distribution is a power-law we use the bootstrap p function
  ## This may take a while
  ## Use the mle to estimate the parameters
  #bs_p = bootstrap_p(data_pl, no_of_sims=1000, threads=4)
  #pvalue <- bs_p$p
  # A p-value - bs_p$p. For this example, p = 0.6738 which indicates that we can not rule out the power law model. 
  # See section 4.2 of the Clauset paper for further details.
  #hist(bs_p$bootstraps[,2], breaks="fd")
  #hist(bs_p$bootstraps[,3], breaks="fd")
  #plot(jitter(bs_p$bootstraps[,2], factor=1.2), bs_p$bootstraps[,3])
  
  #plot(myData)
  
  # calculate tail fraction
  ls <- length(subsample$count)
  lt <- length(tailsample$count)
  # tail/head
  tratio <- (lt)/(ls)
  
  # calculate ratio of tail counts / full sample counts
  tratio1 <- (lt)/(length(myData$count))
  
  # calculate tail ratio from estimated ntail, this is a variable tail-ratio
  lmyData <- length(myData$count)
  ntail <- get_ntail(data_pl)
  n <- get_n(data_pl)
  xtratio <- ntail/n
  #tfrac <- (sumtailsample)/(sumsubsample)
  #tfrac <- (l-x)/x
  #paretovalue1 <- pareto_test(myData$count)
  #paretovalue <- paretovalue1$`p-value`
  
  # the normalisation C=(alpha-1)*x_min^(alpha-1) from the continuous case!
  normalisation <- (est$pars-1)*est$xmin^(est$pars-1)
  # calculate sum of richness with less than 6 or 11 counts
  no5count <- sum(myData$count <= 5)
  no10count <- sum(myData$count <= 10)
  
  tmp_est <- cbind(param2, param4, est$xmin, est$pars, maxrichness, tratio, tratio1, xtratio, normalisation, no5count, no10count, Ntail, iNtail, ltratio)
  pl_estimates <- rbind(pl_estimates, tmp_est)
}

# get propper column names
names(pl_estimates)[1] <- "los"
names(pl_estimates)[2] <- "trans"
names(pl_estimates)[3] <- "R_min"
names(pl_estimates)[4] <- "alpha"
names(pl_estimates)[5] <- "max_R"
names(pl_estimates)[6] <- "t_ratio"
names(pl_estimates)[7] <- "t_ratio1"
names(pl_estimates)[8] <- "v_tratio"
names(pl_estimates)[9] <- "C"
names(pl_estimates)[10] <- "5_counts"
names(pl_estimates)[11] <- "10_counts"
names(pl_estimates)[12] <- "tail_ratio"
names(pl_estimates)[13] <- "tail_ratio1"
names(pl_estimates)[14] <- "tail_ratioR"

write.csv(pl_estimates, paste0(outputfolder, "poweRlaw_fit_new.csv"), row.names = F)

toc()

                                         