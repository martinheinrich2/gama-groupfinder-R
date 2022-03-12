# Power-law fit following the instructions in the poweRlaw package

# calculate x_min and alpha for all frequency files
# xmin is the lower cut-off x_min
# pars is a vector of parameter values (alpha)

library (tools)
library (tictoc)
library (dplyr)
library (ptsuite)
library (qdapTools)
library (ggplot2)
library (ggpubr)
tic()

# set path to data
resultfolder <- "./data/z_cmb/work_data/results/"
outputfolder <- "./data/z_cmb/work_data/results/"

shortname <- "poweRlaw_fit.csv"
filename <- paste(resultfolder, shortname, sep="")
MyData <- read.csv(filename, header=TRUE, sep=",")


names(MyData)[3] <- "R"
# create data blocks smaller 600
#df600 <- MyData %>% filter(los <= 600)


#dfx5 <- df750 %>% filter(X5_counts < 10)
#dfx10 <- df750 %>% filter(X10_counts < 10)

# plot los, trans and alpha as heatmap
ggplot(data = MyData, mapping = aes(x = los, y = trans, fill = alpha)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and scaling parameter alpha")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and normalisation constant as heatmap
ggplot(data = MyData, mapping = aes(x = los, y = trans, fill = C)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths in km/s and Mpc and normalisation constant C")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and R_min as richness
ggplot(data = MyData, mapping = aes(x = los, y = trans, fill = R)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and R_min as richness")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and tratio as heatmap
# number of R with X <= 10 divided by number of R with x > 10
ggplot(data = MyData, mapping = aes(x = los, y = trans, fill = t_ratio)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and tail ratio, number of R with X<=10/number of R with X>10")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and tratio as heatmap
# Ntratio = number of counts X where R has X<=10 divided by number of X where R has X>10 
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = Ntratio)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and tail ratio, number of X where R has X<=10/number of X where R has X>10")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and tratio as heatmap
# number of groups with R>10/number of groups with R <=10
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = tratio)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and number of groups with R > 10 divided by number of groups with R <=10")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and tratio as heatmap
# number of groups with R>10/number of all groups
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = tratio1)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and number of groups with R > 10 divided by number of all groups")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and tratio as heatmap
# number of R>10/number of groups within range 5--10
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = tratio2)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and number of groups with R > 10 divided by number of groups within range 5--10")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and normalisation as heatmap
ggplot(data = df900, mapping = aes(x = los, y = trans, fill = normalisation)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths in km/s and Mpc and normalisation constant C")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and Rmin as heatmap
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = R)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and richness R for which the number of groups is R_min")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and maxrichness as heatmap
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = (max_R-R))) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and maximum richness - Rmin")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and maxrichness as heatmap
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = max_R)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and maximum richness")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and sum of richness where counts are less than 6 as heatmap
ggplot(data = df800, mapping = aes(x = los, y = trans, fill = X5_counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and # of counts >= 5")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and sum of richness where counts are less than 6 as heatmap
ggplot(data = MyData, mapping = aes(x = los, y = trans, fill = X5_counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and maximum richness")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and sum of richness where counts are less than 6 as heatmap
ggplot(data = dfx5, mapping = aes(x = los, y = trans, fill = X5_counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and # of counts >= 5")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and sum of richness where counts are less than 6 as heatmap
ggplot(data = df900, mapping = aes(x = los, y = trans, fill = X5_counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and # of counts <= 6")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and sum of richness where counts are less than 11 as heatmap
ggplot(data = dfx10, mapping = aes(x = los, y = trans, fill = X10_counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and # of counts <= 10")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")

# plot los, trans and sum of richness where counts are less than 11 as heatmap
ggplot(data = df900, mapping = aes(x = los, y = trans, fill = X10_counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and # of counts <= 10")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")


# plot los, trans and xmin as maxrichness
ggplot(data = MyData, mapping = aes(x = los, y = trans, fill = max_R)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_tile() +
  ggtitle(paste("Linking lengths and maxrichness")) +
  labs(x="line-of-sight in km/s", y="transverse in Mpc")
