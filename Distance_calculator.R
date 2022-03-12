# Distancecalculator
# Calculate line-of-sight comoving and transverse comoving distance for each galaxy in data file and add columns with result
# Convert RA and DEC from degree to radians and add columns with result

# Load 'celestial' and 'tictoc' packages
library (celestial)
library (tictoc)

tic("processing data")
# set folder to data/
datafolder <- "./data/z_cmb/work_data/"

# Read csv into R as dataframe "MyData"
print("Reading data ...")
MyData <- read.csv(file=paste0(datafolder, "volume_limited_sample_z03.csv"), header=TRUE, sep=",")

#add new column "CoDist" for the line-of-sight (i.e. radial) comoving distance in Mpc
# use concordance '737 cosmology' with OmegaL = 0.7, OmegaM = 0.3, H0=70 
print("Calculating line-of-sight comoving distance ...")
MyData$CoDist <- cosdistCoDist(z=MyData$Z_CMB, H0=70, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1, wprime = 0, 737)

#add new column "CoDistTran" for the transverse comoving distance in Mpc
print("Calculating transverse comoving distance ...")
MyData$CoDistTran <- cosdistCoDistTran(z=MyData$Z_CMB, H0=70, OmegaM=0.3, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0 = -1, wprime = 0, 737)

# z_cmb is the redshift of the observed galaxies in the CMB frame
# define c as object in km/s and add new column for line of velocity v = c*ln(1+z_cmb)
# this is actually v = c * zeta
# for logarithmic shift zeta = ln(1+z) see also Baldry 2018a (https://arxiv.org/abs/1812.05135)
# log is natural logarithm (ln) in R
print("Calculating line-of-sight velocity ...")
c <- 299792.458
MyData$line_v <- (c*log(1+MyData$Z_CMB))

# caclulate zeta = ln(1+z_cmb)
print("Calculating zeta ...")
MyData$zeta <- (log(1+MyData$Z_CMB))

# convert RA from degree to radians and add new column
print("Converting RA from degree to radians ...")
MyData$RA_rad <- deg2rad(MyData$RA)

# convert DEC from degree to radians and add new column
print("Converting DEC from degree to radians ...")
MyData$DEC_rad <- deg2rad(MyData$DEC)

# copy selected columns to reduced sample
reducedData <- MyData[,c("CATAID","RA","DEC","Z_CMB","CoDist","CoDistTran","line_v","zeta","RA_rad","DEC_rad")]

# Write dataframe "MyData" to file as csv
print("Write data to disc ...")
write.csv(MyData, paste0(datafolder, "converted_sample_z03_zeta_new.csv"), row.names = F)
write.csv(reducedData, paste0(datafolder, "reduced_converted_sample_z03_zeta_new.csv"), row.names = F)
toc()