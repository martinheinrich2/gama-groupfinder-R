# Galaxy-pair-finder
# find galaxy pairs within given line-of-sight and transverse linking lengths
# and write link-number, CATAID1, CATAID2 in file
# Load 'celestial' and 'tictoc' packages
library (compiler)
enableJIT(3)

library (celestial)
library (tictoc)
library (plyr)
library (data.table)

rm()
# set linking parameters line-of-sight in km/s and transverse in Mpc
ll_los <- 500
# loop over sequence of transverse separations in Mpc
for (ll_trans in seq(0.1,2,0.1)){    
    # start clock
    tic("processing data")
    # set folder to data/retuls
    datafolder <- "./data/z_cmb/work_data/"
    resultfolder <- "./data/z_cmb/work_data/galaxy-links/"
    
    # Read CSV into R as dataframe and make temporary copy
    print(paste("Reading reduced converted sample ... processing data:"))
    MyData <- read.csv(file=paste0(datafolder, "reduced_converted_sample_z03_zeta.csv"), header=TRUE, sep=",")
    workData <- MyData
    tmpData <- workData
    
    # create dataframe for links
    links <- data.frame()
    linksDT <- data.table()
    linklist <- list()
    
    # assign values to i and link
    i <- 1
    link <- 0
    
    # get number of rows in dataframe
    rowcount <- nrow(MyData)
    
    # Function to calculate line-of-sight separation
    CalcSepLos <- function(line_v1, line_v2){
      sep_los <- abs(line_v1-line_v2)
      return(sep_los)
    }
    
    # Function to calculate sky separation (angular distance) of two objects as angle theta in degrees
    # DEC and RA must be in radians
    CalcSkySep <- function(DEC1,DEC2,RA1, RA2){
      # use pmin() to ensure that the precision loss leads to a value no higher than exactly 1, see:
      # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
      sky_sep <- acos(pmin((sin(DEC1)*sin(DEC2))+(cos(DEC1)*cos(DEC2))*cos(RA1-RA2),1))
      return(sky_sep)
    }
    
    # Function to add separations and test linking conditions, returns dataframe 'tmpData'
    ConCheck <- function(x){
      # Add column of line-in-sight separation
      tmpData$seplos <- (CalcSepLos(tmpData$line_v[x],tmpData$line_v))
      # separation = theta * (CoDistTran[x] + CoDistTran)/2
      # theta <- acos(sin(DEC1)sin(DEC2)+cos(DEC1)cos(DEC2)cos(RA1-RA2))
      # theta <- CalcSkySep(tmpData$DEC_rad[x],tmpData$DEC_rad,tmpData$RA_rad[x],tmpData$RA_rad)
      tmpData$septrans <- (CalcSkySep(tmpData$DEC_rad[x],tmpData$DEC_rad,tmpData$RA_rad[x],tmpData$RA_rad))*(tmpData$CoDistTran[x]+tmpData$CoDistTran)/2
      # Add column to check if both conditions are satisfied (if only one true result, then there are no neighbours)
      tmpData$truefalse <- (tmpData$seplos <= ll_los & tmpData$septrans <= ll_trans)
      #SumTrueFalse <- sum(tmpData$truefalse)
      return(tmpData)
    }
    
    # loop over all rows in the dataframe and find galaxies within the linking lenght separations
    while (i <= rowcount) {
      candidateID <- tmpData$CATAID[i]
      # Calculate separations and test linking condition
      tmpData <- (ConCheck(i))
      trueID <- which(tmpData$truefalse == 1)
      CATAID <- tmpData$CATAID[trueID]
      IDlength <- length(CATAID)
      # ignore single galaxies
      if (IDlength >= 2) {
        j <- 2
        while (j <= IDlength) {
          link <- link+1
          list_element <- as.list(append(link, (CATAID[c(1,j)])))
          # use rbindlist from data.table package instead of rbind for faster code
          linksDT <- rbindlist(list(linksDT, list_element))
          j <- j+1
        } 
      }
      cat("\r", "links found:", link, sep="\t" ,"galaxies checked:", i, "of ", rowcount, "  ")
      i <- i+1
    }
    links <- as.data.frame(linksDT)
    #links$truefalse <- NULL
    # get propper column names
    names(links)[1] <- "link"
    names(links)[2] <- "CATAID1"
    names(links)[3] <- "CATAID2"
    resultnameL <- paste0(resultfolder, "links_", ll_los, "_", sprintf("%.1f", ll_trans), ".csv")
    write.csv(links, paste0(resultnameL), row.names = F)
    print("")
    # display elapsed time
    toc()  
}
