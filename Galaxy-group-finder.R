# find unique groups from all links files
# this version uses CATAID1 and CATAID2 of the first row to find groups

# read list of all files in 'links' folder and find groups for each file
library(tictoc)
library(dplyr)
library(tools)
library(compiler)

enableJIT(3)
rm()
tic("processing data")
# set folder to data/results
datafolder <- "./data/z_cmb/work_data/links/"
rfolder <- "./data/z_cmb/work_data/richness/"
ffolder <- "./data/z_cmb/work_data/frequency/"
gfolder <- "./data/z_cmb/work_data/groups/"

# initialise empty vector and dataframe
fgroups <- data.frame()
groupelement <- data.frame()

# function to split filename into parts, remove "candidates"
SplitName <- function(name)  {
  name1 <- strsplit(name, "_")
  part1    <- unlist(name1)[2*(1:length(name))-1]
  part2    <- unlist(name1)[2*(1:length(name))  ]
  part3    <- unlist(name1)[2*(1:length(name))+1]
  SplitName <- paste0(part2, "_", part3)
  return(SplitName)
}

# function to load data
loadData <- function(filename) {
  filename <- paste(datafolder, filename, sep="")
  links <- read.csv(filename, header=TRUE, sep=",")
  return(links)
}

# load dataframe
# get list of candidate files
temp = list.files(datafolder, pattern = ".csv")
# set groupnumber to 1
g <- 1
# loop over list of files to find groups
for (i in 1:length(temp)) {
    tic()
  #i <- 5
  g <- 1
  fgroups <- data.frame()
  #print(temp[i])
  filename <- temp[i]
  # split filename
  name <- (file_path_sans_ext(basename(filename)))
  name1 <- strsplit(filename, "_")
  param1    <- unlist(name1)[2*(1:length(filename))-1]
  param2    <- as.numeric(unlist(name1)[2*(1:length(filename))  ])
  param3    <- unlist(name1)[2*(1:length(filename))+1]
  param4    <- as.numeric(strsplit(param3, ".csv", fixed = TRUE))
  newname <- paste0(param2, "_", sprintf("%.1f",param4))
  rm(param1, param2, param3, param4)
  # load list of links in dataframe (i)
  links <- loadData(filename)
  
  # call function to extract parameters from filename
  #name1 <- SplitName (name)
  
  # get unique elements (remove double entries) and count number of entries
  uniquelinks <- dplyr::distinct(links, links$CATAID1, links$CATAID2)
  # set proper colum names
  names(uniquelinks)[1] <- "CATAID1"
  names(uniquelinks)[2] <- "CATAID2"
  unilinks<- arrange(uniquelinks, uniquelinks$CATAID1)
  rowcount <- nrow(unilinks)
  
  while (nrow(unilinks > 0)) {
  
    # test stuff to find groups
    # filter all identical rows with CATAID1 is equal CATAID1(row 1)
    uCATAID1 <- unilinks$CATAID1[1]
    uCATAID2 <- unilinks$CATAID2[1]
    neighbours1 <- unilinks %>% filter(CATAID1 == uCATAID1 | CATAID1 == uCATAID2)
    neighbours2 <- unilinks %>% filter(CATAID2 == uCATAID1 | CATAID2 == uCATAID2)
    # delete CATAID1 and CATAID2 that is in neighbours1 and neighbours2
    if (nrow(unilinks) >0) {
    unilinks <- unilinks[!(unilinks$CATAID1 %in% neighbours1$CATAID1),]
    }
    if (nrow(unilinks) >0) {
    unilinks <- unilinks[!(unilinks$CATAID2 %in% neighbours2$CATAID1),]
    }
    if (nrow(unilinks) ==0) {break()}
    # put all elements of dataframe into vector
    vneighbours1 <- unlist(neighbours1)
    vneighbours2 <- unlist(neighbours2)
    vneighbours1 <- unique(vneighbours1)
    vneighbours2 <- unique(vneighbours2)
    # combine all neighbours into one vector
    vallneighbours <- unique(c(vneighbours1, vneighbours2))
    # put vector into dataframe
    dallneighbours <- as.data.frame(vallneighbours)
    names(dallneighbours)[1] <- "CATAID"
    # delete neighbours from unilinks
    if (nrow(unilinks) > 0) {
    unilinks <- unilinks[!(unilinks$CATAID1 %in% dallneighbours$CATAID),]
    }
    if (nrow(unilinks) > 0) {
    unilinks <- unilinks[!(unilinks$CATAID2 %in% dallneighbours$CATAID),]
    }
      if (length(vallneighbours != 0)) {
        neighbours11 <- data.frame()
        neighbours12 <- data.frame()
        for (k in length(vallneighbours)) {
          neighbours11 <- c(neighbours11, (unilinks %>% filter(CATAID1 == vallneighbours[k])))
          neighbours12 <- c(neighbours12, (unilinks %>% filter(CATAID2 == vallneighbours[k])))
          vneighbours11 <- unlist(neighbours11)
          vneighbours12 <- unlist(neighbours12)
        }
      }
    
    vallneighbours <- unique(c(vallneighbours, vneighbours11, vneighbours12))
    groupelement <- as.data.frame(vallneighbours)
    if (nrow(unilinks) > 0) {
    unilinks <- unilinks[!(unilinks$CATAID1 %in% groupelement$CATAID),]
    }
    if (nrow(unilinks) > 0) {
    unilinks <- unilinks[!(unilinks$CATAID2 %in% groupelement$CATAID),]
    }
    names(groupelement)[1] <- "CATAID"
    groupelement$group <- g
    #if (nrow(unilinks) ==0) {break()}
    rm(neighbours1, neighbours2, neighbours11, neighbours12, vneighbours1, vneighbours11, vneighbours12, vneighbours2, vallneighbours)
    # increment groupcounter
    #if (nrow(unilinks) == 0) {break()}
    fgroups <- rbind(fgroups, groupelement)
    #fgroupsDT <- rbindlist(list(fgroupsDT, groupelement))
    
    g <- g+1
  }
  
  # calculate number of links per group and frequency
  groupcount <- fgroups %>% count(group, sort = FALSE)
  names(groupcount)[1] <- "group"
  names(groupcount)[2] <- "members"
  
  freqcount <- groupcount %>% count(members, sort = FALSE)
  names(freqcount)[1] <- "richness"
  names(freqcount)[2] <- "count"
  
  print(filename)
  toc()
  
  write.csv(groupcount, paste0(ffolder, "frequency_", newname, ".csv"), row.names = F)
  write.csv(freqcount, paste0(rfolder, "richness_", newname, ".csv"), row.names = F)
  write.csv(fgroups, paste0(gfolder, "groups_", newname, ".csv"), row.names = F)
  #break()
  rm(links, unilinks, uniquelinks, groupcount, fgroups, uCATAID1)
}

toc()