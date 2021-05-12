###  Prepare empirical data on migratory connectivity to validate ORSIM   ###

library(dggridR)
library(rworldmap)
library(rgdal)
library(rgeos)


##  Construct an hexagon grid covering the Americas  ##
hexgrid <- dgconstruct(projection="ISEA", topology="HEXAGON", res=7, metric=T)
hexgrid_center <- dgSEQNUM_to_GEO(hexgrid, 1:21872) # 65612 / 590492
hexgrid_centroids <- cbind(hexgrid_center$lon_deg, hexgrid_center$lat_deg)
hex_sel <- which(hexgrid_centroids[,1] < -30 & hexgrid_centroids[,1] > -170 & hexgrid_centroids[,2] > -60 & hexgrid_centroids[,2] < 80)
hexgrid2 <- dgcellstogrid(hexgrid, hex_sel, frame=F,wrapcells=TRUE)
hexgrid2_centroids <- matrix(unlist(lapply(hexgrid2@polygons, function(x) x@labpt)), byrow=T, ncol=2)
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, proj4string(hexgrid2))
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- gBuffer(newmap, byid=TRUE, width=0)
newmap <- gUnaryUnion(newmap, id=newmap@data$world)
hexgrid3 <- gIntersection(hexgrid2, newmap, byid=T)
hexgrid3_centroids <- matrix(unlist(lapply(hexgrid3@polygons, function(x) x@labpt)), byrow=T, ncol=2)
sr <- proj4string(hexgrid3)

##  Lists of species names for the various datasets  ##
species.names <- read.csv("resources/species_names.csv")

##  Seasonal abundance estimates  ##
summer.abundance <- read.csv("results/output/seasonalAbundance_BR.csv", header=T)
winter.abundance <- read.csv("results/output/seasonalAbundance_NB.csv", header=T)

##  Function to get banding data, which have been downloaded and saved on my local disk  ##
banding_data_fct <- function(x){
  ptsBR <- vector()
  ptsNB <- vector()
  for(i in 1:length(x$GISBLONG)){
    if(x$B_SEASON[i] == "B"){
      ptsBR = rbind(ptsBR, c(x$GISBLONG[i], x$GISBLAT[i]))
    }else{
      ptsBR = rbind(ptsBR, c(x$GISRLONG[i], x$GISRLAT[i]))
    }	
  }
  for(i in 1:length(x$GISBLONG)){
    if(x$R_SEASON[i] == "B"){
      ptsNB = rbind(ptsNB, c(x$GISBLONG[i], x$GISBLAT[i]))
    }else{
      ptsNB = rbind(ptsNB, c(x$GISRLONG[i], x$GISRLAT[i]))
    }	
  }
  pts <- list(ptsBR, ptsNB)
  names(pts) <- c("ptsBR", "ptsNB")
  return(pts)
}

##  Function to get hexagons containing the empirical breeding and wintering locations of individuals  ##
empirical_hexagons_fct <- function(x, spp_num){
  dfBR = data.frame(a = 1:length(x$ptsBR[,1]))
  sp.ptsBR <- SpatialPointsDataFrame(x$ptsBR, dfBR)
  proj4string(sp.ptsBR) <- sr
  empiricalData_breedingHexagons <- over(sp.ptsBR, hexgrid3[which(summer.abundance[,spp_num]>0),])
  empiricalData_breedingHexagons2 <- unique(empiricalData_breedingHexagons[which(is.na(empiricalData_breedingHexagons)==FALSE)])
  dfNB = data.frame(a = 1:length(x$ptsNB[,1]))
  sp.ptsNB <- SpatialPointsDataFrame(x$ptsNB, dfNB)
  proj4string(sp.ptsNB) <- sr
  empiricalData_nonbreedingHexagons <- over(sp.ptsNB, hexgrid3[which(winter.abundance[,spp_num]>0),])
  empiricalData_nonbreedingHexagons2 <- unique(empiricalData_nonbreedingHexagons[which(is.na(empiricalData_nonbreedingHexagons)==FALSE)])
  
  ##  Only keep individuals for which both breeding and wintering locations fall within an hexagon occupied by the species (i.e. with a seasonal relative abundance > 0) ##
  toKeep <- which(is.na(empiricalData_breedingHexagons) == F & is.na(empiricalData_nonbreedingHexagons) == F)
  empiricalData_breedingHexagons3 <- empiricalData_breedingHexagons[toKeep]
  empiricalData_nonbreedingHexagons3 <- empiricalData_nonbreedingHexagons[toKeep]
  x$ptsBR <- x$ptsBR[toKeep,]
  x$ptsNB <- x$ptsNB[toKeep,]
  x$type <- x$type[toKeep]
  x$ptsBR_hex <- empiricalData_breedingHexagons3
  x$ptsNB_hex <- empiricalData_nonbreedingHexagons3
  return(x)
}


###  Process data for all species investigated  ###

##  1 - Wood Thrush  ##
spp_number <- 1
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/WOTH_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "wood thrush"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
#pts$ptsBR <- apply(pts$ptsBR, 2, function(x) as.numeric(as.character(x)))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
#pts$ptsNB <- apply(pts$ptsNB, 2, function(x) as.numeric(as.character(x)))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type, pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex)

##  2 - Swainson's Thrush  ##
spp_number <- 2
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/SWTH_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "swainson's thrush"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  3 - Hermit Thrush  ##
spp_number <- 3
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/HETH_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "hermit thrush"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  4 - American Robin  ##
spp_number <- 4
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/AMRO_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- rep("banding", nrow(bandingData))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  5 - Red winged Blackbird  ##
spp_number <- 5
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/RWBL_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- rep("banding", nrow(bandingData))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  6 - Yellow Warbler  ##
spp_number <- 6
spp <- as.character(species.names$ebird_code[spp_number])
# Genetic data
wintering.assignments <- readRDS("resources/Genoscapes/YellowWarbler/YWAR_cleaned.rds")
wintering.assignments <- wintering.assignments[which(wintering.assignments$Stage=="Wintering"),]
breeding.surfaces.origen <- readRDS("resources/Genoscapes/YellowWarbler/BreedingSurfaces.rds")
wintering.surfaces.origen <- readRDS("resources/Genoscapes/YellowWarbler/WinterProbSurfaces.rds") 
breeding_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
}
breeding_origen <- breeding_origen[-c(1,2)]
wintering_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  if(length(which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])) > 0){
    wintering_origen[k] <- which(wintering.assignments$Sample == colnames(wintering.surfaces.origen)[k])
  }else{
    wintering_origen[k] <- NA
  }
}
wintering_origen <- wintering_origen[-c(1,2)]
toRemove <- which(is.na(wintering_origen) == T)
wintering_origen <- wintering_origen[-toRemove]
breeding_origen <- breeding_origen[-toRemove]
pts <- list()
pts$ptsBR = apply(wintering.surfaces.origen[,1:2][breeding_origen,], 2, function(x) as.numeric(as.character(x)))
pts$ptsNB <- as.matrix(wintering.assignments[,3:4][wintering_origen,])
pts$type <- rep("genetics", nrow(pts$ptsNB))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  7 - Common Yellowthroat  ##
spp_number <- 7
spp <- as.character(species.names$ebird_code[spp_number])
# Genetic data
wintering.surfaces.origen <- readRDS("resources/Genoscapes/CommonYellowthroat/COYE_WinterProbSurfaces.rds")
winter.data <- read.csv("resources/Genoscapes/CommonYellowthroat/WinterData.csv")
breeding_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
}
breeding_origen <- breeding_origen[-c(1,2)]
wintering_origen <- winter.data[c(5,6)]
pts <- list()
pts$ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
pts$ptsNB <- winter.data[c(5,6)]
pts$type <- rep("genetics", nrow(pts$ptsNB))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                as.matrix(cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex)))


##  8 - Wilson's Warbler  ##
spp_number <- 8
spp <- as.character(species.names$ebird_code[spp_number])
# Genetic data
wintering.assignments <- readRDS("resources/Genoscapes/WilsonsWarbler/WIWA.WinterAssignment.rds")
wintering.surfaces.origen <- readRDS("resources/Genoscapes/WilsonsWarbler/WIWA_WinterProbSurfaces.rds")
breeding_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
}
breeding_origen <- breeding_origen[-c(1,2)]
wintering_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
    wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
  }else{
    wintering_origen[k] <- NA
  }
}	
wintering_origen <- wintering_origen[-c(1,2)]
toRemove <- which(is.na(wintering_origen) == T)
if(length(toRemove)>0){
  wintering_origen <- wintering_origen[-toRemove]
  breeding_origen <- breeding_origen[-toRemove]
}
pts <- list()
pts$ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
pts$ptsNB <- as.matrix(wintering.assignments[,6:5][wintering_origen,])
pts$type <- rep("genetics", nrow(pts$ptsNB))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                as.matrix(cbind(rep(spp, length(pts$ptsBR[,1])),
                                                paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                                pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex)))


##  9 - Blue winged Warbler  ##
spp_number <- 9
spp <- as.character(species.names$ebird_code[spp_number])
# Tracking data
trackingData <- read.csv("resources/Tracking-data/vermNBlatWKbk.csv", header=T)
trackingData$depLON <- -trackingData$depLON
trackingData$nbLON <- -trackingData$nbLON
pts<- list()
pts$ptsBR <- as.matrix(trackingData[,3:2][which(trackingData$sp == 2),])
pts$ptsNB <- as.matrix(trackingData[,4:5][which(trackingData$sp == 2),])
pts$type <- c(rep("tracking", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  10 - Ovenbird  ##
spp_number <- 10
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/OVEN_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "ovenbird"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  11 - Willow Flycatcher  ##
spp_number <- 11
spp <- as.character(species.names$ebird_code[spp_number])
# Genetic data
wintering.assignments <- readRDS("resources/Genoscapes/WillowFlycatcher/winter-migrant-assignments.rds")
wintering.surfaces.origen <- readRDS("resources/Genoscapes/WillowFlycatcher/WIFL_WinterProbSurfaces.rds")
breeding_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  breeding_origen[k] <- which(wintering.surfaces.origen[,k] == max(wintering.surfaces.origen[,k]))
}
breeding_origen <- breeding_origen[-c(1,2)]
wintering_origen <- vector()
for(k in 3:length(colnames(wintering.surfaces.origen))){
  if(length(which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])) > 0){
    wintering_origen[k] <- which(wintering.assignments$indiv == colnames(wintering.surfaces.origen)[k])
  }else{
    wintering_origen[k] <- NA
  }
}
wintering_origen <- wintering_origen[-c(1,2)]
toRemove <- which(is.na(wintering_origen) == T)
if(length(toRemove)>0){
  wintering_origen <- wintering_origen[-toRemove]
  breeding_origen <- breeding_origen[-toRemove]
}
pts <- list()
pts$ptsBR = wintering.surfaces.origen[,1:2][breeding_origen,]
pts$ptsNB <- as.matrix(wintering.assignments[,13:12][wintering_origen,])
pts$type <- rep("genetics", nrow(pts$ptsNB))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                as.matrix(cbind(rep(spp, length(pts$ptsBR[,1])),
                                                paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                                pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex)))


## 12 - Grasshopper Sparrow  ##
spp_number <- 12
spp <- as.character(species.names$ebird_code[spp_number])
# Tracking data
trackingData <- read.csv("resources/Tracking-data/GrasshopperSparrow_data.csv", sep=";", header=T)
pts<- list()
pts$ptsBR <- as.matrix(trackingData[,1:2])
pts$ptsNB <- as.matrix(trackingData[,3:4])
pts$type <- c(rep("tracking", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))

##  13 - White-troated Sparrow  ##
spp_number <- 13
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/WTSP_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  14 - American Goldfinch  ##
spp_number <- 14
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/AMGO_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  15 - Purple Finch  ##
spp_number <- 15
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/PUFI_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  16 - Common Loon  ##
spp_number <- 16
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/COLO_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  17 - Brown Trasher  ##
spp_number <- 17
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/BRTH_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  18 - Gray Catbird  ##
spp_number <- 18
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/GRCA_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "gray catbird"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  19 - Brown-headed Cowbird  ##
spp_number <- 19
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/BHCO_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  20 - Common Grackle  ##
spp_number <- 20
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/COGR_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  21 - Tree Swallow  ##
spp_number <- 21
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/TRES_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  22 - Barn Swallow  ##
spp_number <- 22
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/BARS_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "barn swallow"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  23 - American Kestrel  ##
spp_number <- 23
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/AMKE_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  24 - Burrowing Owl  ##
spp_number <- 24
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/BUOW_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
pts$type <- c(rep("banding", nrow(pts$ptsNB)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))


##  25 - Osprey  ##
spp_number <- 25
spp <- as.character(species.names$ebird_code[spp_number])
# Banding data
bandingData <- read.csv("resources/Banding-data/OSPR_combine.csv", header=T)
pts <- banding_data_fct(bandingData)
# Tracking data
trackingData <- read.csv("resources/Tracking-data/data_Finch_etal_2017.csv", header=T)
trackingData <- trackingData[which(trackingData$species == "osprey"),]
pts$ptsBR <- rbind(pts$ptsBR, cbind(trackingData$blon, trackingData$blat))
pts$ptsNB <- rbind(pts$ptsNB, cbind(trackingData$wlon1, trackingData$wlat1))
pts$type <- c(rep("banding", nrow(bandingData)), rep("tracking", nrow(trackingData)))
##  Get hexagons containing the empirical breeding and wintering locations of individuals  ## 
pts <- empirical_hexagons_fct(pts, spp_number)
##  Add this species to the processed dataset
species.empirical.data <- rbind(species.empirical.data, 
                                cbind(rep(spp, length(pts$ptsBR[,1])),
                                      paste0(spp, "_", 1:length(pts$ptsBR[,1])),
                                      pts$type,pts$ptsBR, pts$ptsNB, pts$ptsBR_hex, pts$ptsNB_hex))

colnames(species.empirical.data) <- c("species", "individual", "data_type", "breeding_lon", "breeding_lat", "wintering_lon", "wintering_lat", "breeding_hex", "wintering_hex")

write.csv(species.empirical.data, "results/output/empirical_connectivity_data.csv", row.names = F)

save(hexgrid3, file = "results/output/hexagon_grid.RData")





