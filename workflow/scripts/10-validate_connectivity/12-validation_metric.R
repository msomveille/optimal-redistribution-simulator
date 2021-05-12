###  Calculate the validation metric for ORSIM  ###

library(geosphere)
library(ade4)
library(sp)

### Get necessary data

##  Get hexagon grid covering the Americas  ##
load("results/output/hexagon_grid.RData")
hexgrid3_centroids <- matrix(unlist(lapply(hexgrid3@polygons, function(x) x@labpt)), byrow=T, ncol=2)
sr <- proj4string(hexgrid3)

##  Lists of species names for the various datasets  ##
species.names <- read.csv("resources/species_names.csv")

##  Seasonal abundance estimates  ##
summer.abundance <- read.csv("results/output/seasonalAbundance_BR.csv", header=T)
winter.abundance <- read.csv("results/output/seasonalAbundance_NB.csv", header=T)

## Get empirical connectivity data
empirical.data <- read.csv("results/output/empirical_connectivity_data.csv")

##  Load ORSIM output (simulated migratory connectivity)  ##
ORSIM_results <- list()
for(j in 1:length(species.names$ebird_code)){
  ORSIM_results[[j]] <- read.csv(paste("results/output/ORSIM_results_", species.names$ebird_code[j], ".csv", sep=""), header=FALSE)
}
names(ORSIM_results) <- species.names$ebird_code
  


### Compute validation metric

ORSIM.rM <- vector()
for(j in 1:length(species.names$ebird_code)){
  breeding_hex <- empirical.data$breeding_hex[which(empirical.data$species == species.names$ebird_code[j])]
  wintering_hex <- empirical.data$wintering_hex[which(empirical.data$species == species.names$ebird_code[j])]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][wintering_hex,]
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distmat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distmat.empirical <- rbind(distmat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distmat.simulated <- vector()
  for(i in 1:length(breeding_hex)){
    distmat.simulated <- rbind(distmat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distmat.simulated <- distmat.simulated[-toRemove,]
    distmat.empirical <- distmat.empirical[-toRemove,-toRemove]
  }
  # Compute the Mantel correlation coefficients
  ORSIM.rM[j] <- mantel.rtest(as.dist(distmat.simulated, upper=T), as.dist(distmat.empirical, upper=T), nrepet=0)
  print(j)
}



### Compute validation metric using longitude only

ORSIM.rM.longitude <- vector()
for(j in 1:length(species.names$ebird_code)){
  breeding_hex <- empirical.data$breeding_hex[which(empirical.data$species == species.names$ebird_code[j])]
  wintering_hex <- empirical.data$wintering_hex[which(empirical.data$species == species.names$ebird_code[j])]
  longitude_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),1][breeding_hex]
  longitude_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),1][wintering_hex]
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distmat.empirical <- vector()
  for(i in 1:length(longitude_breeding_hex)){
    distmat.empirical <- rbind(distmat.empirical, longitude_breeding_hex[i] - longitude_wintering_hex)
  }
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),1]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- mean(winter.destinations[[i]])
    }
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distmat.simulated <- vector()
  for(i in 1:length(breeding_hex)){
    distmat.simulated <- rbind(distmat.simulated, longitude_breeding_hex[i] - unlist(winter.destinations))
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distmat.simulated <- distmat.simulated[-toRemove,]
    distmat.empirical <- distmat.empirical[-toRemove,-toRemove]
  }
  # Compute the Mantel correlation coefficients
  ORSIM.rM.longitude[j] <- mantel.rtest(as.dist(distmat.simulated, upper=T), as.dist(distmat.empirical, upper=T), nrepet=0)
  print(j)
}


### Compute validation metric using longitude only

ORSIM.rM.latitude <- vector()
for(j in 1:length(species.names$ebird_code)){
  breeding_hex <- empirical.data$breeding_hex[which(empirical.data$species == species.names$ebird_code[j])]
  wintering_hex <- empirical.data$wintering_hex[which(empirical.data$species == species.names$ebird_code[j])]
  latitude_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),2][breeding_hex]
  latitude_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),2][wintering_hex]
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distmat.empirical <- vector()
  for(i in 1:length(latitude_breeding_hex)){
    distmat.empirical <- rbind(distmat.empirical, latitude_breeding_hex[i] - latitude_wintering_hex)
  }
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),2]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- mean(winter.destinations[[i]])
    }
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distmat.simulated <- vector()
  for(i in 1:length(breeding_hex)){
    distmat.simulated <- rbind(distmat.simulated, latitude_breeding_hex[i] - unlist(winter.destinations))
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distmat.simulated <- distmat.simulated[-toRemove,]
    distmat.empirical <- distmat.empirical[-toRemove,-toRemove]
  }
  # Compute the Mantel correlation coefficients
  ORSIM.rM.latitude[j] <- mantel.rtest(as.dist(distmat.simulated, upper=T), as.dist(distmat.empirical, upper=T), nrepet=0)
  print(j)
}


ORSIM_validation_metric <- cbind(as.character(species.names$ebird_code), ORSIM.rM, ORSIM.rM.longitude, ORSIM.rM.latitude)
colnames(ORSIM_validation_metric) <- c("species", "ORSIM_rM", "ORSIM_rM_longitude", "ORSIM_rM_latitude")

write.csv(ORSIM_validation_metric, "results/output/ORSIM_validation_metric.csv", row.names=F)
