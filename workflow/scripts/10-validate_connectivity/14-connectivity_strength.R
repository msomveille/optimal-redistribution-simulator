###  Calculate the empirical and simulated strength of migratory connectivity  ####

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


## Calculate strength of migratory connectivity

MC.empirical <- MC.simulated <- vector()
for(j in 1:length(species.names$ebird_code)){
  breeding_hex <- empirical.data$breeding_hex[which(empirical.data$species == species.names$ebird_code[j])]
  wintering_hex <- empirical.data$wintering_hex[which(empirical.data$species == species.names$ebird_code[j])]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][wintering_hex,]
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove) > 0){
    winter.destinations <- winter.destinations[-toRemove]
    coords_breeding_hex <- coords_breeding_hex[-toRemove,]
    coords_wintering_hex <- coords_wintering_hex[-toRemove,]
  }
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical.BR <- distanceMat.empirical.NB <- distanceMat.simulated.NB <- vector() 
  for(i in 1:length(winter.destinations)){	
    distanceMat.empirical.BR <- rbind(distanceMat.empirical.BR, distHaversine( coords_breeding_hex[i,], coords_breeding_hex ) / 1000)
    distanceMat.empirical.NB <- rbind(distanceMat.empirical.NB, distHaversine( coords_wintering_hex[i,], coords_wintering_hex ) / 1000)
    distanceMat.simulated.NB <- rbind(distanceMat.simulated.NB, distHaversine( matrix(unlist(winter.destinations), ncol=2, byrow=T)[i,] , matrix(unlist(winter.destinations), ncol=2, byrow=T) ) / 1000)
  }
  MC.empirical[j] <- mantel.rtest(as.dist(distanceMat.empirical.BR, upper=T), as.dist(distanceMat.empirical.NB, upper=T), nrepet=0)
  MC.simulated[j] <- mantel.rtest(as.dist(distanceMat.empirical.BR, upper=T), as.dist(distanceMat.simulated.NB, upper=T), nrepet=0)
}

migratory_connectivity_strength <- cbind(as.character(species.names$ebird_code), MC.empirical, MC.simulated)
colnames(migratory_connectivity_strength) <- c("species", "MC_empirical", "MC_simulated")

write.csv(migratory_connectivity_strength, "results/output/migratory_connectivity_strength.csv", row.names=F)
