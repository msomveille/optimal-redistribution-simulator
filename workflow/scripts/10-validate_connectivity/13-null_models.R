###  Compare ORSIM performance to null models  ###

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

##  Load ORSIM validation data
ORSIM_validation_metric <- read.csv("results/output/ORSIM_validation_metric.csv")

##  Compute the null models for each species  ##	

nullModel1.rM <- nullModel2.rM <- nullModel3.rM <- matrix(nrow=1000, ncol=nrow(species.names))
nullModel1.rM.longitude <- nullModel2.rM.longitude <- nullModel3.rM.longitude <- matrix(nrow=1000, ncol=nrow(species.names))
nullModel1.rM.latitude <- nullModel2.rM.latitude <- nullModel3.rM.latitude <- matrix(nrow=1000, ncol=nrow(species.names))
for(j in 1:length(species.names$ebird_code)){
  breeding_hex <- empirical.data$breeding_hex[which(empirical.data$species == species.names$ebird_code[j])]
  wintering_hex <- empirical.data$wintering_hex[which(empirical.data$species == species.names$ebird_code[j])]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][wintering_hex,]
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- distanceMat.empirical.longitude <- distanceMat.empirical.latitude <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
    distanceMat.empirical.longitude <- rbind(distanceMat.empirical.longitude, coords_breeding_hex[i,1] - coords_wintering_hex[,1])
    distanceMat.empirical.latitude <- rbind(distanceMat.empirical.latitude, coords_breeding_hex[i,2] - coords_wintering_hex[,2])
  }
  ##  Compute pariwise distance between empirical breeding sites and all wintering sites 	##	
  distanceMat.simulated.all <- distanceMat.simulated.all.longitude <- distanceMat.simulated.all.latitude <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.simulated.all <- rbind(distanceMat.simulated.all, distHaversine( coords_breeding_hex[i,], hexgrid3_centroids[which(winter.abundance[,j]>0),] ) / 1000)
    distanceMat.simulated.all.longitude <- rbind(distanceMat.simulated.all.longitude, coords_breeding_hex[i,1] - hexgrid3_centroids[which(winter.abundance[,j]>0),1] )
    distanceMat.simulated.all.latitude <- rbind(distanceMat.simulated.all.latitude, coords_breeding_hex[i,2] - hexgrid3_centroids[which(winter.abundance[,j]>0),2] )
  }

  ##  Run null models
  for(k in 1:1000){
  
    ##  Null model 1: randomizing migratory destinations  ##
    # Sampling winter destinations at random across the wintering ground
    winter.destinations.null <- sample(which(winter.abundance[,j]>0), length(breeding_hex), replace=T)
    # Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
    distanceMat.null <- distanceMat.null.longitude <- distanceMat.null.latitude <- vector()
    for(i in 1:length(breeding_hex)){
      distanceMat.null <- rbind(distanceMat.null, distHaversine(coords_breeding_hex[i,], hexgrid3_centroids[winter.destinations.null,] ) / 1000)
      distanceMat.null.longitude <- rbind(distanceMat.null.longitude, coords_breeding_hex[i,1] - hexgrid3_centroids[winter.destinations.null,1] )
      distanceMat.null.latitude <- rbind(distanceMat.null.latitude, coords_breeding_hex[i,2] - hexgrid3_centroids[winter.destinations.null,2] )
    }
    # Compute the Mantel correlation coefficients
    nullModel1.rM[k,j] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
    nullModel1.rM.longitude[k,j] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
    nullModel1.rM.latitude[k,j] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
  
    
    ##  Null model 2: selecting migratory destinations solely based on minimizing migration cost  ##
    # Sampling winter destinations across the wintering ground with a probability that decreases with migration distance
    winter.destinations.null <- apply(distanceMat.simulated.all, 1, function(x) which(x == sample(x,1,prob=(1/exp(0.01*x)))))
    # Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
    distanceMat.null <- distanceMat.null.longitude <- distanceMat.null.latitude <- vector()
    for(i in 1:length(breeding_hex)){
      distanceMat.null <- rbind(distanceMat.null, distHaversine(coords_breeding_hex[i,], hexgrid3_centroids[winter.destinations.null,] ) / 1000)
      distanceMat.null.longitude <- rbind(distanceMat.null.longitude, coords_breeding_hex[i,1] - hexgrid3_centroids[winter.destinations.null,1] )
      distanceMat.null.latitude <- rbind(distanceMat.null.latitude, coords_breeding_hex[i,2] - hexgrid3_centroids[winter.destinations.null,2] )
    }
    # Compute the Mantel correlation coefficients
    nullModel2.rM[k,j] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
    nullModel2.rM.longitude[k,j] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
    nullModel2.rM.latitude[k,j] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
  
    
    ##  Null model 3: selecting migratory destinations solely based on maximizing energy acquisition  ##
    # Sampling winter destinations across the wintering ground with a probability that increases with relative abundance
    winter.destinations.null <- sample(which(winter.abundance[,j]>0), length(breeding_hex), replace=T, prob=winter.abundance[,j][which(winter.abundance[,j]>0)])
    # Compute pairwise distances between empirical breeding sites and corresponding sampled wintering sites
    distanceMat.null <- distanceMat.null.longitude <- distanceMat.null.latitude <- vector()
    for(i in 1:length(breeding_hex)){
      distanceMat.null <- rbind(distanceMat.null, distHaversine(coords_breeding_hex[i,], hexgrid3_centroids[winter.destinations.null,] ) / 1000)
      distanceMat.null.longitude <- rbind(distanceMat.null.longitude, coords_breeding_hex[i,1] - hexgrid3_centroids[winter.destinations.null,1] )
      distanceMat.null.latitude <- rbind(distanceMat.null.latitude, coords_breeding_hex[i,2] - hexgrid3_centroids[winter.destinations.null,2] )
    }
    # Compute the Mantel correlation coefficients
    nullModel3.rM[k,j] <- mantel.rtest(as.dist(distanceMat.null, upper=T), as.dist(distanceMat.empirical, upper=T), nrepet=0)
    nullModel3.rM.longitude[k,j] <- mantel.rtest(as.dist(distanceMat.null.longitude, upper=T), as.dist(distanceMat.empirical.longitude, upper=T), nrepet=0)
    nullModel3.rM.latitude[k,j] <- mantel.rtest(as.dist(distanceMat.null.latitude, upper=T), as.dist(distanceMat.empirical.latitude, upper=T), nrepet=0)
  }
  print(paste0("Null models for ", species.names$common_name[j], " done (", j, "/25)"))
}

save(list = c("nullModel1.rM", "nullModel1.rM.longitude", "nullModel1.rM.latitude", "nullModel2.rM", "nullModel2.rM.longitude", "nullModel2.rM.latitude", "nullModel3.rM", "nullModel3.rM.longitude", "nullModel3.rM.latitude"),
     file = "results/output/null_models_results.RData")


