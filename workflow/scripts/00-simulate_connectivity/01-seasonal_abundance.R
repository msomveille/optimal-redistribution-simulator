###  Download and process seasonal abundance estimates from eBird based on spatio-tempororal exploratory models  ###

library(ebirdst)
library(dggridR)
library(rworldmap)
library(rgdal)
library(rgeos)
library(raster)
library(geosphere)

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

##  Lists of species names for the various datasets  ##
species.names <- read.csv("resources/species_names.csv") 

##  Load and process the rasters of seasonal abundances  ##
seasonalAbundance_BR <- matrix(nrow = length(hexgrid3), ncol = length(species.names$ebird_code))
seasonalAbundance_NB <- matrix(nrow = length(hexgrid3), ncol = length(species.names$ebird_code))
dir.create("results/output/seas_abund_rasters")
for(i in 1:length(species.names$ebird_code)){
  sp_path <- ebirdst_download(as.character(species.names$ebird_code[i]),
                              path = "results/output/seas_abund_rasters")
  abunds <- load_raster("abundance_seasonal", path = sp_path)
  # Project the hexagonal grid using the same projection than the eBird STEM rasters
  if(i==1){
    srr <- proj4string(abunds[[1]])
    hexgrid3 <- spTransform(hexgrid3, srr)
  }
  # Extract breeding abundance onto the hexagon grid
  abunds_BR <- abunds[[1]]
  seasonalAbundance_BR[,i] <- raster::extract(abunds_BR, hexgrid3, fun=mean, na.rm=T)
  # Extract wintering abundance onto the hexagon grid
  abunds_NB <- abunds[[2]]
  seasonalAbundance_NB[,i] <- raster::extract(abunds_NB, hexgrid3, fun=mean, na.rm=T)
  unlink(paste0("results/output/seas_abund_rasters/", list.files("results/output/seas_abund_rasters")), recursive=TRUE)
}
unlink("results/output/seas_abund_rasters", recursive=TRUE)
seasonalAbundance_BR[which(is.na(seasonalAbundance_BR)==T, arr.ind=T)] <- 0
seasonalAbundance_NB[which(is.na(seasonalAbundance_NB)==T, arr.ind=T)] <- 0
colnames(seasonalAbundance_BR) <- colnames(seasonalAbundance_NB) <- species.names$ebird_code
write.csv(seasonalAbundance_BR, "results/output/seasonalAbundance_BR.csv", row.names=F)
write.csv(seasonalAbundance_NB, "results/output/seasonalAbundance_NB.csv", row.names=F)  


##  Compute matrix of pairwise distance between every pair of hexagons on the grid  ##
distance.matrix <- vector()
for(i in 1:length(hexgrid3_centroids[,1])){	
  distance.matrix <- rbind(distance.matrix, distHaversine(hexgrid3_centroids[i,] , hexgrid3_centroids)/1000)
}
write.table(distance.matrix, "results/output/distanceMatrix.csv", row.names=F, col.names=F, sep=";")



