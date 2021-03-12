


##  Construct an hexagon grid covering the Americas  ##

hexgrid <- dgconstruct(projection="ISEA", topology="HEXAGON", res=7, metric=T)
hexgrid_center <- dgSEQNUM_to_GEO(hexgrid, 1:21872) # 65612 / 590492
hexgrid_centroids <- cbind(hexgrid_center$lon_deg, hexgrid_center$lat_deg)
hex_sel <- which(hexgrid_centroids[,1] < -30 & hexgrid_centroids[,1] > -170 & hexgrid_centroids[,2] > -60 & hexgrid_centroids[,2] < 80)
hexgrid2_stem <- dgcellstogrid(hexgrid, hex_sel, frame=F,wrapcells=TRUE)
hexgrid2_stem_centroids <- matrix(unlist(lapply(hexgrid2_stem@polygons, function(x) x@labpt)), byrow=T, ncol=2)
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, proj4string(hexgrid2_stem))
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- gBuffer(newmap, byid=TRUE, width=0)
newmap <- gUnaryUnion(newmap, id=newmap@data$world)
hexgrid3_stem <- gIntersection(hexgrid2_stem, newmap, byid=T)
hexgrid3_stem_centroids <- matrix(unlist(lapply(hexgrid3_stem@polygons, function(x) x@labpt)), byrow=T, ncol=2)
sr <- proj4string(hexgrid3_stem)
