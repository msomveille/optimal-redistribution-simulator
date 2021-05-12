###  Plot figures for Somveille et al. A general theory of avian migratory connectivity  ###

library(geosphere)
library(ade4)
library(sp)
library(tidyverse)
library(maps)

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

##  Load migratory connectivity strenght
migratory_connectivity_strength <- read.csv("results/output/migratory_connectivity_strength.csv")

##  Load results of null models
load("results/output/null_models_results.RData")



##################
## Plot Figures ##
##################


##  FIGURE 2: ORSIM is capturing empirical migratory connectivity patterns  ##

png("results/figures/Figure_2.png", width=8, height=6.8, unit="in", res=300)
par(mfrow=c(3,3), mar=c(2.2,0.1,0.1,0.1), oma=c(0,0,0,0), mgp=c(1,0.6,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)")

spp_sel <- c(10,15,6)

for(j in 1:length(spp_sel)){
  
  empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[spp_sel[j]]),]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,spp_sel[j]] > 0),][empirical.data.spp$breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,spp_sel[j]] > 0),][empirical.data.spp$wintering_hex,]
  # Empirical migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,spp_sel[j]]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,spp_sel[j]]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,spp_sel[j]]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,spp_sel[j]]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_lon)){
    spl = SpatialLines(list(
      Lines(Line(
        rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
              c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
    plot(spl, add=T, col="black", lwd=1.5)
  }
  mtext(species.names$common_name[spp_sel[j]], side=1, line=-1.8, at=-140, cex=1)
  mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 1){title("Empirical connectivity", line=1.6, cex.main=1.5)}
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[spp_sel[j]]][empirical.data.spp$breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,spp_sel[j]]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Simulated migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,spp_sel[j]]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,spp_sel[j]]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,spp_sel[j]]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,spp_sel[j]]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    ORSIM_flows_sel <- ORSIM_flows[i,]
    if(sum(ORSIM_flows_sel)>0){
      for(k in 1:length(which(ORSIM_flows_sel > 0))){
        spl = SpatialLines(list(
          Lines(Line(
            rbind(hexgrid3_centroids[which(winter.abundance[,spp_sel[j]]>0),][which(ORSIM_flows_sel > 0)[k],],
                  coords_breeding_hex[i,] )),ID="a")))
        plot(spl, add=T, col="black", lwd=1.5, cex=2)
      }
    }	
  }
  mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 1){title("Simulated connectivity", line=1.6, cex.main=1.5)}
  
  # Simulated versus empirical
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distanceMat.simulated <- vector()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distanceMat.simulated <- distanceMat.simulated[-toRemove,]
    distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
  }
  plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,9500), ylim=c(0,9500), cex.lab=1.5)
  axis(side=1); axis(side=2)
  abline(a = 0, b = 1, col="red", lwd=2)
  mtext("Empirical", side=1, line=1.6, at=4750, cex=0.85)
  mtext("Simulated", side=2, line=1.6, at=4750, cex=0.85)
  mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.2)
  mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[spp_sel[j]],3))), side=1, line=-1.4, at=7700, cex=1.1)
  if(j == 1){title("Simulated vs. empirical", line=1.6, cex.main=1.5)}
}
dev.off()	




##  FIGURE 3: Energy efficiency drives the seasonal redistribution of individuals within migatory species  ##

# Compute the p-values (i.e. fraction of null simulations with a Mantel correlation higher than ORSIM) 
# and effect size (i.e. ORSIM - mean(null))
pval.nullModel1 <- pval.nullModel2 <- pval.nullModel3 <- vector()
pval.nullModel1.longitude <- pval.nullModel2.longitude <- pval.nullModel3.longitude <- vector()
pval.nullModel1.latitude <- pval.nullModel2.latitude <- pval.nullModel3.latitude <- vector()
es.nullModel1 <- es.nullModel2 <- es.nullModel3 <- vector()
es.nullModel1.longitude <- es.nullModel2.longitude <- es.nullModel3.longitude <- vector()
es.nullModel1.latitude <- es.nullModel2.latitude <- es.nullModel3.latitude <- vector()
for(i in 1:length(ORSIM_validation_metric$ORSIM_rM)){
  pval.nullModel1[i] <- length(which(nullModel1.rM[,i] > ORSIM_validation_metric$ORSIM_rM[i])) / length(nullModel1.rM[,i])
  pval.nullModel2[i] <- length(which(nullModel2.rM[,i] > ORSIM_validation_metric$ORSIM_rM[i])) / length(nullModel2.rM[,i])
  pval.nullModel3[i] <- length(which(nullModel3.rM[,i] > ORSIM_validation_metric$ORSIM_rM[i])) / length(nullModel3.rM[,i])
  es.nullModel1[i] <- ORSIM_validation_metric$ORSIM_rM[i] - mean(nullModel1.rM[,i])
  es.nullModel2[i] <- ORSIM_validation_metric$ORSIM_rM[i] - mean(nullModel2.rM[,i])
  es.nullModel3[i] <- ORSIM_validation_metric$ORSIM_rM[i] - mean(nullModel3.rM[,i])
  pval.nullModel1.longitude[i] <- length(which(nullModel1.rM.longitude[,i] > ORSIM_validation_metric$ORSIM_rM_longitude[i])) / length(nullModel1.rM.longitude[,i])
  pval.nullModel2.longitude[i] <- length(which(nullModel2.rM.longitude[,i] > ORSIM_validation_metric$ORSIM_rM_longitude[i])) / length(nullModel2.rM.longitude[,i])
  pval.nullModel3.longitude[i] <- length(which(nullModel3.rM.longitude[,i] > ORSIM_validation_metric$ORSIM_rM_longitude[i])) / length(nullModel3.rM.longitude[,i])
  es.nullModel1.longitude[i] <- ORSIM_validation_metric$ORSIM_rM_longitude[i] - mean(nullModel1.rM.longitude[,i])
  es.nullModel2.longitude[i] <- ORSIM_validation_metric$ORSIM_rM_longitude[i] - mean(nullModel2.rM.longitude[,i])
  es.nullModel3.longitude[i] <- ORSIM_validation_metric$ORSIM_rM_longitude[i] - mean(nullModel3.rM.longitude[,i])
  pval.nullModel1.latitude[i] <- length(which(nullModel1.rM.latitude[,i] > ORSIM_validation_metric$ORSIM_rM_latitude[i])) / length(nullModel1.rM.latitude[,i])
  pval.nullModel2.latitude[i] <- length(which(nullModel2.rM.latitude[,i] > ORSIM_validation_metric$ORSIM_rM_latitude[i])) / length(nullModel2.rM.latitude[,i])
  pval.nullModel3.latitude[i] <- length(which(nullModel3.rM.latitude[,i] > ORSIM_validation_metric$ORSIM_rM_latitude[i])) / length(nullModel3.rM.latitude[,i])
  es.nullModel1.latitude[i] <- ORSIM_validation_metric$ORSIM_rM_latitude[i] - mean(nullModel1.rM.latitude[,i])
  es.nullModel2.latitude[i] <- ORSIM_validation_metric$ORSIM_rM_latitude[i] - mean(nullModel2.rM.latitude[,i])
  es.nullModel3.latitude[i] <- ORSIM_validation_metric$ORSIM_rM_latitude[i] - mean(nullModel3.rM.latitude[,i])
}

# Format the results for plotting
results.null.models <- data.frame(
  species = ORSIM_validation_metric$species,
  ORSIM_rM = ORSIM_validation_metric$ORSIM_rM,
  pval_null1 = pval.nullModel1,
  es_null1 = es.nullModel1,
  pval_null2 = pval.nullModel2,
  es_null2 = es.nullModel2,
  pval_null3 = pval.nullModel3,
  es_null3 = es.nullModel3,
  MC_empirical = migratory_connectivity_strength$MC_empirical,
  MC_simulated = migratory_connectivity_strength$MC_simulated
)

# Plot the figure
pdf("results/figures/Figure_3.pdf", width=8, height=2.7, pointsize = 0.5)
theme_set(theme_classic())
p1 <- ggplot(data = results.null.models, aes(ORSIM_rM)) +
  geom_density(alpha = 1, fill = "dark grey", color = "black") +
  xlim(-1,1) +
  xlab(bquote('r'["M"])) +
  ylab("Number of species") +
  geom_text(x=0, y=3.3, label="(a)")
p2 <- ggplot(data = results.null.models) +
  geom_density(aes(es_null1, fill = "Null model 1"), alpha = 0.4) +
  geom_density(aes(es_null2, fill = "Null model 2"), alpha = 0.4) +
  geom_density(aes(es_null3, fill = "Null model 3"), alpha = 0.4) +
  xlim(-2,2) +
  xlab("Effect size") +
  ylab("Number of species") +
  theme(legend.position = c(0.255,0.75), legend.title = element_blank()) +
  scale_color_manual(values = c("Null model 1" = "yellow", "Null model 2" = "red", "Null model 3" = "blue")) +
  geom_text(x=0, y=1.61, label="(b)")
p3 <- ggplot(data = results.null.models, aes(x=MC_simulated, y=MC_empirical)) +
  geom_point() +
  geom_smooth(method = "lm", se=T, col="black") +
  xlab("MC simulated") +
  ylab("MC empirical") +
  geom_text(x=0.7, y=0.87, label="(c)")
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

print(paste0("Correlation between empirical and simulated MC: ", cor(results.null.models$MC_simulated, results.null.models$MC_empirical)))
print(paste0("P-value of the correlation between empirical and simulated MC: ", summary(lm(results.null.models$MC_empirical ~ results.null.models$MC_simulated))$coefficients[2,4]))

print(paste0("t value of the correlation between empirical MC and ORSIM rM: ", summary(lm(results.null.models$MC_empirical ~ results.null.models$ORSIM_rM))$coefficients[2,3]))
print(paste0("P-value of the correlation between empirical MC and ORSIM rM: ", summary(lm(results.null.models$MC_empirical ~ results.null.models$ORSIM_rM))$coefficients[2,4]))




##  FIGURE S1: species' relative abundance distribution  ##

png("results/figures/Figure_S1.png", width=8, height=9, unit="in", res=300)
par(mfrow=c(5,5), mar=c(2.2,2.2,1,1), mgp=c(1.3,0.5,0), bg="white")	
for(j in 1:length(species.names$ebird_code)){
  maps::map("world", fill=T, col="light grey", border="light grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -30), ylim=c(-40,70))
  rbPal.BR <- colorRampPalette(c("yellow2","red3"))
  rbPal.NB <- colorRampPalette(c("skyblue","dark blue"))
  datcol.BR <- rbPal.BR(6)[as.numeric(cut(summer.abundance[,j], breaks=seq(range(summer.abundance[,j], na.rm=T)[1], range(summer.abundance[,j], na.rm=T)[2], range(summer.abundance[,j], na.rm=T)[2]/6)))]
  datcol.BR[which(summer.abundance[,j] == range(summer.abundance[,j], na.rm=T)[2])] <- "red3"
  datcol.BR[which(summer.abundance[,j] == "NaN" | is.na(summer.abundance[,j]) == T | summer.abundance[,j] == 0)] <- "light grey"
  datcol.NB <- rbPal.NB(6)[as.numeric(cut(winter.abundance[,j], breaks=seq(range(winter.abundance[,j], na.rm=T)[1], range(winter.abundance[,j], na.rm=T)[2], range(winter.abundance[,j], na.rm=T)[2]/6)))]
  datcol.NB[which(winter.abundance[,j] == range(winter.abundance[,j], na.rm=T)[2])] <- "dark blue"
  datcol.NB[which(winter.abundance[,j] == "NaN" | is.na(winter.abundance[,j]) == T | winter.abundance[,j] == 0)] <- "light grey"
  datcol <- datcol.NB
  datcol[which(datcol.NB == "light grey" & datcol.BR != "light grey")] <- datcol.BR[which(datcol.NB == "light grey" & datcol.BR != "light grey")]
  plot(hexgrid3, col=datcol, border=datcol, add=T)
  legend("bottomleft", inset=0.06, bg="white", box.col="white", legend=round(rev(seq(range(summer.abundance[,j], na.rm=T)[1], range(summer.abundance[,j], na.rm=T)[2], range(summer.abundance[,j], na.rm=T)[2]/6))[1:6],2), fill=rev(rbPal.BR(6)), cex=0.55, border=rev(rbPal.BR(6)), title="Relative\nabundance")
  legend("bottomleft", inset=c(0.3,0.06), bg="white", box.col="white", legend=round(rev(seq(range(winter.abundance[,j], na.rm=T)[1], range(winter.abundance[,j], na.rm=T)[2], range(winter.abundance[,j], na.rm=T)[2]/6))[1:6],2), fill=rev(rbPal.NB(6)), cex=0.55, border=rev(rbPal.NB(6)))
  title(main=species.names$common_name[j], line=0.3)
}
dev.off()



##  FIGURE S5: Empirical versus simulated migratory connectivity for 5 species: Wood Thrush, Swainson's Thrush, Hermit Thrush, American Robin, and Ovenbird  ##

#pdf("Manuscript/Science/Figures/FigS2.pdf", width=12, height=20)
png("results/figures/Figure_S5.png", width=8, height=11.3, unit="in",res=300)
par(mfrow=c(5,3), mar=c(2.2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)")

for(j in 1:5){
  
  empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[j]),]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][empirical.data.spp$breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][empirical.data.spp$wintering_hex,]
  # Empirical migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_lon)){
    spl = SpatialLines(list(
      Lines(Line(
          rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
                c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
    plot(spl, add=T, col="black", lwd=1.5)
  }
  mtext(species.names$common_name[j], side=1, line=-3, at=-135, cex=0.9)
  mtext(panel[(j+(2*(j-1)))], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 1){title("Empirical connectivity", line=1.6, cex.main=1.5)}
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][empirical.data.spp$breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Simulated migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    ORSIM_flows_sel <- ORSIM_flows[i,]
    if(sum(ORSIM_flows_sel)>0){
      for(k in 1:length(which(ORSIM_flows_sel > 0))){
        spl = SpatialLines(list(
          Lines(Line(
            rbind(hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows_sel > 0)[k],],
                  coords_breeding_hex[i,] )),ID="a")))
        plot(spl, add=T, col="black", lwd=1.5, cex=2)
      }
    }	
  }
  mtext(panel[(j+(2*(j-1))+1)], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 1){title("Simulated connectivity", line=1.6, cex.main=1.5)}
  
  # Simulated versus empirical
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distanceMat.simulated <- vector()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distanceMat.simulated <- distanceMat.simulated[-toRemove,]
    distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
  }
  plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
  axis(side=1); axis(side=2)
  abline(a = 0, b = 1, col="red", lwd=2)
  mtext("Empirical", side=1, line=1.6, at=3500, cex=0.85)
  mtext("Simulated", side=2, line=1.6, at=3500, cex=0.85)
  mtext(panel[(j+(2*(j-1))+2)], side=3, line=0.5, at=-500, cex=1.2)
  mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[j],3))), side=1, line=-1.3, at=5700, cex=1.1)
  if(j == 1){title("Simulated vs. empirical", line=1.6, cex.main=1.5)}
}
dev.off()	



##  FIGURE S6: Empirical versus simulated migratory connectivity for 5 species: Yellow Warbler, Common Yellowthroat, Wilson's Warbler, Blue winged Warbler, and White-throated Sparrow  ##

png("results/figures/Figure_S6.png", width=8, height=11.3, unit="in",res=300)
par(mfrow=c(5,3), mar=c(2.2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)")

for(j in 6:10){
  
  empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[j]),]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][empirical.data.spp$breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][empirical.data.spp$wintering_hex,]
  # Empirical migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_lon)){
    spl = SpatialLines(list(
      Lines(Line(
        rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
              c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
    plot(spl, add=T, col="black", lwd=1.5)
  }
  mtext(species.names$common_name[j], side=1, line=-3, at=-135, cex=0.9)
  mtext(panel[(j-5+(2*(j-6)))], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 6){title("Empirical connectivity", line=1.6, cex.main=1.5)}
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][empirical.data.spp$breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Simulated migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    ORSIM_flows_sel <- ORSIM_flows[i,]
    if(sum(ORSIM_flows_sel)>0){
      for(k in 1:length(which(ORSIM_flows_sel > 0))){
        spl = SpatialLines(list(
          Lines(Line(
            rbind(hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows_sel > 0)[k],],
                  coords_breeding_hex[i,] )),ID="a")))
        plot(spl, add=T, col="black", lwd=1.5, cex=2)
      }
    }	
  }
  mtext(panel[(j-5+(2*(j-6))+1)], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 6){title("Simulated connectivity", line=1.6, cex.main=1.5)}
  
  # Simulated versus empirical
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distanceMat.simulated <- vector()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distanceMat.simulated <- distanceMat.simulated[-toRemove,]
    distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
  }
  plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
  axis(side=1); axis(side=2)
  abline(a = 0, b = 1, col="red", lwd=2)
  mtext("Empirical", side=1, line=1.6, at=3500, cex=0.85)
  mtext("Simulated", side=2, line=1.6, at=3500, cex=0.85)
  mtext(panel[(j-5+(2*(j-6))+2)], side=3, line=0.5, at=-500, cex=1.2)
  mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[j],3))), side=1, line=-1.3, at=5700, cex=1.1)
  if(j == 6){title("Simulated vs. empirical", line=1.6, cex.main=1.5)}
}
dev.off()	



##  FIGURE S7: Empirical versus simulated migratory connectivity for 5 species: Willow Flycatcher, Grasshopper Sparrow, White-throated Sparrow, American Goldfinch, and Purple Finch  ##

png("results/figures/Figure_S7.png", width=8, height=11.3, unit="in",res=300)
par(mfrow=c(5,3), mar=c(2.2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)")

for(j in 11:15){
  
  empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[j]),]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][empirical.data.spp$breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][empirical.data.spp$wintering_hex,]
  # Empirical migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_lon)){
    spl = SpatialLines(list(
      Lines(Line(
        rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
              c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
    plot(spl, add=T, col="black", lwd=1.5)
  }
  mtext(species.names$common_name[j], side=1, line=-3, at=-135, cex=0.9)
  mtext(panel[(j-10+(2*(j-11)))], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 11){title("Empirical connectivity", line=1.6, cex.main=1.5)}
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][empirical.data.spp$breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Simulated migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    ORSIM_flows_sel <- ORSIM_flows[i,]
    if(sum(ORSIM_flows_sel)>0){
      for(k in 1:length(which(ORSIM_flows_sel > 0))){
        spl = SpatialLines(list(
          Lines(Line(
            rbind(hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows_sel > 0)[k],],
                  coords_breeding_hex[i,] )),ID="a")))
        plot(spl, add=T, col="black", lwd=1.5, cex=2)
      }
    }	
  }
  mtext(panel[(j-10+(2*(j-11))+1)], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 11){title("Simulated connectivity", line=1.6, cex.main=1.5)}
  
  # Simulated versus empirical
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distanceMat.simulated <- vector()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distanceMat.simulated <- distanceMat.simulated[-toRemove,]
    distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
  }
  plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
  axis(side=1); axis(side=2)
  abline(a = 0, b = 1, col="red", lwd=2)
  mtext("Empirical", side=1, line=1.6, at=3500, cex=0.85)
  mtext("Simulated", side=2, line=1.6, at=3500, cex=0.85)
  mtext(panel[(j-10+(2*(j-11))+2)], side=3, line=0.5, at=-500, cex=1.2)
  mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[j],3))), side=1, line=-1.3, at=5700, cex=1.1)
  if(j == 11){title("Simulated vs. empirical", line=1.6, cex.main=1.5)}
}
dev.off()	



##  FIGURE S8: Empirical versus simulated migratory connectivity for 5 species: Common Loon, Brown Thrasher, Gray Catbird, Brown-headed Cowbird and Common Grackle  ##

png("results/figures/Figure_S8.png", width=8, height=11.3, unit="in", res=300)
par(mfrow=c(5,3), mar=c(2.2,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)")

for(j in 16:20){
  
  empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[j]),]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][empirical.data.spp$breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][empirical.data.spp$wintering_hex,]
  # Empirical migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_lon)){
    spl = SpatialLines(list(
      Lines(Line(
        rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
              c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
    plot(spl, add=T, col="black", lwd=1.5)
  }
  mtext(species.names$common_name[j], side=1, line=-3, at=-135, cex=0.9)
  mtext(panel[(j-15+(2*(j-16)))], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 16){title("Empirical connectivity", line=1.6, cex.main=1.5)}
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][empirical.data.spp$breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Simulated migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    ORSIM_flows_sel <- ORSIM_flows[i,]
    if(sum(ORSIM_flows_sel)>0){
      for(k in 1:length(which(ORSIM_flows_sel > 0))){
        spl = SpatialLines(list(
          Lines(Line(
            rbind(hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows_sel > 0)[k],],
                  coords_breeding_hex[i,] )),ID="a")))
        plot(spl, add=T, col="black", lwd=1.5, cex=2)
      }
    }	
  }
  mtext(panel[(j-15+(2*(j-16))+1)], side=3, line=0.5, at=-175, cex=1.2)
  if(j == 16){title("Simulated connectivity", line=1.6, cex.main=1.5)}
  
  # Simulated versus empirical
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distanceMat.simulated <- vector()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distanceMat.simulated <- distanceMat.simulated[-toRemove,]
    distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
  }
  plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
  axis(side=1); axis(side=2)
  abline(a = 0, b = 1, col="red", lwd=2)
  mtext("Empirical", side=1, line=1.6, at=3500, cex=0.85)
  mtext("Simulated", side=2, line=1.6, at=3500, cex=0.85)
  mtext(panel[(j-15+(2*(j-16))+2)], side=3, line=0.5, at=-500, cex=1.2)
  mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[j],3))), side=1, line=-1.3, at=5700, cex=1.1)
  if(j == 16){title("Simulated vs. empirical", line=1.6, cex.main=1.5)}
}
dev.off()	



##  FIGURE S9: Empirical versus simulated migratory connectivity for 5 species: Tree Swallow, Barn Swallow, American Kestrel, Burrowing Owl and Osprey  ##

png("results/figures/Figure_S9.png", width=8, height=12.8, unit="in", res=300)
par(mfrow=c(5,3), mar=c(6,0.1,0.1,0.1), mgp=c(1.9,0.7,0), bg="white")

panel <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)")

j = 21

empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[j]),]
coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][empirical.data.spp$breeding_hex,]
coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][empirical.data.spp$wintering_hex,]
# Empirical migratory connectivity
maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
plot(hexgrid3[which(winter.abundance[,j]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
for(i in 1:length(empirical.data.spp$breeding_lon)){
  spl = SpatialLines(list(
    Lines(Line(
      rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
            c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
  plot(spl, add=T, col="black", lwd=1.5)
}
mtext(species.names$common_name[j], side=1, line=-3, at=-135, cex=0.9)
mtext(panel[(j-20+(2*(j-21)))], side=3, line=0.5, at=-175, cex=1.2)
title("Empirical connectivity", line=1.6, cex.main=1.5)
# Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
ORSIM_flows <- ORSIM_results[[j]][empirical.data.spp$breeding_hex,]
winter.destinations <- list()
for(i in 1:length(empirical.data.spp$breeding_hex)){
  winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
  if(length(which(ORSIM_flows[i,] > 0)) > 1){
    winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
  }
}
# Simulated migratory connectivity
maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-15,70))
plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
plot(hexgrid3[which(winter.abundance[,j]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
for(i in 1:length(empirical.data.spp$breeding_hex)){
  ORSIM_flows_sel <- ORSIM_flows[i,]
  if(sum(ORSIM_flows_sel)>0){
    for(k in 1:length(which(ORSIM_flows_sel > 0))){
      spl = SpatialLines(list(
        Lines(Line(
          rbind(hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows_sel > 0)[k],],
                coords_breeding_hex[i,] )),ID="a")))
      plot(spl, add=T, col="black", lwd=1.5, cex=2)
    }
  }	
}
mtext(panel[(j-20+(2*(j-21))+1)], side=3, line=0.5, at=-175, cex=1.2)
title("Simulated connectivity", line=1.6, cex.main=1.5)

# Simulated versus empirical
# Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
distanceMat.empirical <- vector()
for(i in 1:nrow(coords_breeding_hex)){
  distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
}
# Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
distanceMat.simulated <- vector()
for(i in 1:length(empirical.data.spp$breeding_hex)){
  distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
}
toRemove <- which(lapply(winter.destinations, length) == 0)
if(length(toRemove)>0){
  distanceMat.simulated <- distanceMat.simulated[-toRemove,]
  distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
}
plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,7000), ylim=c(0,7000), cex.lab=1.5)
axis(side=1); axis(side=2)
abline(a = 0, b = 1, col="red", lwd=2)
mtext("Empirical", side=1, line=1.6, at=3500, cex=0.85)
mtext("Simulated", side=2, line=1.6, at=3500, cex=0.85)
mtext(panel[(j-20+(2*(j-21))+2)], side=3, line=0.5, at=-500, cex=1.2)
mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[j],3))), side=1, line=-1.3, at=5700, cex=1.1)
title("Simulated vs. empirical", line=1.6, cex.main=1.5)

for(j in 22:25){
  
  empirical.data.spp <- empirical.data[which(empirical.data$species == species.names$ebird_code[j]),]
  coords_breeding_hex <- hexgrid3_centroids[which(summer.abundance[,j] > 0),][empirical.data.spp$breeding_hex,]
  coords_wintering_hex <- hexgrid3_centroids[which(winter.abundance[,j] > 0),][empirical.data.spp$wintering_hex,]
  # Empirical migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-45,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][unique(empirical.data.spp$wintering_hex),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_lon)){
    spl = SpatialLines(list(
      Lines(Line(
        rbind(c(empirical.data.spp$breeding_lon[i], empirical.data.spp$breeding_lat[i]), 
              c(empirical.data.spp$wintering_lon[i], empirical.data.spp$wintering_lat[i]))), ID="a")), proj4string = CRS(proj4string(hexgrid3)))
    plot(spl, add=T, col="black", lwd=1.5)
  }
  mtext(species.names$common_name[j], side=1, line=-3, at=-135, cex=0.9)
  mtext(panel[(j-20+(2*(j-21)))], side=3, line=0.5, at=-175, cex=1.2)
  # Get wintering sites (hexagons) where simulated individuals from the empirical breeding sites (hexagons containing empirical breeding locations) are predicted to migrate
  ORSIM_flows <- ORSIM_results[[j]][empirical.data.spp$breeding_hex,]
  winter.destinations <- list()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    winter.destinations[[i]] <- hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows[i,] > 0),]
    if(length(which(ORSIM_flows[i,] > 0)) > 1){
      winter.destinations[[i]] <- apply(winter.destinations[[i]], 2, mean)
    }
  }
  # Simulated migratory connectivity
  maps::map("world", fill=T, col="grey", border="grey", bg="white", mar=c(0,0,0,0), xlim=c(-180, -20), ylim=c(-45,70))
  plot(hexgrid3[which(summer.abundance[,j]>0),], col="orange", border="orange", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),], col="light blue", border="light blue", add=T)
  plot(hexgrid3[which(summer.abundance[,j]>0),][unique(empirical.data.spp$breeding_hex),], col="red", border="red", add=T)
  plot(hexgrid3[which(winter.abundance[,j]>0),][which(apply(ORSIM_flows, 2, sum) > 0),], col="blue", border="blue", add=T)
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    ORSIM_flows_sel <- ORSIM_flows[i,]
    if(sum(ORSIM_flows_sel)>0){
      for(k in 1:length(which(ORSIM_flows_sel > 0))){
        spl = SpatialLines(list(
          Lines(Line(
            rbind(hexgrid3_centroids[which(winter.abundance[,j]>0),][which(ORSIM_flows_sel > 0)[k],],
                  coords_breeding_hex[i,] )),ID="a")))
        plot(spl, add=T, col="black", lwd=1.5, cex=2)
      }
    }	
  }
  mtext(panel[(j-20+(2*(j-21))+1)], side=3, line=0.5, at=-175, cex=1.2)

  # Simulated versus empirical
  # Compute pairwise distances between the sets of breeding and wintering locations of individuals in the empirical datasets
  distanceMat.empirical <- vector()
  for(i in 1:nrow(coords_breeding_hex)){
    distanceMat.empirical <- rbind(distanceMat.empirical, distHaversine(coords_breeding_hex[i,], coords_wintering_hex) / 1000)
  }
  # Compute pairwise distances between empirical breeding sites and corresponding simulated wintering sites
  distanceMat.simulated <- vector()
  for(i in 1:length(empirical.data.spp$breeding_hex)){
    distanceMat.simulated <- rbind(distanceMat.simulated, distHaversine(coords_breeding_hex[i,], matrix(unlist(winter.destinations), ncol=2, byrow=T)) / 1000)
  }
  toRemove <- which(lapply(winter.destinations, length) == 0)
  if(length(toRemove)>0){
    distanceMat.simulated <- distanceMat.simulated[-toRemove,]
    distanceMat.empirical <- distanceMat.empirical[-toRemove,-toRemove]
  }
  plot(as.vector(distanceMat.empirical), as.vector(distanceMat.simulated), pch=20, cex=0.8, axes=F, main=NULL, ylab="", xlab="", xlim=c(0,12000), ylim=c(0,12000), cex.lab=1.5)
  axis(side=1); axis(side=2)
  abline(a = 0, b = 1, col="red", lwd=2)
  mtext("Empirical", side=1, line=1.6, at=6000, cex=0.85)
  mtext("Simulated", side=2, line=1.6, at=6000, cex=0.85)
  mtext(panel[(j-20+(2*(j-21))+2)], side=3, line=0.5, at=-500, cex=1.2)
  mtext(bquote('r'["M"]*' = '* .(round(ORSIM_validation_metric$ORSIM_rM[j],3))), side=1, line=-1.3, at=9700, cex=1.1)
}
dev.off()	


