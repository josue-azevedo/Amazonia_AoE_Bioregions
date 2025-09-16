# Run Rcluster"

setwd("~/Biogeografia_Bioregions_2021/Bioregions_All_species_July21")

library(recluster)
library(phytools)
library(geiger) 
library(dendextend)
library(raster)
library(sp)
library(maptools)
library(spdep)
library(prabclus)
library(geosphere)
library(speciesgeocodeR)
library(biomapME)
library(rgdal)
library(viridis)
library(maps)
library(betapart)
library(RColorBrewer)
library(plyr)
library(devtools)
library(recluster)
library(dendextend)

ocorrencias=read.table("ALL_July_21.txt", header = T)
head(ocorrencias)
colnames(ocorrencias)= c("species", "decimallongitude", "decimallatitude")
head(ocorrencias)
unique(ocorrencias$species)
length(unique(ocorrencias$species))
ocorrencias_p_plot = ocorrencias
coordinates(ocorrencias_p_plot) = ocorrencias_p_plot[,c("decimallongitude", "decimallatitude")]
proj4string(ocorrencias_p_plot) = "+proj=longlat +datum=WGS84 +no_defs"

mapa_amazonia = rgdal::readOGR("Amazonia.shp")
proj4string(mapa_amazonia)
grid_amazonia=rgdal::readOGR("grid_amazonia.shp")
proj4string(grid_amazonia) = proj4string(mapa_amazonia)
grid_amazonia_raster = raster("grid_amazonia.tif")

par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1) ##voltar ao parametros originais do plot-visual
plot(mapa_amazonia)
plot(ocorrencias_p_plot, pch=16, add = T) #pch simbolos
plot(grid_amazonia, border="red", add = T)

trueCentroids = centroid(grid_amazonia)
points(trueCentroids, col = "blue", pch=16, cex=0.2)
head(ocorrencias)
str(trueCentroids)
trueCentroids <- as.data.frame(trueCentroids)
str(trueCentroids)
trueCentroids$species <- rep("zzzzz",nrow(trueCentroids))
head(trueCentroids)
head(ocorrencias)
trueCentroids<- trueCentroids[,c(3,1,2)]
colnames(trueCentroids) <- colnames(ocorrencias)
head(trueCentroids)
pontos_and_grids_centroides <- rbind(ocorrencias,trueCentroids)

head(pontos_and_grids_centroides)
tail(pontos_and_grids_centroides)
nrow(pontos_and_grids_centroides)

# Default plot settings
par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)

#Use recluster.region procedure. Abaixo pode usar tanto os dados brutos quanto checar pela influencia do sampling gap usando os polígonos
#pres_abs6 = pres_abs5[,which(colSums(pres_abs5)>2)]
#dissim_matrix = pres_abs6 # raw data

# Diss matrix from polygonos 
outp_all <- SpGeoCod(ocorrencias, grid_amazonia,areanames="layer")
spp_ranges = CalcRange(outp_all,method = "pseudospherical", terrestrial = TRUE, rare = "buffer", buffer.width = 100000)
# calc presence absence
crossing<- raster::intersect(spp_ranges,grid_amazonia)
crossing_dat<-crossing[,c("species","layer")]
crossing_dat<-na.exclude(crossing_dat@data)
head(crossing_dat)
grid_coord<- centroid(grid_amazonia)
grid_coord<- as.data.frame(grid_coord)
grid_coord$layer<-grid_amazonia@data$layer
colnames(grid_coord) = c("x","y","layer")
head(grid_coord)
crossing_dat$long <- grid_coord$x[match(crossing_dat$layer,grid_coord$layer)]
crossing_dat$lat <- grid_coord$y[match(crossing_dat$layer,grid_coord$layer)]
tail(crossing_dat)
crossing_dat = crossing_dat[,c("species","long","lat")]
colnames(crossing_dat) <- colnames(ocorrencias)
tail(crossing_dat)
grid_coord$species <- rep("zzzzz",nrow(grid_coord))
head(grid_coord)
head(ocorrencias)
grid_coord<- grid_coord[,c("species","x","y")]
head(grid_coord)
colnames(grid_coord) <- colnames(ocorrencias)
head(grid_coord)
pontos_and_grids_centroides_pol <- rbind(crossing_dat,grid_coord)
SpGeoCod_object_pol = SpGeoCod(pontos_and_grids_centroides_pol, grid_amazonia, areanames = "layer")
pres_abs_p = SpGeoCod_object_pol$spec_table
pres_abs2_p <-pres_abs_p[,order(as.numeric(colnames(pres_abs_p)),decreasing = FALSE ) ]
pres_abs3_p = pres_abs2_p[!rownames(pres_abs2_p) %in% "zzzzz",]
pres_abs4_p = pres_abs3_p[,!colnames(pres_abs3_p) %in% "not_classified"]
pres_abs5_p = pres_abs4_p
pres_abs5_p<-1*(pres_abs5_p>0)
pres_abs6_p = pres_abs5_p[,which(colSums(pres_abs5_p)>1)]
dissim_matrix = pres_abs6_p # polygons 

#write.csv(manuscript_table, "manuscript_table_node_age.csv")

dissim_matrix = read.csv("dissim_matrix.csv", row.names = 1, header = TRUE)

colnames(dissim_matrix) = as.numeric(gsub("X","", colnames(dissim_matrix)))

dissim_matrix[1:5,1:5]

# Run Rcluster
library(recluster)
turn_cl<-recluster.region(t(dissim_matrix),dist="simpson",tr=50,rettree=TRUE,mincl=22,maxcl=22,method="average")
#saveRDS(turn_cl, "turn_cl_All_species.rds") #comando para salvar o arquivo turn gerado, substituir os nomes de acordo com o banco de dados usado
turn_cl=readRDS("turn_cl_All_species.rds") #para abrir o arquivo salvo

# Setting Best K
plot(turn_cl$solutions[,1],turn_cl$solutions[,3], xlab="Number of clusters", ylab="Explained dissimilarity")
abline(h=.95,col="red")

plot(turn_cl$solutions[1:30,1],turn_cl$solutions[1:30,3], xlab="Number of clusters", ylab="Explained dissimilarity")
abline(h=.95,col="red")

turn_cl$solutions
turn_cl$grouping
# Set best K according to figure
best_k = 22

#Select solution with three cluster and plot the tree.
plot(turn_cl$tree[[best_k]]) # Horrible, but we will change it
turn_cl$grouping[,best_k-1]














## Preparing a UPGMA dendrogram for plotting ####
dend <- as.dendrogram(turn_cl$tree[[best_k]])
## Getting cluster membership classification
cluster_membership = as.factor(turn_cl$grouping[,best_k-1]) # -1 because column name starts in 2

# Check cluster membership and colnames os diss_m
names(cluster_membership) %in% colnames(dissim_matrix)
names(cluster_membership) == colnames(dissim_matrix)

# Change cluster membership label order according to the dend label order
cluster_membership2=cluster_membership[order.dendrogram(dend)]
head(cluster_membership2)

# Set the names of dendrogram labels to be iqual to cluster_membership2
labels(dend) <- names(cluster_membership2)

# Checking the order of the colnames ###
names(cluster_membership2)
labels(dend)

# Check if cluster_membership is in the same order as dend labels
names(cluster_membership2) == labels(dend) # Not in the same order

# Here is to match again information in dend and cluster_membership2
dend2=branches_attr_by_clusters(dend,cluster_membership2)
names(cluster_membership2) == labels(dend2)

# This is to visualize to which branches the clusters are assigned (needed only when preparing the script)
#dend2 %>% set("branches_lwd",2.5) %>% ladderize(FALSE) %>%  plot(horiz  = TRUE,labels=FALSE,axes=FALSE) # still very ugly, lets improve a bit

# Now, note that when transformed to phylo, the tip.labels are in an increasing order as in cluster_membership
tree = ladderize(as.phylo(dend2))

tree$tip.label == names(cluster_membership)

## Grid names (or numbers)
tip.label<- names(cluster_membership[!duplicated(cluster_membership)])

# Cluster names (or numbers) #aqui de da o nome de cada clado, se manter somente n?meros, deleta "Bioregion", indo direto para unique
clade.label<- paste0("Bioregion ", unique(cluster_membership))

# N of grids per cluster membership
N = rep(10,length(unique(cluster_membership)))

# If you wanted equal sized triangles in plot (see what I mean by triangle later when plotting), set N as below:

#N = rep(10,length(unique(cluster_membership)))

# Now create a subtree from the tip.label names
sub_tree = ladderize(keep.tip(tree,tip.label))
plot(sub_tree)

# Reset tip labels of the subtree
tip.label = sub_tree$tip.label

## set crown node depth to 1/2 the maximum depth
depth<-sapply(tip.label,function(x,y) 0.5*y$edge.length[which(y$edge[,2]==which(y$tip.label== x))],y=sub_tree)
trans<-data.frame(tip.label,clade.label,N,depth)
rownames(trans)<-NULL
rm(tip.label,clade.label,N,depth)
trans
tt<-phylo.toBackbone(sub_tree,trans)

## plot
plot(tt)

# You probably can't see it, clean up the graphics device
dev.off()
plot(tt)

# If you still can't see the names, change the cex parameter
# Triangle width corresponds to the number of cells in a cluster
plot(tt, cex=.5)

# Or change other plotting parameters
par(fig = c(0, 1, 0, 1)) # change each number per time to see the best adjust
plot(tt, cex=.5)
#
par(fig = c(0.2, 1, 0, 1)) # change each number per time to see the best adjust
plot(tt, cex=.5)
#
par(fig = c(0.8, 1, 0, 1)) # change each number per time to see the best adjust
plot(tt, cex=.3)

# save plot as pdf for modifying colors in another program

# Plot dendrogram in colors # But needs to set similar colors together
plot(tt,col=viridis(best_k), cex=.3) #cex=.5 ? o tamanho da fonte

# Putting similar colors together
# Getting the order of the tips in the three
sub_tree$tip.label = trans$clade.label

is_tip <- sub_tree$edge[,2] <= length(sub_tree$tip.label)
ordered_tips <- sub_tree$edge[is_tip, 2]
order_in_tree = sub_tree$tip.label[ordered_tips]
order_in_tree

# Ordering the trans dataframe to add vector of colors
trans2 = trans
trans3 = trans2[order(match(trans2$clade.label,order_in_tree)),]

# Now add the colors of the palette (ex, viridis below)
trans3$col = viridis(best_k) # or create a vector of colors, eg, c("red","lightred","ect...")

# Or a combination (if 12,16,20,24,28,32,36... clusters)
#trans3$col = c(brewer.pal(best_k/4, "Blues")[-1],brewer.pal(best_k/4, "Greens")[-1],brewer.pal(best_k/4, "Greys")[-1],brewer.pal(best_k/4, "Oranges")[-1],brewer.pal(best_k/4, "Purples")[-1])

# Here using brewer.pal divergent colours (max 11 clusters)
#trans3$col = brewer.pal(best_k, "RdYlBu")

# Back to the original order
trans4 = trans3[order(match(trans3$clade.label, trans$clade.label)),]

# Now the plot color order was changed!
plot(tt,col=trans4$col, cex=.4)

# Below we'll use the same colors for the grid cells
# Passa para o grid
classified_areas = grid_amazonia[grid_amazonia$layer %in% names(turn_cl$grouping[,best_k-1]),]
classified_areas$class = turn_cl$grouping[,best_k-1]

classified_areas$col =  plyr::mapvalues(classified_areas$class, from=unique(turn_cl$grouping[,best_k-1]), to= trans4$col) 

# Change plot settings
# Default plot settings
# Default plot settings
par(mfrow=c(1,1))
par(mar=c(1,1,5,5)+0.1)
#plot(classified_areas) # Mude os números acima se quiser mudar a posição do plot. obs. teste um número de cada vez!

#para salvar como TIFF
#tiff("Clados_mapa_All.tiff", units="in", width=5, height=5, res=300) ###salvar imagem em tiff


plot(grid_amazonia_raster, col="white", legend=FALSE,axes=FALSE, box=FALSE)
plot(classified_areas, col = classified_areas$col, add=TRUE)
plot(grid_amazonia,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Bioregions rcluster"))
text(rasterize(classified_areas,grid_amazonia_raster,"class"), col="white",cex=.6)

# Marcado com # abaixo foi apenas para teste
#plot(classified_areas)
#c(0, 1, 0, 1)
# Não usei os parametros abaixo
#(x1, y1) = (0, 0) lower-left corner
#(x2, y2) = (1, 1) upper-right corner
#par(fig = c(0, 1, 0, 1),new = T) # whole area for test
par(fig = c(0.65,1, 0, 1),new = T) # Se plotar novamente, tem que ser desde o último comando par() acima
plot(tt,col=trans4$col,cex=0.4)
dev.off()

par(fig = c(0.5, 1, 0, 1),new = T)
plot(tt,col=trans4$col, # Cor dos tri?ngulos.
     cex=0.5, # Tamanho das letras
     lwd = 0.2, # Largura das linhas
     fixed.height=TRUE # Largura dos tri?ngulos (fixo=TRUE ou proporcional=FALSE)
)

#cria??o vetor cores #html notation no qgis na ordem que aparece no plot lateral
cores=c("#e6550d", "#a63603", "#005120", "#26fa05", "#fff701", "#08519c", "#006d2c", "#756bb1", "#969696", "#252525", "#636363", "#a50f15", "#54278f", "#de2d26", "#3182bd", "#31a354", "#ff1201", "#fbb4b9", "#feebe2", "#00fed8", "#fb6a4a", "#00c6a8") #substituir trans4$col por cores criado nesta linha

# Choose one of the par() parameters above before saving
pdf("UPGMA_rcluster_All_species.pdf")
par(fig = c(0.5, 1, 0, 1),new = T)
plot(tt,col=cores, # Cor dos tri?ngulos.
     cex=0.6, # Tamanho das letras
     lwd = 0.2, # Largura das linhas
     fixed.height=TRUE # Largura dos tri?ngulos (fixo=TRUE ou proporcional=FALSE)
)

plot(tt)
dev.off()

tiff("clados_All_species.tiff", units="in", width=5, height=5, res=300) ###salvar imagem em tiff
par(fig = c(0.5, 1, 0, 1),new = T)
plot(tt,col=cores, # Cor dos tri?ngulos.
     cex=0.6, # Tamanho das letras
     lwd = 0.2, # Largura das linhas
     fixed.height=TRUE # Largura dos tri?ngulos (fixo=TRUE ou proporcional=FALSE)
)
dev.off()


#dev.off() ###para fechar a imagem iniciada acima

# Salva o grid
writeOGR(classified_areas,"grid_rcluster_All_species.shp",driver="ESRI Shapefile",layer =classified_areas@data$class )

# Pega as coordenadas derivadas dos pol?gonos das spp e salva com um novo nome
pts_grids_cents_pol_2 = pontos_and_grids_centroides_pol
# Retira o elemento zzzzz
pts_grids_cents_pol_2 = pts_grids_cents_pol_2[!pts_grids_cents_pol_2$species %in% "zzzzz",]
# Transforma em um spatial points dataframe
coordinates(pts_grids_cents_pol_2) = pts_grids_cents_pol_2[,c("decimallongitude","decimallatitude")]
# Cruza os pontos com as ?reas classificadas por bioregion
spp_per_bioregion = raster::intersect(pts_grids_cents_pol_2,classified_areas)
# Prepara o dataframe
spp_per_bioregion_2 = unique(spp_per_bioregion@data[,c("species","class")])
spp_per_bioregion_2 = spp_per_bioregion_2[order(spp_per_bioregion_2$class),]
# Aqui est?o as classifica??es
spp_per_bioregion_2
# Salva em csv
write.csv(spp_per_bioregion_2, "spp_per_bioregion_All_species.csv")


# Getting additional colors when there is an excessive number of clusters (só para não perder o script)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
classified_areas
