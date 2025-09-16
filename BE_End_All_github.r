# Running BE analysis, as well as Species Richness calculations

setwd("~/Biogeografia_BE_polygons/All_End/")
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
#install.packages("Rtools", dependencies = TRUE) #BAIXA TODOS OS ASSOCIADOS

ocorrencias=read.table("END_All_July_21.txt", header = T)
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

# Produzindo a matriz de presença e ausência #######
# Na nova versão do SpeciesGeoCoderR, a função agora se chama SpGeoCod
### essa parte nova ? referente a species por grid antiga, qndo se criava arquivos em outra pasta

# calc ranges #210721 analise usando poligonos linhas 64 a 85
outp_all <- SpGeoCod(ocorrencias, grid_amazonia,areanames="layer")
spp_ranges = CalcRange(outp_all)
plot(spp_ranges)
# calc presence absence
crossing<- raster::intersect(spp_ranges,grid_amazonia)
crossing_dat<-crossing[,c("species","layer")]
crossing_dat<- na.exclude(crossing_dat)

# Getting coordinates from the raster centroids
grid_coord<- centroid(grid_amazonia)
grid_coord<- as.data.frame(grid_coord)
grid_coord$layer<-grid_amazonia@data$layer
colnames(grid_coord) = c("x","y","layer")

crossing_dat$long <- grid_coord$x[match(crossing_dat$layer,grid_coord$layer)]
crossing_dat$lat <- grid_coord$y[match(crossing_dat$layer,grid_coord$layer)]

spp_coords<- as.data.frame(crossing_dat[,c(1,3,4)])
spp_coords2 = spp_coords
colnames(spp_coords2) = colnames(trueCentroids)
pontos_and_grids_centroides <- rbind(spp_coords2,trueCentroids)

write.csv(pontos_and_grids_centroides, "coordenadas_from_shapes_End_All.csv")


SpGeoCod_object = SpGeoCod(pontos_and_grids_centroides, grid_amazonia, areanames = "layer")

# Veja que agora ele não salva a tabela no seu diretório, mas sim produz um objeto dentro do R que eu chamei de SpGeoCod_object. Agora eu preciso apenas da tabela de presença e ausência
pres_abs = SpGeoCod_object$spec_table

# Para visualizar as 5 primeiras linhas e 5 primeiras colunas
pres_abs[1:5,1:5]

pres_abs2 <-pres_abs[,order(as.numeric(colnames(pres_abs)),decreasing = FALSE ) ]

nrow(pres_abs2)
ncol(pres_abs2)
colnames(pres_abs2)
rownames(pres_abs2)

pres_abs3 = pres_abs2[!rownames(pres_abs2) %in% "zzzzz",]
pres_abs4 = pres_abs3[,!colnames(pres_abs3) %in% "not_classified"]
rownames(pres_abs4)
colnames(pres_abs4)

# Análise de elementos bióticos hierárquico ####
trueCentroids = centroid(grid_amazonia)
xx <- poly2nb(grid_amazonia,row.names = grid_amazonia@data$layer)
nrow(trueCentroids)
ff<- coord2dist(coordmatrix=trueCentroids[,2:1], file.format="decimal2",output.dist=TRUE, neighbors=FALSE)

pres_abs5 = pres_abs4
pres_abs5<-1*(pres_abs5>0)

#comando para rodar PAE no TNT, linhas 119 a 126
t(pres_abs5[1:10,1:10])
TNT=t(pres_abs5)
row.names(TNT)=paste0("grid", row.names(TNT))
write.table(TNT, "presabsTNT4.txt")
#dim(pres_abs5)
row.names(pres_abs5)
write.table(row.names(pres_abs5), "tabela_spp.txt")
colnames(pres_abs5)
class(pres_abs5)
which(pres_abs5[row.names(pres_abs5)=="Bachia_remota", ]==1)
#plotar grids especificos dentro do grid_amazonia
#par(mfrow=c(1,1))
#par(mar=c(5,4,4,2)+0.1) ##voltar ao parametros originais do plot-visual
#plot(mapa_amazonia)
#plot(ocorrencias_p_plot, pch=16, add = T) #pch simbolos
#plot(grid_amazonia, border="red")
#plot(grid_amazonia[grid_amazonia@data$layer=="896",], add=TRUE, col="black")
#plot(grid_amazonia[grid_amazonia@data$layer%in% c("896","895","1049","1011","974","973","1010","972"),], add=TRUE, col="black")

x$prab[1:100,1:10]

ncol(x$prab)

x <- prabinit(prabmatrix=pres_abs5,rows.are.species = FALSE, neighborhood=xx, geodist=(ff), distance="geco", gtf=0.1)

# 

library(prabclus)


#BE_significance1 = prabtest(x,times = 100,sf.sim=FALSE, sf.const=FALSE, 
                           #prange = c(0.3,0.5),
                           #pd = 0.9
#                           teststat = "nn")

#saveRDS(BE_significance1, "Desktop/BE_significance1.rds")
                           
#BE_significance2 = prabtest(x,times = 100,sf.sim=FALSE, sf.const=FALSE,
                           #prange = c(0.3,0.5),
                           #pd = 0.9
                           #teststat = "isovertice")

#BE_significance1 = readRDS("Desktop/BE_significance1.rds")

#summary(BE_significance)

#min(BE_significance1$results)
#max(BE_significance1$results)

x$regperspec[1:5,1:5]
x

summary(x)

library(prabclus)
#cutdist 0.5 ? o ponto de corte defaut dos grupamentos
d<-(hprabclust(x, cutdist=0.4, cutout=1, method="average", nnout=2, mdsplot=TRUE, mdsmethod="classical"))

#spp_per_be <- cbind.data.frame(d$rclustering, rownames(pres_abs3))
spp_per_be <- cbind.data.frame(d$rclustering, names(x$regperspec))
colnames(spp_per_be) = c("BE_number", "species")

spp_per_be <- spp_per_be[order(spp_per_be$BE_number),]
head(spp_per_be)
tail(spp_per_be)

write.csv(spp_per_be, "End_All_BE_results.csv")

# AQUI O BE não hierárquico (publicação original). Basta desmarcar os passos abaixo e seguir para os plots.

#d <- (prabclust(x, mdsmethod = "classical", mdsdim = 4, nnk =ceiling(x$n.species/40), nclus = 0:9, modelid = "all", permutations=10))
#spp_per_be <- cbind.data.frame(d$clustering, rownames(pres_abs3))
#colnames(spp_per_be) = c("BE_number", "species")
#spp_per_be <- spp_per_be[order(spp_per_be$BE_number),]
#spp_per_be
#write.csv(spp_per_be, "my_be2_results.csv")

# Plots e salvando shapefiles #####

# Aqui eu vou adicionar um loop para tentar facilitar a vida e produzir todas as figuras e escrever os SHPs de uma vez. Então, vai demorar um tempo para rodar e os resultados irão ser salvos no diretório atual.
# Se você quiser rodar um por um depois, pule a parte do loop (função for) e coloque o número do BE desejado no i=1,2,3... em seguida que está marcado com um #

# Seleciona o número de cada BE menos o 0 (noise)
BEs_for_loop = unique(spp_per_be$BE_number)[-1]
BEs_for_loop
# Inicia o plot em PDF
# Se quiser salvar plots com nomes diferentes, basta mudar o nome entre parentesis abaixo
pdf(paste0(gsub(":", "-", Sys.time()),"BE_cut04_End_All.pdf"))

# Inicia o loop
for(i in BEs_for_loop){
# i = 2 # o i ser? alterado no processo passo a passo, qndo se cria BE um por um, somente substituindo a numeracao no i. O zero nao significa nada
spp_BE_0 <- spp_per_be[spp_per_be$BE_number==i,]
spp_BE_0 <- as.character(spp_BE_0$species)
outp_BE_0 <- SpGeoCod(pontos_and_grids_centroides[pontos_and_grids_centroides$species %in% spp_BE_0,], grid_amazonia,areanames="layer")

# criando o grid de riqueza para os BEs. Mesmo que acima porém servirá para produzir nossos polígonos
map_BE_0<- RichnessGrid(outp_BE_0, ras = grid_amazonia_raster, type = "spnum")
proj4string(map_BE_0) = proj4string(grid_amazonia)
#plot(map_BE_0) #qndo fazendo loop, nao pode estar ativo a funcao, assim como o i acima

# Produzindo grid para plotar
# Alguns ajustes extras para o raster ficar na extensão correta

map_BE_0 = mask(map_BE_0,grid_amazonia_raster)

# Aqui para definir se os BEs tem uma ou duas spp por grid
cut_val = ifelse(maxValue(map_BE_0)>1,2,1)

map_BE_0[map_BE_0 <cut_val] = NA

# Aqui serve para plotar os mapas lado a lado
#The mai argument specifies the margin by sides, i.e., c(bottom, left, top, right).
par(mfrow=c(2,2), tcl=-0.5, family="serif", mai=c(0.3,0.2,0.5,0.3))

# Se quiser voltar ao normal, basta usar:
#par(mfrow=c(1,1))
#par(mar=c(5,4,4,2)+0.1)

# Plotando apenas os polígonos com mais de duas spp
#plot(map_BE_0, col="white", legend=FALSE,axes=FALSE, box=FALSE)
#plot(outp_BE_0$polygons, border="lightgray")
#plot(map_BE_0, col = "red",add=T)
#plot(outp_BE_0$polygons,border="lightgray",add=T)
#maps::map('world',add=TRUE)
#title(paste0("Biotic Element ", i, "\n(two or more species)"))
# Adiciona legenda para os pontos da spp

# Default plot settings
#par(mfrow=c(1,1))
#par(mar=c(5,4,4,2)+0.1)
# Selecionando species points e colorindo
spp_plot_map_BE = ocorrencias_p_plot[ocorrencias_p_plot@data$species %in% spp_BE_0,]

#if (length(spp_BE_0)<3){
#  cores_points = brewer.pal(length(spp_BE_0), "Set3")[-1]
#} else {
#  cores_points = brewer.pal(length(spp_BE_0), "Set3")
#}

spp_plot_map_BE$col =  plyr::mapvalues(spp_plot_map_BE$species, from=spp_BE_0, to=magma(length(spp_BE_0))) 

plot(map_BE_0, col="white", legend=FALSE,axes=FALSE, box=FALSE)
plot(spp_plot_map_BE, pch=16, col = spp_plot_map_BE$col, add=TRUE)
plot(outp_BE_0$polygons,border="lightgray",add=T)
plot(rasterToPolygons(map_BE_0), border="black", lwd=2, add=T)
maps::map('world',add=TRUE)
title(paste0("Biotic Element ", i, "\n(",cut_val," or more species)"))

# Utilizando os ranges das spp apenas para melhor delimitar os BEs
#map_BE_0_ranges = CalcRange(ocorrencias[ocorrencias$species %in% spp_BE_0,], method = "pseudospherical", terrestrial = TRUE, rare = "buffer", buffer.width = 100000)

#map_BE_0_range = RangeRichness(map_BE_0_ranges,ras = grid_amazonia_raster, terrestrial = TRUE,buffer=100000)
#map_BE_0_range[map_BE_0_range <1] = NA
#map_BE_0_range = resample(map_BE_0_range,grid_amazonia_raster,method="ngb")
#map_BE_0_range = mask(map_BE_0_range,grid_amazonia_raster)

plot(map_BE_0, col="white", legend=FALSE,axes=FALSE, box=FALSE)
#plot(outp_BE_0$polygons, border="lightgray")
#plot(map_BE_0_range, col = rev(magma(maxValue(map_BE_0)+1)[-(maxValue(map_BE_0)+1)]),add=T)
#plot(map_BE_0_range, col = rev(viridis(maxValue(map_BE_0))),add=T)
#plot(outp_BE_0$polygons,border="lightgray",add=T)
#maps::map('world',add=TRUE)
#title(paste0("Biotic Element ", i, "\n(species richness- range maps)"))
legend("topleft", legend = gsub("_"," ", spp_BE_0), pch=16,col= unique(spp_plot_map_BE$col), cex = .6, bty = "n",text.font = 3,ncol = 2) #tiramos as legendas de especies, caso decida inclui-la, desmarcar o #

# Plotando apenas os grids com 50% das spp
map_BE_0_2 = map_BE_0
map_BE_0_2[map_BE_0_2 < maxValue(map_BE_0_2)/2] = NA
plot(map_BE_0, col="white", legend=FALSE,axes=FALSE, box=FALSE)
#plot(outp_BE_0$polygons, border="lightgray")
plot(map_BE_0_2, col = "red",add=T,legend=FALSE)
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Biotic Element ", i, "\n(half of species)"))

# Plotando com riqueza
plot(map_BE_0, col="white", legend=FALSE,axes=FALSE, box=FALSE)
#plot(outp_BE_0$polygons, border="lightgray")
#plot(map_BE_0, col = rev(magma(maxValue(map_BE_0)+1)[-(maxValue(map_BE_0)+1)]),add=T)
plot(map_BE_0, col = rev(viridis(maxValue(map_BE_0))),add=T)

plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Biotic Element ", i, "\n(species richness)"))




#abaixo para salvar apenas o BE discriminado
map_BE_0_pol = rasterToPolygons(map_BE_0)
writeSpatialShape(map_BE_0_pol, paste0("BE_End_All_", i, ".shp"))


map_BE_0_2_pol = rasterToPolygons(map_BE_0_2)
writeSpatialShape(map_BE_0_2_pol, paste0("BE_End_All_50p_", i, ".shp"))


map_BE_0_rich = rasterToPolygons(map_BE_0)
writeSpatialShape(map_BE_0_rich, paste0("BE_End_All_rich_", i, ".shp"))

#writeRaster(map_BE_0_2, paste0("BE_End_All_50p_", i, ".tif"), overwrite=TRUE)
#writeRaster(map_BE_0, paste0("BE_End_All_rich_", i, ".tif"), overwrite=TRUE)


}


# Fim do loop, vai demorar um tempo para finalizar. Não esquece de rodar a próxima linha antes de procurar pelo PDF no seu diretório

# Salva o pdf com múltiplas páginas (verifique no diretório) #APOS O LOOP PARA FECHAR E SALVAR O PDF
dev.off()

# Se quiser mudar as cores da legenda de riqueza, substitua o magma por
# magma(), plasma(), inferno(), cividis(), mako(), rocket()
#####

#
##############################################
# New! Mapeando a riqueza de espécies e tendencia de coleta ####
#############################################

# Default plot settings
par(mfrow=c(1,2))

outp_all <- SpGeoCod(ocorrencias, grid_amazonia,areanames="layer")

# Riqueza com o raw data
pdf("Riqueza_End_All.pdf")#somente para salvar direto, sem aparecer no R
richness_all<- RichnessGrid(outp_all,ras = grid_amazonia_raster, type = "spnum")
# Cortando para incluir apenas a área dentro do grid
#richness_all = resample(richness_all,grid_amazonia_raster,method="ngb")
richness_all = mask(richness_all,grid_amazonia_raster)
plot(richness_all, col = rev(viridis(maxValue(richness_all))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Species richness", "\n(Raw point localities)"))

# Sampling effort (number of records per grid)
sampling_all<- RichnessGrid(outp_all, ras = grid_amazonia_raster, type = "abu")
# Cortando para incluir apenas a área dentro do grid
#sampling_all = resample(sampling_all,grid_amazonia_raster,method="ngb")
sampling_all = mask(sampling_all,grid_amazonia_raster)
plot(sampling_all, col = rev(viridis(maxValue(sampling_all))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Sampling effort", "\n(Raw point localities)"))
dev.off() #s? para salvar em pdf

# Correlação species richness e sampling effort (não funcionou bem aqui, pular)
#plot(corLocal(richness_all,richness_all,ngb=3,test=TRUE),col = rev(viridis(6)))

# Default plot settings
par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)

# Riqueza polígonos
pdf("Riqueza_Poly_End_All.pdf")#somente para salvar direto, sem aparecer no R
All_ranges = CalcRange(outp_all, method = "pseudospherical", terrestrial = FALSE, rare = "buffer", buffer.width = 100000)
all_range = RangeRichness(All_ranges, ras = grid_amazonia_raster, terrestrial = FALSE)
#all_range = resample(all_range,grid_amazonia_raster,method="ngb")
all_range = mask(all_range,grid_amazonia_raster)
plot(all_range, col = rev(viridis(maxValue(all_range))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
maps::map('world',add=TRUE)
title(paste0("Species richness", "\n(range polygons)"))
#writeRaster(all_range, paste0("Riqueza_Pol_End_All.tif"), overwrite=TRUE)

all_range_pol = rasterToPolygons(all_range) #para criar shape, primeiro precisa transformar de raster para shape, depois abaixo criar o shape
writeOGR(all_range_pol, ".", "Richness_End_All_Poly.shp", driver="ESRI Shapefile")

par(mfrow=c(1,2))
plot(sampling_all, col = rev(viridis(maxValue(sampling_all))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Sampling effort", "\n(Raw point localities)"))
plot(all_range, col = rev(viridis(maxValue(all_range))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
maps::map('world',add=TRUE)
title(paste0("Species richness", "\n(range polygons)"))
dev.off() #s? para salvar em pdf
######

#
# Calculando Weighted Endemism raw data ####

install.packages("remotes")
remotes::install_version("SDMTools", "1.1-221")
library(biomapME)
devtools::install_github("GregGuerin/biomapME", dependencies=TRUE)


WE_reptiles_raw = weighted.endemism(species_records = ocorrencias, records = "single", species = "species", longitude = "decimallongitude", latitude = "decimallatitude",frame.raster= grid_amazonia_raster, type = "weighted",plot.raster = FALSE,  weight.type = "cell", geo.type = "cell", outlier_pct = 100, verbose = TRUE)

WE_reptiles_raw_plot = WE_reptiles_raw$WE_raster
WE_reptiles_raw_plot = resample(WE_reptiles_raw_plot,grid_amazonia_raster,method="ngb")
WE_reptiles_raw_plot = mask(WE_reptiles_raw_plot,grid_amazonia_raster)
# Default plot settings
par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)
plot(WE_reptiles_raw_plot, col = rev(viridis(100)))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Weighted endemism", "\n(Raw point localities)"))

# Species richness and weighted endemism
pdf("Weighted_End_End_All.pdf")#somente para salvar direto, sem aparecer no R
par(mfrow=c(1,2))
plot(richness_all, col = rev(viridis(maxValue(richness_all))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Species richness", "\n(Raw point localities)"))
plot(WE_reptiles_raw_plot, col = rev(viridis(100)))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Weighted endemism", "\n(Raw point localities)"))
dev.off() #s? para salvar em pdf
#######

# WE range polygons #####
# calc ranges

spp_ranges = CalcRange(outp_all)
# calc presence absence
crossing<- raster::intersect(spp_ranges,grid_amazonia)
crossing_dat<-crossing[,c("species","layer")]
crossing_dat<- na.exclude(crossing_dat)
#rm(crossing)
# Getting coordinates from the raster centroids
grid_coord<- centroid(grid_amazonia)
grid_coord<- as.data.frame(grid_coord)
grid_coord$layer<-grid_amazonia@data$layer
colnames(grid_coord) = c("x","y","layer")
#head(grid_coord)
crossing_dat$long <- grid_coord$x[match(crossing_dat$layer,grid_coord$layer)]
crossing_dat$lat <- grid_coord$y[match(crossing_dat$layer,grid_coord$layer)]
#tail(crossing_dat)
spp_coords<- as.data.frame(crossing_dat[,c(1,3,4)])

# Agora calcula o WE a partir dos ranges para comparar com os pols
pdf("Weighted_End_Poly_End_All.pdf")#somente para salvar direto, sem aparecer no R
WE_reptiles_span = weighted.endemism(species_records = spp_coords, records = "single", species = "species", longitude = "long", latitude = "lat",frame.raster= grid_amazonia_raster, type = "weighted",plot.raster = FALSE,  weight.type = "cell", geo.type = "cell", outlier_pct = 100, verbose = TRUE)

WE_reptiles_span_plot = WE_reptiles_span$WE_raster
WE_reptiles_span_plot = resample(WE_reptiles_span_plot,grid_amazonia_raster,method="ngb")
WE_reptiles_span_plot = mask(WE_reptiles_span_plot,grid_amazonia_raster)
# Plot species ranges
plot(all_range, col = rev(viridis(maxValue(all_range))))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Species richness", "\n(range polygons)"))
plot(WE_reptiles_span_plot, col = rev(viridis(100)))
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Weighted endemism", "\n(range polygons)"))
dev.off() #s? para salvar em pdf

# Default plot settings
par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)

# Testes para Weighted endemism (maior ou menor que o esperado em comparação com a riqueza de espécies?)

# Para o raw data
pdf("Weighted_End_RawRanges_test_End_All.pdf")#somente para salvar direto, sem aparecer no R
WE_raw_test = endemism.null.test(WE_reptiles_raw)

# Espera terminar acima antes de rodar esses passos

WE_raw_test_plot = resample(WE_raw_test$out.above.below.raster,grid_amazonia_raster,method="ngb")
WE_reptiles_span_plot = mask(WE_raw_test_plot,grid_amazonia_raster)
plot(WE_reptiles_span_plot,col = rev(viridis(10)), legend=FALSE)
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Significant WE areas", "\n(raw data)"))
legend(x=c(-57,-57), y=c(20,-10), legend = c("higher","lower"), pch=15, col = c("#440154FF","#FDE725FF"), cex = 1.2, bty = "n")

# Para range data
WE_span_test = endemism.null.test(WE_reptiles_span)

# Espera terminar para plotar abaixo
WE_span_test_plot = resample(WE_span_test$out.above.below.raster,grid_amazonia_raster,method="ngb")
WE_span_test_plot = mask(WE_span_test_plot,grid_amazonia_raster)
plot(WE_span_test_plot,col = rev(viridis(10)), legend=FALSE)
plot(outp_BE_0$polygons,border="lightgray",add=T)
maps::map('world',add=TRUE)
title(paste0("Significant WE areas", "\n(range polygons)"))
legend(x=c(-57,-57), y=c(20,-10), legend = c("higher","lower"), pch=15, col = c("#440154FF","#FDE725FF"), cex = 1.2, bty = "n")
dev.off()
######

# Vamos produzir uma matrix de disimilaridade para regionalização seguindo um script to Mario Moura para a Mata Altantica ####
# Precisamos de uma matrix de dissimilaridade apenas entre os grids com relativamente um bom número de spp. Vou usar a Sorensen pair-wise dissimilarity.

# Precisamos apenas da matrix de presença e ausência
pres_abs5 = pres_abs4
# Transformar em presença e ausência
pres_abs5<-1*(pres_abs5>0)

# Seleciona sites com três espécies ou mais.
# Aqui, para a análise seguinte, eu recomendo usar um número de espécies que seria considerado relativamente real para uma comunidade de cada grupo de répteis. Eu coloquei 3 amphisbaenas (>2 abaixo).
pres_abs6 = pres_abs5[,which(colSums(pres_abs5)>5)]


all_spp.betapair<-beta.pair(t(pres_abs6), index.family="sor")#objeto criado da matrix de dissimilaridade usando o presenca ausencia 6
nrow(as.matrix(all_spp.betapair$beta.sor))
ncol(as.matrix(all_spp.betapair$beta.sor))

all_spp.betapair.beta.turnover <- as.data.frame(as.matrix(all_spp.betapair$beta.sor))
all_spp.betapair.beta.turnover[1:10,1:10]

#daqui pra baixo se referem ao cluster, acima a planilha pro proximo script
# Aqui eu vou usar a dissimilaridade total ($beta.sor na linha abaixo), mas você pode dividir em turnover ($beta.sim) ou nesdtedness (beta.sne) e produzir três mapas no usando o script do Moura
# 

# Adiciona long, lat e layer name
all_spp.betapair.beta.turnover$long <- grid_coord$x[match(rownames(all_spp.betapair.beta.turnover),grid_coord$layer)]
all_spp.betapair.beta.turnover$lat <- grid_coord$y[match(rownames(all_spp.betapair.beta.turnover),grid_coord$layer)]
all_spp.betapair.beta.turnover$layer <- grid_coord$layer[match(rownames(all_spp.betapair.beta.turnover),grid_coord$layer)]

# Agora mude para o script Moura_et_al_S1

# Site para gerar dados para o programa NDM/VNDM http://gex.mfuhlendorf.com/

##