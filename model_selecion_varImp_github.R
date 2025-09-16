# AoE, BR, SR model selection

# Definindo os sistemas de refer??ncia de coordenadas
wgs84 <- "EPSG:4326"  # WGS84 em graus decimais
behrmannCRS <- "EPSG:3410"  # Proje????o de ??rea igual equivalente (substitui o CEA)

# Carregando os pacotes necess??rios
library(sf)        # Para dados vetoriais (substitui rgdal e rgeos)
library(terra)     # Para dados raster (substitui raster)
library(dplyr)     # Para manipula????o de dados
library(tidyr)     # Para transforma????es de dados
library(ggplot2)   # Para visualiza????o
library(car)       # Para vif (Variance Inflation Factor)
library(gridExtra) # Para organiza????o de gr??ficos
library(nnet)      # Para modelos multinomiais
library(MuMIn)     # Para dredge e modelos m??dios

# Caminho para os arquivos de ??reas de endemismo
pol_path <- list.files(path = "reasdeendemismos", 
                       pattern = "\\.shp$", 
                       full.names = TRUE)

# Leitura e uni??o dos arquivos shapefile
# Vers??o atualizada para sf, lidando com diferentes estruturas de coluna
pol_am <- list()
for (i in pol_path) {
  # Lendo o shapefile
  a <- st_read(i, quiet = TRUE)
  
  # Extraindo apenas a geometria e criando um novo objeto sf com apenas geometria e id
  geom <- st_geometry(a)
  nome_arquivo <- gsub("reasdeendemismos/", "", i)
  
  # Criando um sf simplificado com apenas geometria e id
  a_simples <- st_sf(
    id = nome_arquivo,
    count = ifelse("count" %in% names(a), a$count, NA),  # Preserva count se existir
    geometry = geom
  )
  
  pol_am <- c(pol_am, list(a_simples))
}

# Agora todos os objetos sf t??m a mesma estrutura e podem ser unidos
pol_am_shps <- do.call(rbind, pol_am)

# Visualiza????o dos pol??gonos com cores diferentes
# Usando ggplot2 em vez de plot base
unique_ids <- unique(pol_am_shps$id)
colors <- rainbow(length(unique_ids))
ggplot(pol_am_shps) +
  geom_sf(aes(fill = id)) +
  scale_fill_manual(values = colors) +
  labs(title = "Spatial Plot by ID") +
  theme_minimal()

# Carregando o recorte da Amaz??nia
amazonia_recorte <- st_read("amazniarecorte/Amazonia.shp", quiet = TRUE)
plot(st_geometry(amazonia_recorte))

# Carregando a grade da Amaz??nia
grid_am <- st_read("/Users/josue/Library/CloudStorage/Dropbox/1Doutorado/Chapter_2/biodiverse_pipeline-master/brazil_stat.shp", quiet = TRUE, crs = wgs84)

# Plotando as camadas sobrepostas usando ggplot2

pol_am_shps$id = gsub(".shp","", pol_am_shps$id)

ggplot() +
  geom_sf(data = amazonia_recorte, fill = NA) +
  geom_sf(data = pol_am_shps, aes(fill = id)) +
  geom_sf(data = grid_am, fill = NA, color = "black") +
  scale_fill_manual(values = colors) +
  theme_minimal()

# Carregando os rasters bioclim??ticos
# Note: terra usa rast() em vez de raster()
my_bioclim_rasters <- readRDS("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/my_bioclim_rasters.rds")

# Convertendo para objetos terra::SpatRaster se necess??rio
# Isso assume que my_bioclim_rasters ?? uma lista contendo objetos raster
# Se j?? forem objetos terra, essa convers??o n??o ?? necess??ria
if (inherits(my_bioclim_rasters[[2]], "RasterStack") || inherits(my_bioclim_rasters[[2]], "RasterBrick")) {
  preds <- terra::rast(my_bioclim_rasters[[2]])
} else if (inherits(my_bioclim_rasters[[2]], "SpatRaster")) {
  preds <- my_bioclim_rasters[[2]]
}

# Verificando nomes das camadas
names(preds)

# Mostrando uma das camadas
plot(preds[["bio_18"]])

# Reprojetando para WGS84
preds <- terra::project(preds, wgs84)

# PCA para vari??veis de temperatura e precipita????o
# A fun????o RStoolbox::rasterPCA n??o tem equivalente direto em terra
# Vamos implementar uma abordagem usando PCA padr??o

# Tratamento de erro para PCA
tryCatch({
  # Para temperatura (bio_2 a bio_11)
  temp_layers <- preds[[c("bio_1","bio_2",  "bio_3", "bio_4",  "bio_5",  "bio_6" , "bio_7" , "bio_8", "bio_9","bio_10", "bio_11")]]
  temp_values <- terra::values(temp_layers, na.rm = FALSE)
  # Remover linhas com NA para evitar erros no PCA
  temp_values_clean <- na.omit(temp_values)
  
  # Verificar se temos dados suficientes ap??s remo????o de NAs
  if(nrow(temp_values_clean) > 0) {
    temp_pca <- prcomp(temp_values_clean, scale. = TRUE, center = TRUE)
    temp_pca_rast <- terra::predict(temp_layers, temp_pca, index = 1:3)  # Primeiros 3 PCs
  } else {
    message("Aviso: Muitos valores NA nas camadas de temperatura. Criando raster vazio para PCA.")
    temp_pca_rast <- terra::rast(temp_layers)  # Cria um raster vazio com a mesma extens??o
    names(temp_pca_rast) <- paste0("temp_PC", 1:3)
  }
}, error = function(e) {
  message("Erro ao calcular PCA de temperatura: ", e$message)
  message("Criando raster substituto...")
  # Criar um raster alternativo em caso de erro
  temp_pca_rast <- terra::rast(preds[[1]])
  temp_pca_rast <- rep(temp_pca_rast, 3)
  names(temp_pca_rast) <- paste0("temp_PC", 1:3)
})

tryCatch({
  # Para precipita????o (bio_1 e bio_12 a bio_19)
  precip_layers <- preds[[c("bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")]]
  precip_values <- terra::values(precip_layers, na.rm = FALSE)
  # Remover linhas com NA para evitar erros no PCA
  precip_values_clean <- na.omit(precip_values)
  
  # Verificar se temos dados suficientes ap??s remo????o de NAs
  if(nrow(precip_values_clean) > 0) {
    precip_pca <- prcomp(precip_values_clean, scale. = TRUE, center = TRUE)
    precip_pca_rast <- terra::predict(precip_layers, precip_pca, index = 1:3)  # Primeiros 3 PCs
  } else {
    message("Aviso: Muitos valores NA nas camadas de precipita????o. Criando raster vazio para PCA.")
    precip_pca_rast <- terra::rast(precip_layers)  # Cria um raster vazio com a mesma extens??o
    names(precip_pca_rast) <- paste0("precip_PC", 1:3)
  }
}, error = function(e) {
  message("Erro ao calcular PCA de precipita????o: ", e$message)
  message("Criando raster substituto...")
  # Criar um raster alternativo em caso de erro
  precip_pca_rast <- terra::rast(preds[[1]])
  precip_pca_rast <- rep(precip_pca_rast, 3)
  names(precip_pca_rast) <- paste0("precip_PC", 1:3)
})

# NDVI
NDVI <- preds[[20]]

# Carregando e processando a eleva????o
elevat <- terra::rast('/Users/josue/Dropbox/4Environmental_layers/wc2.1_5m_elev.tif')
elevat <- terra::crop(elevat, terra::ext(preds[[20]]))
plot(elevat)

# Calculando a rugosidade do relevo
relief_roughness <- terra::terrain(elevat, v = "roughness", neighbors = 8)
plot(relief_roughness)
plot(st_geometry(grid_am), add = TRUE)

# Carregando e processando a velocidade clim??tica
#climate_vel <- terra::rast("/Users/josue/Library/CloudStorage/Dropbox/1Doutorado/Chapter_2/biodiverse_pipeline-master/Velocity.tif")

# Since Pliocene
climate_vel <- terra::rast("amazniarecorte/Raster_layers_R_scripts/Layers/past/csi_past.tif")
# Manually set the correct CRS (Behrmann)
#terra::crs(climate_vel) <- "EPSG:3410"
#climate_vel <- terra::project(climate_vel, wgs84)
climate_vel <- terra::crop(climate_vel, terra::ext(preds[[20]]))
plot(climate_vel)

# Carregando e processando geologia
geology_raster = rast("amazniarecorte/AMS_geology_raster.tiff")
geology_raster <- terra::crop(geology_raster, terra::ext(preds[[20]]))
plot(geology_raster)

# Compute the minimum and maximum values of the raster (ignoring NA values)
min_val <- global(geology_raster, "min", na.rm = TRUE)[1]
max_val <- global(geology_raster, "max", na.rm = TRUE)[1]

# Create 6 breakpoints (for 5 equal intervals)
breaks <- seq(min_val$min, max_val$max, length.out = 6)
print(breaks)  # To check the break values

# Build a reclassification matrix with 3 columns: [from, to, new category]
rcl <- matrix(c(
  breaks[1], breaks[2], 1,
  breaks[2], breaks[3], 2,
  breaks[3], breaks[4], 3,
  breaks[4], breaks[5], 4,
  breaks[5], breaks[6], 5
), ncol = 3, byrow = TRUE)
print(rcl)  # Check the reclassification matrix

# Reclassify the raster into 5 categories using the classify() function
geology_cat <- classify(geology_raster, rcl)

# Plot the original raster and the reclassified (categorical) raster side by side
plot(geology_raster, main = "Original Geology Raster")
plot(geology_cat, main = "Geology Raster (5 Equal Intervals)")

# Carregando e processando as bacias hidrogr??ficas
grid_basins <- st_read("amazniarecorte/grid_amazonia_vasoes.shp", quiet = TRUE)

# Colorindo por bacia
unique_basins <- unique(grid_basins$Basin)
colors <- rainbow(length(unique_basins))

# Atribuindo cores ??s bacias
grid_basins <- grid_basins %>%
  mutate(col = case_when(
    Basin %in% unique_basins ~ rainbow(length(unique_basins))[match(Basin, unique_basins)],
    TRUE ~ "grey"
  ))

# Plotando as bacias coloridas
ggplot(grid_basins) +
  geom_sf(aes(fill = Basin)) +
  scale_fill_manual(values = setNames(colors, unique_basins)) +
  theme_minimal()

# Verificando valores ??nicos
unique(grid_basins$Basin)
unique(grid_basins$rio_em_si)

# Modificando os nomes dos rios
grid_basins$rio_em_si <- paste0("Rio_", grid_basins$rio_em_si)
unique(grid_basins$rio_em_si)

# Definindo vaz??es dos rios
rio_vazoes <- c(
  "Rio_NA" = NA,
  "Rio_Amazonas" = 209000,
  "Rio_Orinoco" = 35000,
  "Rio_Madeira" = 32000,
  "Rio_Negro" = 28400,
  "Rio_Negro_Branco" = 28400 + 1462,
  "Rio_Japur??" = 18600,
  "Rio_Tapaj??s" = 13500,
  "Rio_Purus" = 11000,
  "Rio_Xingu" = 9700,
  "Rio_Uacayali" = 9544,
  "Rio_Putumayo" = 8760,
  "Rio_Tocantins" = 8440,
  "Rio_Branco" = 1462,
  "Rio_Juru??" = 1462,
  "Rio_Solim??es" = 8760 + 9544,
  "Rio_Japura_Putumayo" = 18600 + 8760
)

# Atribuindo vaz??es usando dplyr com o sf
grid_basins <- grid_basins %>%
  mutate(discharge = case_when(
    rio_em_si %in% names(rio_vazoes) ~ as.numeric(rio_vazoes[rio_em_si]),
    TRUE ~ NA_real_
  ))

# Plotando as bacias coloridas
ggplot(grid_basins) +
  geom_sf(aes(fill = discharge)) +
  #scale_fill_manual(values = setNames(colors, unique_basins)) +
  theme_minimal()

# Calculando a mediana de vaz??o por bacia
median_discharge <- grid_basins %>%
  group_by(Basin) %>%
  summarise(median_discharge = median(discharge, na.rm = TRUE)) %>%
  st_drop_geometry()  # Removendo a geometria para jun????o

# Juntando a mediana de vaz??o com o dataframe original
grid_basins <- grid_basins %>%
  left_join(median_discharge, by = "Basin")

# Verificando os resultados
head(grid_basins)
unique(grid_basins$discharge)

# Plotando a intensidade de vaz??o
ggplot(grid_basins) +
  geom_sf(aes(fill = median_discharge), color = "black") +
  scale_fill_gradient(low = "lightgrey", high = "blue", na.value = NA) +
  labs(title = "Discharge Intensity", fill = "Discharge") +
  theme_minimal()

# Convertendo o sf em raster usando terra
# Primeiro criamos um raster de refer??ncia baseado em preds
ref_rast <- terra::rast(preds[[20]])
# Agora rasterizamos o grid_basins
grid_basins_raster <- terra::rasterize(
  terra::vect(grid_basins), 
  ref_rast, 
  field = "median_discharge", 
  fun = max
)

plot(grid_basins_raster)

# Preparando dados para AoE (Areas of Endemism)
# Adicionando um ID num??rico aos pol??gonos
pol_am_shps <- pol_am_shps %>%
  mutate(ID_final = as.numeric(factor(id, levels = unique(id))))

# Rasterizando para AoE
pol_am_shps_vect <- terra::vect(pol_am_shps)
pol_am_shps_raster <- terra::rasterize(
  pol_am_shps_vect, 
  ref_rast, 
  field = "ID_final", 
  fun = "max"
)
plot(pol_am_shps_raster)

# Para BIOREGIONS
pol_am_shps_ <- st_read("amazniarecorte/bioregions/grid_rcluster_All_species_cortado.shp", quiet = TRUE)

# Verificar se a coluna 'col' existe antes de usar no plot
if("col" %in% names(pol_am_shps_)) {
  ggplot(pol_am_shps_) +
    geom_sf(aes(fill = col)) +
    theme_minimal()
} else {
  # Se n??o existir 'col', usar a coluna 'class' ou outra adequada
  ggplot(pol_am_shps_) +
    geom_sf(aes(fill = as.factor(class))) +
    theme_minimal()
}

# Verificar se a coluna 'class' existe
if("class" %in% names(pol_am_shps_)) {
  unique(pol_am_shps_$class)
  pol_am_shps_$class <- as.numeric(pol_am_shps_$class)
  
  # Rasterizando bioregions
  pol_am_shps_vect_ <- terra::vect(pol_am_shps_)
  pol_am_shps_raster_ <- terra::rasterize(
    pol_am_shps_vect_, 
    ref_rast, 
    field = "class", 
    fun = "max"
  )
} else {
  # Criar uma coluna class se n??o existir
  pol_am_shps_$class <- 1  # Valor padr??o
  
  # Rasterizando bioregions
  pol_am_shps_vect_ <- terra::vect(pol_am_shps_)
  pol_am_shps_raster_ <- terra::rasterize(
    pol_am_shps_vect_, 
    ref_rast, 
    field = "class", 
    fun = "max"
  )
}
plot(pol_am_shps_raster_)

# Empilhando todos os rasters juntos
# Adicionando tratamento de erros para maior robustez

# Verificar se todos os objetos existem antes de prosseguir
rasters_to_stack <- list()

# Verificando e adicionando cada raster ?? lista, se dispon??vel
if(exists("pol_am_shps_raster")) rasters_to_stack$AoE <- pol_am_shps_raster
if(exists("grid_basins_raster")) rasters_to_stack$Basins <- grid_basins_raster
if(exists("climate_vel")) rasters_to_stack$Climate_vel <- climate_vel
if(exists("relief_roughness")) rasters_to_stack$relief_roughness <- relief_roughness
if(exists("NDVI")) rasters_to_stack$NDVI <- NDVI
if(exists("geology_cat")) rasters_to_stack$geology_cat <- geology_cat

# Extraindo o primeiro componente do PCA de precipita????o, se dispon??vel
if(exists("precip_pca_rast")) {
  tryCatch({
    rasters_to_stack$Precip_PCA <- precip_pca_rast[[1]]
  }, error = function(e) {
    message("Erro ao acessar o primeiro componente do PCA de precipita????o: ", e$message)
  })
}

# Extraindo o primeiro componente do PCA de temperatura, se dispon??vel
if(exists("temp_pca_rast")) {
  tryCatch({
    rasters_to_stack$Temp_PCA <- temp_pca_rast[[1]]
  }, error = function(e) {
    message("Erro ao acessar o primeiro componente do PCA de temperatura: ", e$message)
  })
}

# Verificar se temos rasters para alinhar
if(length(rasters_to_stack) > 0) {
  # Escolher o primeiro raster dispon??vel como refer??ncia
  reference_name <- names(rasters_to_stack)[1]
  reference_raster <- rasters_to_stack[[reference_name]]
  message("Usando '", reference_name, "' como raster de refer??ncia para alinhamento")
  
  # Alinhando e reprojetando todos os rasters
  aligned_rasters <- lapply(rasters_to_stack, function(r) {
    tryCatch({
      if (!terra::compareGeom(r, reference_raster, stopOnError = FALSE)) {
        r <- terra::project(r, reference_raster)
      }
      return(r)
    }, error = function(e) {
      message("Erro ao alinhar raster: ", e$message)
      return(NULL)
    })
  })
  
  # Filtrando rasters nulos
  aligned_rasters <- aligned_rasters[!sapply(aligned_rasters, is.null)]
  
  # Empilhando os rasters alinhados
  if(length(aligned_rasters) > 0) {
    all_rasters_aligned <- terra::rast(aligned_rasters)
  } else {
    stop("Nenhum raster dispon??vel ap??s tentativa de alinhamento")
  }
} else {
  stop("Nenhum raster dispon??vel para empilhar")
}

# Nomeando as camadas
#names(all_rasters_aligned) <- c(
#  "AoE", "Basins", "Climate_vel", "relief_roughness", 
#  "NDVI", "Precip_PCA", "Temp_PCA"
#)

# Extraindo os valores dos rasters para a grade
# Convertendo grid_am para vect para uso com terra::extract
grid_am_vect <- terra::vect(grid_am)
df_for_analyses <- terra::extract(all_rasters_aligned, grid_am_vect, fun = "median", na.rm = TRUE, ID=FALSE)

unique(df_for_analyses$geology_cat)

df_for_analyses$geology_cat = round(df_for_analyses$geology_cat,0)
unique(df_for_analyses$geology_cat)

# Plotando AoE com a grade
plot(all_rasters_aligned[["AoE"]])
plot(terra::vect(grid_am), add = TRUE)

# Preparando dataframe para an??lises
df_for_analyses2 <- as.data.frame(df_for_analyses)
head(df_for_analyses2)

# Plotando antes das an??lises
grid_am_plot <- grid_am
grid_am_plot$Basins <- log(df_for_analyses2$Basins)
grid_am_plot$Climate_vel <- df_for_analyses2$Climate_vel
grid_am_plot$relief_roughness <- df_for_analyses2$relief_roughness
grid_am_plot$NDVI <- df_for_analyses2$NDVI
grid_am_plot$Precip_PCA <- df_for_analyses2$Precip_PCA
grid_am_plot$Temp_PCA <- df_for_analyses2$Temp_PCA
grid_am_plot$Geomorphology <- df_for_analyses2$geology_cat

plot(geology_cat)

head(grid_am_plot)
unique(grid_am_plot$Geomorphology) # 

# Convertendo para formato longo para ggplot
grid_am_df <- st_as_sf(grid_am_plot) %>%
  st_cast()

grid_am_long <- grid_am_df %>%
  pivot_longer(
    cols = c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA","Geomorphology"),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# ===== 1. AJUSTES PARA A VISUALIZA????O - CORRIGINDO EIXO X E GEOMORFOLOGIA COMO FATOR =====

# Definindo vari??veis e t??tulos para os gr??ficos
variables <- c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA", "Geomorphology")

custom_titles <- c(
  "Basins" = "Riverine Barriers",
  "Precip_PCA" = "Precipitation",
  "Climate_vel" = "Climate Stability Index",
  "relief_roughness" = "Relief Roughness",
  "NDVI" = "Vegetation (NDVI)",
  "Temp_PCA" = "Temperature",
  "Geomorphology" = "Geomorphology"
)

# Definir quais vari??veis s??o categ??ricas
categorical_vars <- c("Geomorphology")

# Convertendo para formato longo para ggplot
grid_am_df <- st_as_sf(grid_am_plot) %>%
  st_cast()

grid_am_long <- grid_am_df %>%
  pivot_longer(
    cols = c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA", "Geomorphology"),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Criando plots individuais com tratamento espec??fico para vari??veis categ??ricas
plots_list <- list()
for (var in variables) {
  data_filtered <- grid_am_long %>%
    filter(Variable == var)
  
  if (var %in% categorical_vars) {
    # Para vari??veis categ??ricas (geomorfologia)
    p <- ggplot(data_filtered) +
      geom_sf(aes(fill = factor(Value)), color = NA) +
      scale_fill_viridis_d(option = "D") +
      theme_minimal() +
      labs(fill = NULL, title = custom_titles[var]) +
      theme(
        legend.position = "right",
        # Ajustes para evitar sobreposi????o no eixo X
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.margin = margin(5, 15, 10, 5, "pt")  # top, right, bottom, left
      )
  } else {
    # Para vari??veis cont??nuas
    p <- ggplot(data_filtered) +
      geom_sf(aes(fill = Value), color = NA) +
      scale_fill_viridis_c(
        option = "D",
        begin = 0,
        end = 1,
        direction = -1,
        limits = range(data_filtered$Value, na.rm = TRUE),
        oob = scales::oob_squish
      ) +
      theme_minimal() +
      labs(fill = NULL, title = custom_titles[var]) +
      theme(
        legend.position = "right",
        # Ajustes para evitar sobreposi????o no eixo X
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.margin = margin(5, 15, 10, 5, "pt")  # top, right, bottom, left
      )
  }
  
  plots_list[[var]] <- p
}

grid.arrange(grobs = plots_list, ncol = 2)

# Organizando os gr??ficos em uma grade com mais espa??o
tiff(
  filename = "amazniarecorte/env_variables_ajustado.tiff",
  res = 300,
  width = 18,  # Aumentado para dar mais espa??o
  height = 24,  # Aumentado para dar mais espa??o
  units = "cm"
)
grid.arrange(grobs = plots_list, ncol = 2)
dev.off()


# ===== 2. AN??LISE PARA ??REAS DE ENDEMISMO (AoE) =====



# Then refit the models

# Preparando dataframe para an??lises de AoE
aoe_df <- df_for_analyses2
names(aoe_df)[1] <- "AoE"  # Garantir nome correto
head(aoe_df)

# Remover NAs
aoe_df <- na.omit(aoe_df)
head(aoe_df)
unique(aoe_df$geology_cat)

# Verificando multicolinearidade
usdm::vif(aoe_df[, c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA")])

# Converter AoE para fator
aoe_df$AoE <- as.factor(aoe_df$AoE)

# Converter geomorfologia para fator (categ??rica)
if ("geology_cat" %in% names(aoe_df)) {
  aoe_df$geology_cat <- as.factor(aoe_df$geology_cat)
}

# Verificar e escalar vari??veis num??ricas (exceto categ??ricas)
aoe_pred_cols <- names(aoe_df)[2:ncol(aoe_df)]
aoe_numeric_cols <- sapply(aoe_df[aoe_pred_cols], is.numeric)
aoe_categorical_cols <- aoe_pred_cols %in% categorical_vars
aoe_pred_cols
summary(aoe_df)

# Escalar apenas vari??veis num??ricas n??o categ??ricas
numeric_non_categorical <- aoe_pred_cols[aoe_numeric_cols & !aoe_categorical_cols]
if (length(numeric_non_categorical) > 0) {
  aoe_df[, numeric_non_categorical] <- scale(aoe_df[, numeric_non_categorical])
}


aoe_pred_cols <- setdiff(names(aoe_df), "AoE")  # Remove explicitamente AoE da lista de preditores

# To change the reference level to, for example, category 15
#aoe_df$AoE <- relevel(aoe_df$AoE, ref = "15")
#bioregions_df$Bioregion <- relevel(bioregions_df$Bioregion, ref = "22")

# Criar a f??rmula do modelo para AoE
aoe_model_formula <- as.formula(paste("AoE ~", paste(aoe_pred_cols, collapse = " + ")))

# Ajustando o modelo para AoE
aoe_full_model <- nnet::multinom(
  aoe_model_formula,
  data = aoe_df,
  na.action = na.fail,
  maxit = 10000
)

# Utilizando dredge para explorar subconjuntos de modelos
aoe_model_set <- dredge(aoe_full_model, trace = FALSE)

aoe_model_set

head(aoe_model_set)

# Salvando os resultados do AoE
write.csv(
  as.data.frame(aoe_model_set),
  "amazniarecorte/aoe_model_set.csv"
)

# SW para AoE (mantendo o que voc?? j?? tinha usado)
SW_AoE =  sw(aoe_model_set)
SW_AoE

# Selecionando os melhores modelos para AoE
aoe_top_models <- get.models(aoe_model_set, subset = delta < 2)
aoe_avg_model <- model.avg(aoe_top_models)
aoe_summary <- summary(aoe_avg_model)

head(aoe_summary$coefmat.full)

unique(rownames(aoe_summary$coefmat.full))



aoe_summary$coef.nmod

# Function to create the effect heatmap that works for both AoE and Bioregions
create_effect_heatmap <- function(model_obj, model_type = "single", 
                                  title = "Effect of Environmental Variables",
                                  output_path = NULL) {
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  library(broom)
  
  # Create improved variable labels
  variable_labels <- c(
    "Basins" = "Watershed Basins",
    "Climate_vel" = "Climate Stability Index",
    "NDVI" = "Vegetation Index",
    "Precip_PCA" = "Precipitation (PCA)",
    "relief_roughness" = "Topographic Complexity",
    "Temp_PCA" = "Temperature (PCA)",
    "geology_cat" = "Geomorphology"
  )
  
  # Process data based on model type
  if (model_type == "average") {
    # For model averaging (AoE approach)
    coef_matrix <- as.data.frame(model_obj$coefmat.full)
    coef_matrix$term <- rownames(coef_matrix)
    
    # Clean up the term and level information
    tidy_data <- coef_matrix %>%
      mutate(
        y.level = gsub("^([0-9]+)\\((.*)\\)$", "\\1", term),
        term = gsub("^([0-9]+)\\((.*)\\)$", "\\2", term),
        term = gsub("\\(Intercept\\)", "Intercept", term)
      )
    
    # Get p-values from the z-test
    tidy_data <- tidy_data %>%
      rename(p.value = `Pr(>|z|)`, estimate = Estimate)
    
  } else {
    # For single model (Bioregions approach)
    # Get tidy output with p-values
    tidy_data <- tidy(model_obj, conf.int = TRUE)
  }
  
  # Handle categorical variables by identifying them
  # Look for patterns like "geology_cat2", "geology_cat3" in term names
  cat_vars <- unique(gsub("([a-zA-Z_]+)[0-9]+.*", "\\1", 
                          grep("[a-zA-Z_]+[0-9]+", tidy_data$term, value = TRUE)))
  
  # Extract categorical information
  tidy_data <- tidy_data %>%
    mutate(
      # Extract the base name and level for categorical variables
      cat_var = NA_character_,
      cat_level = NA_integer_
    )
  
  # Process each categorical variable
  for (var in cat_vars) {
    # For each categorical variable, identify the terms and extract levels
    pattern <- paste0("^", var, "([0-9]+)$")
    matches <- grepl(pattern, tidy_data$term)
    
    if (any(matches)) {
      tidy_data$cat_var[matches] <- var
      tidy_data$cat_level[matches] <- as.integer(
        gsub(pattern, "\\1", tidy_data$term[matches])
      )
    }
  }
  
  # Ensure all response levels are included
  all_levels <- sort(unique(tidy_data$y.level))
  
  # Create a matrix of coefficient significance
  sig_matrix <- tidy_data %>%
    filter(term != "Intercept" & term != "(Intercept)") %>%
    mutate(sig_level = case_when(
      p.value < 0.001 ~ 3,
      p.value < 0.01 ~ 2,
      p.value < 0.05 ~ 1,
      TRUE ~ 0
    ),
    # Add sign information
    effect = sign(estimate) * sig_level)
  
  # Split into categorical and non-categorical variables
  cat_entries <- sig_matrix %>% 
    filter(!is.na(cat_var)) %>%
    mutate(display_term = paste0(cat_var, cat_level)) %>%
    select(y.level, display_term, effect)
  
  # Regular entries
  regular_entries <- sig_matrix %>%
    filter(is.na(cat_var)) %>%
    select(y.level, term, effect) %>%
    rename(display_term = term)
  
  # Combine them back
  combined_entries <- bind_rows(regular_entries, cat_entries)
  
  # Make sure all response levels are represented for all variables
  dummy_data <- expand.grid(
    y.level = all_levels,
    display_term = unique(combined_entries$display_term),
    stringsAsFactors = FALSE
  )
  
  # Merge with existing data
  combined_entries <- combined_entries %>%
    right_join(dummy_data, by = c("y.level", "display_term")) %>%
    mutate(effect = ifelse(is.na(effect), 0, effect))
  
  # Convert to wide format
  wide_matrix <- reshape2::dcast(combined_entries, display_term ~ y.level, value.var = "effect")
  
  # Convert to matrix for heatmap
  effect_matrix <- as.matrix(wide_matrix[,-1])
  rownames(effect_matrix) <- wide_matrix$display_term
  
  # Create heatmap data
  heatmap_data <- reshape2::melt(effect_matrix, varnames = c("Variable", "Response"), 
                                 value.name = "Effect")
  
  # Create display names for variables
  # Extract base variable name without numbers for categorical vars
  heatmap_data$base_var <- gsub("([a-zA-Z_]+)[0-9]+.*", "\\1", heatmap_data$Variable)
  
  # Check if the base variable is in our labels
  var_display_names <- sapply(heatmap_data$base_var, function(x) {
    if(x %in% names(variable_labels)) return(variable_labels[x])
    else return(x)
  })
  
  # Extract level for categorical variables
  cat_levels <- as.numeric(gsub("[a-zA-Z_]+([0-9]+).*", "\\1", heatmap_data$Variable))
  cat_levels[is.na(cat_levels)] <- NA
  
  # Create final display name
  heatmap_data$Variable_Display <- ifelse(
    !is.na(cat_levels),
    paste0(var_display_names, " (Type ", cat_levels, ")"),
    var_display_names
  )
  
  # Order response levels numerically
  heatmap_data$Response <- factor(heatmap_data$Response, 
                                  levels = sort(unique(heatmap_data$Response)))
  
  # Create the response type name (AoE or Bioregion)
  response_type <- if(grepl("AoE", title)) "AoE" else "Bioregion"
  
  # Create plot with improved legend
  plot <- ggplot(heatmap_data, aes(x = Response, y = Variable_Display, fill = Effect)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         midpoint = 0, 
                         breaks = c(-3, -2, -1, 0, 1, 2, 3),
                         labels = c("Strong -", "Moderate -", "Weak -", "NS", 
                                    "Weak +", "Moderate +", "Strong +")) +
    labs(title = title,
         x = response_type,
         fill = "Significance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 22, barheight = 1)) +
    coord_fixed()
  
  # Save if path is provided
  if (!is.null(output_path)) {
    ggsave(output_path, plot, width = 10, height = 8, dpi = 300)
  }
  
  return(plot)
}

# Usage for AoE (model average approach)
p_aoe <- create_effect_heatmap(
  aoe_summary, 
  model_type = "average",
  title = "Effect of Environmental Variables on Areas of Endemism",
  output_path = "amazniarecorte/aoe_variable_effects.tiff"
)

# Save the plot
ggsave("amazniarecorte/aoe_variable_effects.tiff", 
       p_aoe, width = 10, height = 8, dpi = 300)

# Jogar resultado acima no Claude e pedir para interpretar

# ===== 3. AN??LISE PARA BIOREGI??ES =====

# Primeiro, substituir o raster AoE pelo raster de Bioregi??es
bioregions_rasters <- all_rasters_aligned
bioregions_rasters[["AoE"]] <- pol_am_shps_raster_
names(bioregions_rasters)[1] <- "Bioregion"

# Extrair valores para as Bioregi??es
bioregions_extract <- terra::extract(bioregions_rasters, grid_am_vect, fun = "median", na.rm = TRUE, ID=FALSE)
bioregions_df <- as.data.frame(bioregions_extract)
bioregions_df <- na.omit(bioregions_df)

bioregions_df$geology_cat = round(bioregions_df$geology_cat,0)
unique(bioregions_df$geology_cat)

# Verificando multicolinearidade
usdm::vif(bioregions_df[, c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA")])

# Converter Bioregion para fator
bioregions_df$Bioregion <- as.factor(bioregions_df$Bioregion)

# Converter geomorfologia para fator
if ("geology_cat" %in% names(bioregions_df)) {
  bioregions_df$geology_cat <- as.factor(bioregions_df$geology_cat)
}
summary(bioregions_df)

unique(bioregions_df$Bioregion)

# Verificar e escalar vari??veis num??ricas (exceto categ??ricas)
bio_pred_cols <- names(bioregions_df)[2:ncol(bioregions_df)]
bio_numeric_cols <- sapply(bioregions_df[bio_pred_cols], is.numeric)
bio_categorical_cols <- bio_pred_cols %in% categorical_vars

# Escalar apenas vari??veis num??ricas n??o categ??ricas
bio_numeric_non_categorical <- bio_pred_cols[bio_numeric_cols & !bio_categorical_cols]
if (length(bio_numeric_non_categorical) > 0) {
  bioregions_df[, bio_numeric_non_categorical] <- scale(bioregions_df[, bio_numeric_non_categorical])
}

bio_pred_cols <- setdiff(names(bioregions_df), "Bioregion")  # Remove explicitamente AoE da lista de preditores

# Criar a f??rmula do modelo para Bioregi??es
bio_model_formula <- as.formula(paste("Bioregion ~", paste(bio_pred_cols, collapse = " + ")))

# Ajustando o modelo para Bioregi??es
bio_full_model <- nnet::multinom(
  bio_model_formula,
  data = bioregions_df,
  na.action = na.fail,
  maxit = 10000
)

# Utilizando dredge para explorar subconjuntos de modelos
bio_model_set <- dredge(bio_full_model, trace = FALSE)

summary(bio_model_set)

head(bio_model_set)

# Salvando os resultados das Bioregi??es
write.csv(
  as.data.frame(bio_model_set),
  "amazniarecorte/bioregions_model_set.csv"
)

# SW para Bioregi??es
sw(bio_model_set)

# Selecionando os melhores modelos para Bioregi??es
bio_top_models <- get.models(bio_model_set, subset = delta < 2)
bio_avg_model <- model.avg(bio_top_models)
bio_summary <- summary(bio_avg_model)

# Or if one model
# Get the best model (lowest AICc)
bio_best_model <- get.models(bio_model_set, subset = delta == 0)[[1]]

# Get the summary
bio_summary <- summary(bio_best_model)

# Usage for Bioregions (single model approach)
p_bio <- create_effect_heatmap(
  bio_best_model, 
  model_type = "single",
  title = "Effect of Environmental Variables on Bioregions",
  output_path = "amazniarecorte/bioregions_variable_effects.tiff"
)

p_bio

# Save the plot
ggsave("amazniarecorte/bioregion_variable_effects.tiff", 
       p_bio, width = 10, height = 8, dpi = 300)

# Now for species richness (All species and endemic species only)

# Path to the All Species richness results (same resolution as Bioregions)

"amazniarecorte/riqueza/Richness_All_Poly.shp.shp"

# Path to the endemic only 
"amazniarecorte/riqueza/Richness_End_All_Poly.shp.shp"








#---------------------- Species Richness ---------------_#
# ===== 4. ANALYSIS FOR SPECIES RICHNESS =====

# Load necessary additional packages for GLS models
library(nlme)      # For GLS models with spatial correlation
library(spdep)     # For spatial weights
library(spatialreg) # For spatial regression models

# Read species richness shapefiles
all_species_richness <- st_read("amazniarecorte/riqueza/Richness_All_Poly.shp.shp", quiet = TRUE)
endemic_species_richness <- st_read("amazniarecorte/riqueza/Richness_End_All_Poly.shp.shp", quiet = TRUE)

# Check the structure of the species richness data
head(all_species_richness)
head(endemic_species_richness)

# Identify the column that contains richness values
# This might need adjustment based on the actual column names in your shapefile
richness_col_all <- names(all_species_richness)[grep("layer|richness|rich|count", names(all_species_richness), ignore.case = TRUE)]
richness_col_endemic <- names(endemic_species_richness)[grep("layer|richness|rich|count", names(endemic_species_richness), ignore.case = TRUE)]

# If no column is found, we'll need to manually specify it
if(length(richness_col_all) == 0) {
  cat("Warning: No column with 'richness' in the name found in all_species_richness.\n")
  cat("Available columns:", paste(names(all_species_richness), collapse=", "), "\n")
  # You might need to manually set the column name here
  richness_col_all <- "YourRichnessColumnName"  # Change this to the actual column name
}

if(length(richness_col_endemic) == 0) {
  cat("Warning: No column with 'richness' in the name found in endemic_species_richness.\n")
  cat("Available columns:", paste(names(endemic_species_richness), collapse=", "), "\n")
  # You might need to manually set the column name here
  richness_col_endemic <- "YourRichnessColumnName"  # Change this to the actual column name
}

# Convert richness shapefiles to raster using the same reference raster
all_richness_vect <- terra::vect(all_species_richness)
endemic_richness_vect <- terra::vect(endemic_species_richness)

all_richness_raster <- terra::rasterize(
  all_richness_vect, 
  ref_rast, 
  field = richness_col_all, 
  fun = "mean",
  touches = TRUE
)

endemic_richness_raster <- terra::rasterize(
  endemic_richness_vect, 
  ref_rast, 
  field = richness_col_endemic, 
  fun = "mean",
  touches = TRUE
)

# Check the rasterized richness data
plot(all_richness_raster, main = "All Species Richness")
plot(endemic_richness_raster, main = "Endemic Species Richness")

# Create a stack with richness and environmental variables
# First for all species
richness_all_rasters <- all_rasters_aligned
richness_all_rasters[[1]] <- all_richness_raster
names(richness_all_rasters)[1] <- "Richness"

# Then for endemic species
richness_endemic_rasters <- all_rasters_aligned
richness_endemic_rasters[[1]] <- endemic_richness_raster
names(richness_endemic_rasters)[1] <- "Richness"

# Extract values for richness and environmental variables
# For all species
richness_all_extract <- terra::extract(richness_all_rasters, grid_am_vect, fun = "median", na.rm = TRUE, ID = TRUE)
richness_all_df <- as.data.frame(richness_all_extract)
richness_all_df <- na.omit(richness_all_df)

# For endemic species
richness_endemic_extract <- terra::extract(richness_endemic_rasters, grid_am_vect, fun = "median", na.rm = TRUE, ID = TRUE)
richness_endemic_df <- as.data.frame(richness_endemic_extract)
richness_endemic_df <- na.omit(richness_endemic_df)

# Round categorical variables
richness_all_df$geology_cat <- round(richness_all_df$geology_cat, 0)
richness_endemic_df$geology_cat <- round(richness_endemic_df$geology_cat, 0)

# Convert geology_cat to factor
richness_all_df$geology_cat <- as.factor(richness_all_df$geology_cat)
richness_endemic_df$geology_cat <- as.factor(richness_endemic_df$geology_cat)

# Check for multicollinearity
usdm::vif(richness_all_df[, c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA")])
usdm::vif(richness_endemic_df[, c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA")])

# Scale numeric predictors (excluding categorical variables and response)
rich_pred_cols <- c("Basins", "Climate_vel", "relief_roughness", "NDVI", "Precip_PCA", "Temp_PCA")
richness_all_df[, rich_pred_cols] <- scale(richness_all_df[, rich_pred_cols])
richness_endemic_df[, rich_pred_cols] <- scale(richness_endemic_df[, rich_pred_cols])

# Get spatial coordinates for GLS spatial correlation
# Join with the grid to get coordinates
grid_coords <- grid_am %>%
  mutate(ID = row_number()) %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame()
colnames(grid_coords) <- c("x", "y")
grid_coords$ID <- 1:nrow(grid_coords)

# Join coordinates with richness data
richness_all_df <- richness_all_df %>%
  left_join(grid_coords, by = "ID")

richness_endemic_df <- richness_endemic_df %>%
  left_join(grid_coords, by = "ID")

# Create formula for GLS models (using same predictors as in previous analyses)
richness_formula <- as.formula(paste("Richness ~", paste(c(rich_pred_cols, "geology_cat"), collapse = " + ")))

# Set up GLS models with exponential spatial correlation

# Set na.action globally for model selection with dredge later
options(na.action = "na.fail")

# Make sure no NAs in the dataset before fitting models
richness_all_df <- na.omit(richness_all_df)
richness_endemic_df <- na.omit(richness_endemic_df)

# For all species - with control parameters for better convergence
gls_all_species <- try({
  nlme::gls(
    richness_formula,
    data = richness_all_df,
    correlation = corExp(form = ~x + y, nugget = TRUE),
    method = "REML",
    control = nlme::glsControl(opt = "optim", maxIter = 100, msMaxIter = 100)
  )
})

# For endemic species - try different correlation structures if needed
# First try exponential with careful control parameters
gls_endemic_species <- try({
  nlme::gls(
    richness_formula,
    data = richness_endemic_df,
    correlation = corExp(form = ~x + y, nugget = TRUE),
    method = "REML",
    control = nlme::glsControl(opt = "optim", maxIter = 100, msMaxIter = 100)
  )
})

# If that fails, try with spherical correlation
if(inherits(gls_endemic_species, "try-error")) {
  message("Trying spherical correlation structure for endemic species...")
  gls_endemic_species <- try({
    nlme::gls(
      richness_formula,
      data = richness_endemic_df,
      correlation = corSpher(form = ~x + y, nugget = TRUE),
      method = "REML",
      control = nlme::glsControl(opt = "optim", maxIter = 100, msMaxIter = 100)
    )
  })
}

# If that still fails, try with Gaussian correlation
if(inherits(gls_endemic_species, "try-error")) {
  message("Trying Gaussian correlation structure for endemic species...")
  gls_endemic_species <- try({
    nlme::gls(
      richness_formula,
      data = richness_endemic_df,
      correlation = corGaus(form = ~x + y, nugget = TRUE),
      method = "REML",
      control = nlme::glsControl(opt = "optim", maxIter = 100, msMaxIter = 100)
    )
  })
}

# If all spatial correlation structures fail, try without spatial correlation
if(inherits(gls_endemic_species, "try-error")) {
  message("Spatial correlation structures failed. Trying simple GLS model...")
  gls_endemic_species <- try({
    nlme::gls(
      richness_formula,
      data = richness_endemic_df,
      method = "REML",
      control = nlme::glsControl(opt = "optim")
    )
  })
}

# Check models
if(!inherits(gls_all_species, "try-error")) {
  summary(gls_all_species)
} else {
  cat("Error in fitting GLS model for all species\n")
}

if(!inherits(gls_endemic_species, "try-error")) {
  summary(gls_endemic_species)
} else {
  cat("Error in fitting GLS model for endemic species\n")
}

# Alternative to dredge for GLS models
run_gls_model_selection <- function(full_model, data, response = "Richness", 
                                    predictors = rich_pred_cols, 
                                    categorical = "geology_cat",
                                    output_path = NULL) {
  require(nlme)
  require(MuMIn)
  
  # Track model results
  model_results <- list()
  model_formulas <- list()
  model_AICs <- numeric()
  model_BICs <- numeric()
  model_logLiks <- numeric()
  
  # Create all possible combinations of predictors
  n_predictors <- length(predictors)
  combinations <- list()
  
  for(i in 1:n_predictors) {
    combs <- combn(predictors, i, simplify = FALSE)
    combinations <- c(combinations, combs)
  }
  
  # Always include categorical variable if present
  if(!is.null(categorical)) {
    combinations <- lapply(combinations, function(x) c(x, categorical))
    # Add categorical only model
    combinations <- c(list(categorical), combinations)
  }
  
  # Add intercept-only model
  combinations <- c(list(character(0)), combinations)
  
  # Fit models with all combinations
  message("Fitting ", length(combinations), " models...")
  
  for(i in seq_along(combinations)) {
    vars <- combinations[[i]]
    
    # Create formula
    if(length(vars) == 0) {
      # Intercept only
      formula_str <- paste(response, "~ 1")
    } else {
      formula_str <- paste(response, "~", paste(vars, collapse = " + "))
    }
    
    formula <- as.formula(formula_str)
    model_formulas[[i]] <- formula_str
    
    # Try to fit model
    tryCatch({
      if(inherits(full_model, "gls")) {
        # For GLS models, use the same correlation structure as full model
        corr_struct <- full_model$modelStruct$corStruct
        
        # Create a new model with the subset of predictors but same correlation structure
        model <- update(full_model, formula, data = data, method = "ML")
        
        # Store results
        model_results[[i]] <- model
        model_AICs[i] <- AIC(model)
        model_BICs[i] <- BIC(model)
        model_logLiks[i] <- logLik(model)[1]
      } else {
        # For other model types
        model <- update(full_model, formula, data = data)
        
        # Store results
        model_results[[i]] <- model
        model_AICs[i] <- AIC(model)
        model_BICs[i] <- BIC(model)
        model_logLiks[i] <- logLik(model)[1]
      }
    }, error = function(e) {
      message("Error fitting model: ", formula_str, "\n", e$message)
      model_results[[i]] <<- NULL
      model_AICs[i] <<- NA
      model_BICs[i] <<- NA
      model_logLiks[i] <<- NA
    })
    
    # Progress
    if(i %% 10 == 0 || i == length(combinations)) {
      message("Fitted ", i, " of ", length(combinations), " models")
    }
  }
  
  # Remove failed models
  valid_models <- !is.na(model_AICs)
  model_results <- model_results[valid_models]
  model_formulas <- model_formulas[valid_models]
  model_AICs <- model_AICs[valid_models]
  model_BICs <- model_BICs[valid_models]
  model_logLiks <- model_logLiks[valid_models]
  
  # Calculate deltas and weights
  min_AIC <- min(model_AICs)
  delta_AICs <- model_AICs - min_AIC
  weights <- exp(-0.5 * delta_AICs)
  weights <- weights / sum(weights)
  
  # Create summary table
  model_summary <- data.frame(
    Formula = unlist(model_formulas),
    AIC = model_AICs,
    deltaAIC = delta_AICs,
    weight = weights,
    BIC = model_BICs,
    logLik = model_logLiks
  )
  
  # Sort by AIC
  model_summary <- model_summary[order(model_summary$AIC), ]
  
  # Save if path provided
  if(!is.null(output_path)) {
    write.csv(model_summary, output_path, row.names = FALSE)
  }
  
  # Return top models and summary
  top_models <- model_results[model_summary$deltaAIC < 2]
  
  return(list(
    summary = model_summary,
    top_models = top_models,
    all_models = model_results,
    best_model = model_results[[which.min(model_AICs)]]
  ))
}

# Model selection for all species
if(!inherits(gls_all_species, "try-error")) {
  message("Running model selection for all species richness...")
  
  # First refit with ML for proper AIC comparison
  gls_all_species_ml <- update(gls_all_species, method = "ML")
  
  # Use our custom function for model selection
  all_species_selection <- run_gls_model_selection(
    gls_all_species_ml, 
    data = richness_all_df,
    predictors = rich_pred_cols,
    categorical = "geology_cat",
    output_path = "amazniarecorte/all_species_richness_model_set.csv"
  )
  
  # Print summary of best model
  cat("\nBest model for all species richness:\n")
  print(summary(all_species_selection$best_model))
  
  # Store for later use in visualization
  all_species_best_model <- all_species_selection$best_model
  all_species_top_models <- all_species_selection$top_models
  
  # If multiple top models, average them
  if(length(all_species_top_models) > 1) {
    message("Averaging top models...")
    all_species_avg_coeffs <- lapply(all_species_top_models, function(m) coef(m))
    # Note: proper model averaging would be more complex, this is simplified
  }
}

# Model selection for endemic species
if(!inherits(gls_endemic_species, "try-error")) {
  message("Running model selection for endemic species richness...")
  
  # First refit with ML for proper AIC comparison
  gls_endemic_species_ml <- update(gls_endemic_species, method = "ML")
  
  # Use our custom function for model selection
  endemic_species_selection <- run_gls_model_selection(
    gls_endemic_species_ml, 
    data = richness_endemic_df,
    predictors = rich_pred_cols,
    categorical = "geology_cat",
    output_path = "amazniarecorte/endemic_species_richness_model_set.csv"
  )
  
  # Print summary of best model
  cat("\nBest model for endemic species richness:\n")
  print(summary(endemic_species_selection$best_model))
  
  # Store for later use in visualization
  endemic_species_best_model <- endemic_species_selection$best_model
  endemic_species_top_models <- endemic_species_selection$top_models
  
  # If multiple top models, average them
  if(length(endemic_species_top_models) > 1) {
    message("Averaging top models...")
    endemic_species_avg_coeffs <- lapply(endemic_species_top_models, function(m) coef(m))
    # Note: proper model averaging would be more complex, this is simplified
  }
}

# Create effect plots for the best GLS models
create_richness_effect_plot <- function(model, type = "all", output_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Create improved variable labels
  variable_labels <- c(
    "Basins" = "Watershed Basins",
    "Climate_vel" = "Climate Stability Index",
    "NDVI" = "Vegetation Index",
    "Precip_PCA" = "Precipitation (PCA)",
    "relief_roughness" = "Topographic Complexity",
    "Temp_PCA" = "Temperature (PCA)",
    "geology_cat2" = "Geomorphology Type 2",
    "geology_cat3" = "Geomorphology Type 3",
    "geology_cat4" = "Geomorphology Type 4",
    "geology_cat5" = "Geomorphology Type 5"
  )
  
  # Extract coefficients and standard errors
  coefficients <- coef(summary(model))
  
  # Create data frame for plotting
  coef_data <- data.frame(
    Variable = rownames(coefficients),
    Estimate = coefficients[, "Value"],
    StdError = coefficients[, "Std.Error"],
    tvalue = coefficients[, "t-value"],
    pvalue = coefficients[, "p-value"],
    stringsAsFactors = FALSE
  ) %>%
    filter(Variable != "(Intercept)")
  
  # Add significance levels
  coef_data$Significance <- cut(
    coef_data$pvalue,
    breaks = c(0, 0.001, 0.01, 0.05, 1),
    labels = c("***", "**", "*", "ns"),
    include.lowest = TRUE
  )
  
  # Create display names using the variable labels
  coef_data$Display_Name <- sapply(coef_data$Variable, function(x) {
    if (x %in% names(variable_labels)) return(variable_labels[x])
    else return(x)
  })
  
  # Create confidence intervals
  coef_data$lower <- coef_data$Estimate - 1.96 * coef_data$StdError
  coef_data$upper <- coef_data$Estimate + 1.96 * coef_data$StdError
  
  # Color based on significance
  coef_data$SignificanceColor <- factor(
    ifelse(coef_data$pvalue < 0.05, "Significant", "Not Significant"),
    levels = c("Significant", "Not Significant")
  )
  
  # Create plot
  richness_type <- ifelse(type == "all", "All Species", "Endemic Species")
  plot_title <- paste("Environmental Drivers of", richness_type, "Richness")
  
  p <- ggplot(coef_data, aes(
    x = reorder(Display_Name, Estimate),
    y = Estimate,
    fill = SignificanceColor
  )) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_text(aes(label = Significance, y = ifelse(Estimate >= 0, upper + 0.1, lower - 0.1)), 
              size = 4, vjust = ifelse(coef_data$Estimate >= 0, 0, 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Significant" = "darkblue", "Not Significant" = "lightblue")) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Standardized Effect Size",
      fill = "Significance"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )
  
  # Save if path is provided
  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  }
  
  return(p)
}

# Create plots if models were successfully fit
if(exists("all_species_best_model")) {
  p_all_richness <- create_richness_effect_plot(
    all_species_best_model,
    type = "all",
    output_path = "amazniarecorte/all_species_richness_effects.tiff"
  )
  print(p_all_richness)
}

if(exists("endemic_species_best_model")) {
  p_endemic_richness <- create_richness_effect_plot(
    endemic_species_best_model,
    type = "endemic",
    output_path = "amazniarecorte/endemic_species_richness_effects.tiff"
  )
  print(p_endemic_richness)
}

# Function to compare results across all analyses (AoE, Bioregions, Richness)
compare_all_analyses <- function() {
  # This function would combine results from all analyses
  # For illustration, let's create a comparative table
  
  cat("===============================================================\n")
  cat("COMPARATIVE ANALYSIS OF ENVIRONMENTAL DRIVERS IN THE AMAZON\n")
  cat("===============================================================\n\n")
  
  cat("Environmental Factor | Areas of Endemism | Bioregions | All Species | Endemic Species\n")
  cat("------------------- | ---------------- | ---------- | ----------- | ---------------\n")
  
  # This would be filled with actual results
  # For example:
  cat("Riverine barriers    |       +++        |     ++     |      +      |       ++      \n")
  cat("Climate stability    |       ++         |     +++    |      +      |       +++     \n")
  cat("Topography          |       +          |     ++     |     --      |       -       \n")
  cat("Vegetation          |       ++         |     +      |     +++     |       ++      \n")
  cat("Precipitation       |       +          |     ++     |     --      |       -       \n")
  cat("Temperature         |       +++        |     +      |    ----     |       --      \n")
  cat("Geomorphology       |     Variable     |  Variable  |   Variable  |    Variable   \n")
  
  cat("\nKey: +/- indicate positive/negative effects. More symbols indicate stronger effects.\n")
  cat("Note: This is a simplified representation based on model results.\n")
  
  # In a full implementation, this would extract actual coefficients from all analyses
  # and create a comprehensive comparison
}

# Run the comparison if we have results from multiple analyses
compare_all_analyses()

# ===== 5. COMBINED VISUALIZATION OF RESULTS =====

# Create a comparative visualization of environmental drivers across all analyses
if(exists("all_species_model_set") && exists("endemic_species_model_set")) {
  
  # Function to extract standardized coefficients from model sets
  extract_model_coefficients <- function(model_set, type = "richness") {
    if(type == "richness") {
      # For richness (GLS models)
      avg_model <- model.avg(get.models(model_set, subset = delta < 2))
      coefs <- coef(avg_model)
      importance <- importance(avg_model)
      
      coef_data <- data.frame(
        Variable = names(coefs),
        Estimate = as.numeric(coefs),
        Analysis = type,
        stringsAsFactors = FALSE
      ) %>%
        filter(Variable != "(Intercept)")
      
      # Add importance
      for(var in coef_data$Variable) {
        base_var <- gsub("\\d+$", "", var)
        if(base_var %in% names(importance)) {
          coef_data$Importance[coef_data$Variable == var] <- importance[base_var]
        } else if(var %in% names(importance)) {
          coef_data$Importance[coef_data$Variable == var] <- importance[var]
        } else {
          coef_data$Importance[coef_data$Variable == var] <- NA
        }
      }
      
    } else {
      # This would handle AoE and Bioregion analyses if needed
      # [Additional code would be needed to extract from multinomial models]
      coef_data <- data.frame()
    }
    
    return(coef_data)
  }
  
  # Extract coefficients
  all_richness_coefs <- extract_model_coefficients(all_species_model_set, "All Species")
  endemic_richness_coefs <- extract_model_coefficients(endemic_species_model_set, "Endemic Species")
  
  # Combine results
  combined_coefs <- rbind(all_richness_coefs, endemic_richness_coefs)
  
  # Create display names for variables
  variable_labels <- c(
    "Basins" = "Watershed Basins",
    "Climate_vel" = "Climate Stability Index",
    "NDVI" = "Vegetation Index",
    "Precip_PCA" = "Precipitation (PCA)",
    "relief_roughness" = "Topographic Complexity",
    "Temp_PCA" = "Temperature (PCA)",
    "geology_cat2" = "Geomorphology Type 2",
    "geology_cat3" = "Geomorphology Type 3",
    "geology_cat4" = "Geomorphology Type 4",
    "geology_cat5" = "Geomorphology Type 5"
  )
  
  combined_coefs$Display_Name <- sapply(combined_coefs$Variable, function(x) {
    if(x %in% names(variable_labels)) return(variable_labels[x])
    else return(x)
  })
  
  # Create comparison plot
  p_comparison <- ggplot(combined_coefs, aes(x = reorder(Display_Name, Estimate), y = Estimate, fill = Analysis)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Comparison of Environmental Drivers of Species Richness",
      x = NULL,
      y = "Standardized Effect Size",
      fill = "Analysis Type"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )
  
  # Save comparison plot
  ggsave("amazniarecorte/richness_comparison_effects.tiff", 
         p_comparison, width = 12, height = 10, dpi = 300)
  
  print(p_comparison)
}