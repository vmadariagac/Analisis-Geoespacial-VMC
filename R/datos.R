library(sf)
library(dplyr)
library(spdep)
library(spatialreg)
library(CARBayes)

ruta_gpkg <- "G:/Mi unidad/Universidad/POSGRADO/Analisis Geoespacial/Analisis-Geoespacial-VMC/R/Barrios_var.gpkg"

gdf <- st_read(
  ruta_gpkg,
  layer = "variables"
)

dim(gdf)
names(gdf)
head(gdf)
st_crs(gdf)

nrow(gdf)
length(unique(gdf$nombre_bar))

library(sf)
library(spdep)

vecinos_queen <- poly2nb(
  gdf,
  queen = TRUE,
  row.names = gdf$nombre_bar
)

vecinos_queen <- poly2nb(
  gdf,
  queen = TRUE,
  row.names = gdf$nombre_bar,
  snap = 1
)

summary(vecinos_queen)
card(vecinos_queen)
islas <- which(card(vecinos_queen) == 0)
islas

# Para Moran, SAR y SEM
listw_queen <- nb2listw(
  vecinos_queen,
  style = "W",
  zero.policy = FALSE)

# Para CARBayes
W_queen <- nb2mat(
  vecinos_queen,
  style = "B",
  zero.policy = FALSE)
summary(listw_queen)

dim(W_queen)
rowSums(W_queen)[1:10]

gdf$Densidad <- gdf$Movimientos / gdf$Area_km2
gdf$Log_Densidad <- log1p(gdf$Densidad)

hist(
    gdf$Log_Densidad,
    breaks=20,
    col="lightblue",
    main="Log densidad de movimientos",
    xlab="log(1+densidad)")

#OLS
modelo_lm <- lm(
    Log_Densidad ~
        Porc_Artificial +
        Porc_Qca +
        Porc_TPu +
        Pend_mean,
    data=gdf)
summary(modelo_lm)

#Residuos OLS
gdf$residuos <- residuals(modelo_lm)
moran.test(gdf$residuos, listw_queen)

lm.RStests( modelo_lm, listw = listw_queen)
print(resultado_rs)

resultado_rs <- lm.RStests(
  modelo_lm,
  listw = listw_queen,
  test = "all")
summary(resultado_rs)

#SAR
modelo_sar <- lagsarlm(
    Log_Densidad ~
        Porc_Artificial +
        Porc_Qca +
        Porc_TPu +
        Pend_mean,
    data = gdf,
    listw = listw_queen)
    summary(modelo_sar)

library(spatialreg)

impactos <- impacts(modelo_sar, listw = listw_queen, R = 1000)
summary(
    impactos,
    zstats = TRUE,
    short = TRUE)

#SEM
modelo_sem <- errorsarlm(
  Log_Densidad ~
    Porc_Artificial +
    Porc_Qca +
    Porc_TPu +
    Pend_mean,
  data = gdf,
  listw = listw_queen)
  summary(modelo_sem)

#Residuo SEM
  res_sem <- residuals(modelo_sem)
moran.test(
  res_sem,
  listw_queen)

#Comparación

AIC(modelo_lm)
AIC(modelo_sem)
AIC(modelo_sar)

logLik(modelo_lm)
logLik(modelo_sem)
logLik(modelo_sar)

#CAR
summary(gdf$Movimientos)
summary(gdf$Area_km2)

sum(is.na(gdf$Movimientos))
sum(is.na(gdf$Area_km2))

sum(gdf$Movimientos < 0)
sum(gdf$Area_km2 <= 0)

var(gdf$Movimientos)
mean(gdf$Movimientos)

#Matriz
dim(W_queen)

all(W_queen == t(W_queen))
all(diag(W_queen) == 0)

range(rowSums(W_queen))
sum(rowSums(W_queen) == 0)

#Modelo de referencia
modelo_poisson <- glm(
  Movimientos ~
    Porc_Artificial +
    Porc_Qca +
    Porc_TPu +
    Pend_mean +
    offset(log(Area_km2)),
  family = poisson(link = "log"),
  data = gdf)

summary(modelo_poisson)
AIC(modelo_poisson)

dispersion_poisson <- sum(
  residuals(modelo_poisson, type = "pearson")^2
) / df.residual(modelo_poisson)

dispersion_poisson

#CAR Leroux
library(CARBayes)

set.seed(42)

modelo_car_poisson_prueba <- S.CARleroux(
  formula =
    Movimientos ~
      Porc_Artificial +
      Porc_Qca +
      Porc_TPu +
      Pend_mean +
      offset(log(Area_km2)),
  family = "poisson",
  data = st_drop_geometry(gdf),
  W = W_queen,
  burnin = 50000,
  n.sample = 250000,
  thin = 5,
  n.chains = 1)

  print(modelo_car_poisson_prueba)
summary(modelo_car_poisson_prueba)

names(modelo_car_poisson_prueba)
str(modelo_car_poisson_prueba, max.level = 1)
names(modelo_car_poisson_prueba$samples)
dim(modelo_car_poisson_prueba$samples$phi)
phi_media <- apply(
  modelo_car_poisson_prueba$samples$phi,
  2,
  mean)
gdf$Efecto_CAR <- phi_media
summary(gdf$Efecto_CAR)
head(gdf[, c("nombre_bar", "Efecto_CAR")])

library(ggplot2)

ggplot(gdf) +
  geom_sf(
    aes(fill = Efecto_CAR),
    color = "black",
    linewidth = 0.1
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    name = "Efecto CAR"
  ) +
  labs(
    title = "Efecto espacial residual del modelo CAR Poisson"
  ) +
  theme_void()

"INLA" %in% rownames(installed.packages())
install.packages(
    "INLA",
    repos = c(
        getOption("repos"),
        INLA = "https://inla.r-inla-download.org/R/stable"
    ),
    dep = TRUE
)