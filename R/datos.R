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


#Car BN
library(CARBayes)

modelo_car_poisson <- S.CARleroux(
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

print(modelo_car_poisson)
summary(modelo_car_poisson)

modelo_car_poisson$summary.results
pred_car <- modelo_car_poisson$fitted.values

length(pred_car)
length(gdf$Movimientos)

print(pred_car)

rmse_car <- sqrt(mean((gdf$Movimientos - pred_car)^2))

mae_car <- mean(abs(gdf$Movimientos - pred_car))

print(rmse_car)
print(mae_car)

y_obs <- gdf$Movimientos
r2_car <- 1 - sum((y_obs - pred_car)^2) / sum((y_obs - mean(y_obs))^2)

print(r2_car)
logLik(modelo_car_poisson)

#
#
#SAR
pred_sar <- modelo_sar$fitted.values
rmse_sar <- sqrt(mean((gdf$Log_Densidad - pred_sar)^2))

mae_sar <- mean(abs(gdf$Log_Densidad - pred_sar))

print(pred_sar)
print(rmse_sar)
print(mae_sar)





#
#
#CAR

install.packages(
  "INLA",
  repos = c(
    "https://inla.r-inla-download.org/R/stable",
    "https://cloud.r-project.org"
  ),
  dep = TRUE)

  install.packages(
  "INLA",
  repos = c(
    "https://inla.r-inla-download.org/R/testing",
    "https://cloud.r-project.org"
  ),
  dep=TRUE)

library(INLA)

inla.read.graph

"INLA" %in% rownames(installed.packages())
datos <- st_drop_geometry(gdf)
datos$ID <- 1:nrow(datos)
datos$log_area <- log(datos$Area_km2)

#CAR BYM BN
nrow(datos)
dim(W_queen)

library(spdep)
nb <- mat2listw(
  W_queen,
  style="B"
)$neighbours

nb2INLA(
  "barrios.adj", nb)
grafo <- inla.read.graph("barrios.adj")

inla.read.graph
grafo$n

formula_nb <- Movimientos ~
  Porc_Artificial +
  Porc_Qca +
  Porc_TPu +
  Pend_mean +
  offset(log_area)

#BN GLM con INLA
modelo_nb <- inla(
  formula_nb,
  family = "nbinomial",
  data = datos,
  control.predictor = list(
    compute = TRUE),

  control.compute = list(
    dic = TRUE,
    waic = TRUE))


modelo_nb$dic$dic
modelo_nb$waic$waic

formula_nb_esp <- Movimientos ~
  Porc_Artificial +
  Porc_Qca +
  Porc_TPu +
  Pend_mean +
  offset(log_area) +
  f(
    ID,
    model = "bym2",
    graph = grafo
  )

#ModeloBYM BN con INLA
  modelo_nb_esp <- inla(
  formula_nb_esp,
  family = "nbinomial",
  data = datos,
  control.predictor = list(
    compute = TRUE
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE
  ))

print(modelo_nb_esp)
summary(modelo_nb_esp)

modelo_nb_esp$dic$dic
modelo_nb_esp$waic$waic

pred_nb_esp <- modelo_nb_esp$summary.fitted.values$mean

rmse_nb_esp <- sqrt(mean(
  (datos$Movimientos - pred_nb_esp)^2
))

mae_nb_esp <- mean(
  abs(datos$Movimientos - pred_nb_esp)
)

print(rmse_nb_esp)
print(mae_nb_esp)

modelo_nb_esp$summary.fixed
modelo_nb_esp$summary.hyperpar

pred_nb_esp <- modelo_nb_esp$summary.fitted.values$mean

rmse_nb_esp <- sqrt(mean(
  (datos$Movimientos - pred_nb_esp)^2
))

mae_nb_esp <- mean(
  abs(datos$Movimientos - pred_nb_esp)
)

r2_nb_esp <- 1 -
sum((datos$Movimientos-pred_nb_esp)^2)/
sum((datos$Movimientos-mean(datos$Movimientos))^2)

print(rmse_nb_esp)
print(mae_nb_esp)
print(r2_nb_esp)

efecto_espacial <- modelo_nb_esp$summary.random$ID$mean[1:348]

gdf$efecto_espacial <- efecto_espacial
plot(
  gdf["efecto_espacial"]
)

pred_nb_esp <- modelo_nb_esp$summary.fitted.values$mean

gdf$pred_nb_esp <- pred_nb_esp

library(RColorBrewer)

plot(
  gdf["pred_nb_esp"],
  col = brewer.pal(9, "YlOrRd"),
  border = "white",
  main = "Predicción Binomial Negativa espacial (INLA)"
)

library(ggplot2)

ggplot(gdf) +
  geom_sf(
    aes(fill = pred_nb_esp),
    color = "black",
    linewidth = 0.1
  ) +
  scale_fill_gradientn(
  colors = brewer.pal(9,"YlOrRd"),
  limits = c(0,35),
  name="Movimientos esperados BYM BN",
  guide = guide_colorbar(
    barwidth = 1,
    barheight = 15)
) +
  theme_void()

library(classInt)

breaks_pred <- classIntervals(
  gdf$pred_nb_esp,
  n = 5,
  style = "quantile"
)$brks

breaks_pred

gdf$pred_cat <- cut(
  gdf$pred_nb_esp,
  breaks = breaks_pred,
  include.lowest = TRUE
)

install.packages("ggspatial")

library(sf)
library(ggspatial)
library(RColorBrewer)
library(grid)

ggplot(gdf) +

  geom_sf(
    aes(fill = pred_cat),
    color = "black",
    linewidth = 0.1
  ) +

  scale_fill_brewer(
    palette = "YlOrRd",
    name = "Movimientos\nesperados"
  ) +
  labs(
    title = "Predicción del modelo Binomial Negativo CAR BYM"
  ) +

  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering()
  ) +

  annotation_scale(
    location = "bl",
    width_hint = 0.25
  ) +

  coord_sf(
    datum = st_crs(gdf)
  ) +

  theme_bw() +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.text = element_text(
      size = 8
    ),

    axis.title = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      angle = 90,
      hjust = 0.5
    ),

    legend.key.height = unit(1.5, "cm")
  )

modelo_nb_esp$summary.fitted.values$sd

gdf$incertidumbre <- modelo_nb_esp$summary.fitted.values$sd

summary(gdf$incertidumbre)

library(ggplot2)
library(ggspatial)
library(RColorBrewer)
library(grid)
install.packages("viridis")
library(viridis)

ggplot(gdf) +

  geom_sf(
    aes(fill = incertidumbre),
    color = "black",
    linewidth = 0.1
  ) +

  scale_fill_gradientn(
    colors = viridisLite::viridis(9),
    name = "Desviación estándar predictiva con BYM BN",
    guide = guide_colorbar(
      title.position = "right",
      barheight = unit(5, "cm"),
      barwidth = unit(0.5, "cm"))
  ) +

   labs(
    title = "Incertidumbre de la predicción del modelo Binomial Negativo BYM"
  ) +

  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering()
  ) +

  annotation_scale(
    location = "bl",
    width_hint = 0.25
  ) +

  coord_sf(
    datum = st_crs(gdf)
  ) +

  scale_x_continuous(
    breaks = seq(
      4699257,
      4726591,
      by = 5000
    )
  ) +

  scale_y_continuous(
    breaks = seq(
      2239727,
      2263174,
      by = 5000
    )
  ) +

  theme_bw() +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.text = element_text(
      size = 8
    ),

    plot.title = element_text(
      angle = 0,
      hjust = 0.5
    ),

    axis.title = element_blank(),

    legend.position = "right",

    legend.title = element_text(
      angle = 90,
      hjust = 0.5
    ),

    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  )