library(raster)
library(hydro)
library(tidyverse)
library(terra)

bb <- extent(-125, -66, 24, 50)


# historic -------------

# water balance: compute from monthlies
f <- list.files("/Volumes/T7/CHELSA/v2/raw/", full.names = T)
f <- f[grepl("tas|pr", f)]
hst <- f %>%
      stack() %>%
      crop(bb) %>%
      hydro(already_latlong = TRUE, ncores = 8)
layers <- tolower(names(hst))

# mean annual temperature
tmean1 <- f[grepl("_tas_", f)] %>%
      rast() %>% crop(bb) %>% mean() %>%
      setNames("tmean")

# tmax of warmest month (bio5)
tmaxwm1 <- f[grepl("_tasmax_", f)] %>%
      rast() %>% crop(bb) %>% max() %>%
      setNames("tmaxwm")

# tmin of coldest month (bio6)
tmincm1 <- f[grepl("_tasmin_", f)] %>%
      rast() %>% crop(bb) %>% min() %>%
      setNames("tmincm")

clim1 <- c(tmean1, tmaxwm1, tmincm1,
           setNames(rast(hst), layers))


# future ------------

# water balance: compute from monthlies
f <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
f <- f[grepl("tas|pr", f) &
             grepl("2071", f) &
             grepl("ssp585", f) &
             grepl("gfdl-esm4", f)]
fut <- f %>%
      stack() %>%
      crop(bb) %>%
      hydro(already_latlong = TRUE, ncores = 8)

# mean annual temperature
tmean2 <- f[grepl("_tas_", f)] %>%
      rast() %>% crop(bb) %>% mean() %>%
      setNames("tmean")

# tmax of warmest month (bio5)
tmaxwm2 <- f[grepl("_tasmax_", f)] %>%
      rast() %>% crop(bb) %>% max() %>%
      setNames("tmaxwm")

# tmin of coldest month (bio6)
tmincm2 <- f[grepl("_tasmin_", f)] %>%
      rast() %>% crop(bb) %>% min() %>%
      setNames("tmincm")

clim2 <- c(tmean2, tmaxwm2, tmincm2,
           setNames(rast(fut), layers))



# munge -----------------

# select focal variables
vars <- c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")
clim1 <- clim1[[vars]]
clim2 <- clim2[[vars]]

# log-transform precip
clim1$ppt <- log10(clim1$ppt)
clim2$ppt <- log10(clim2$ppt)

# compute, apply, and store scales
scales <- list()
for(var in vars){
      x <- values(clim1[[var]])

      m <- mean(x, na.rm = T)
      s <- sd(x, na.rm = T)

      clim1[[var]] <- (clim1[[var]] - m) / s
      clim2[[var]] <- (clim2[[var]] - m) / s

      scales[[var]]$mean <- m
      scales[[var]]$sd <- s
}
saveRDS(scales, "data/clim_scales.rds")


clim <- c(clim1, clim2) %>%
      setNames(c(paste0(vars, 1),
                 paste0(vars, 2)))

# mask land
land <- rast("/Volumes/T7/CHELSA/v2/raw//CHELSA_ai_1981-2010_V.2.1.tif") %>%
      crop(clim)
clim <- mask(clim, land)

# export
writeRaster(
      clim,
      "data/clim_quantized.tif",
      datatype = "INT2S",
      scale = 0.01,
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
      overwrite = TRUE
)
