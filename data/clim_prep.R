library(raster)
library(hydro)
library(tidyverse)

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

# temperature: load
t1 <- f[grepl("_tas_", f)] %>%
      rast() %>% crop(bb) %>% mean() %>%
      setNames("tmean")
clim1 <- c(t1, setNames(rast(hst), layers))


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

# temperature: compute from monthlies
t2 <- f[grepl("_tas_", f)] %>%
      rast() %>% crop(bb) %>% mean() %>%
      setNames("tmean")
clim2 <- c(t2, setNames(rast(fut), layers))



# munge -----------------

# select focal variables
clim1 <- clim1[[c("tmean", "ppt", "aet", "cwd")]]
clim2 <- clim2[[c("tmean", "ppt", "aet", "cwd")]]

# log-transform precip
clim1$ppt <- log10(clim1$ppt)
clim2$ppt <- log10(clim2$ppt)

# compute, apply, and store scales
scales <- list()
for(var in names(clim1)){
      x <- values(clim1[[var]])

      m <- mean(x, na.rm = T)
      s <- sd(x, na.rm = T)

      clim1[[var]] <- (clim1[[var]] - m) / s
      clim2[[var]] <- (clim2[[var]] - m) / s

      scales[[var]]$mean <- m
      scales[[var]]$sd <- s
}
saveRDS(scales, "data/clim_scales.rds")

# export
clim <- c(clim1, clim2) %>%
      setNames(c(paste0(names(clim1), 1),
                 paste0(names(clim2), 2)))

writeRaster(
      clim,
      "data/clim_quantized.tif",
      datatype = "INT2S",
      scale = 0.01,
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
      overwrite = TRUE
)
