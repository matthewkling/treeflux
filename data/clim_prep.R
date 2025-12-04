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
writeRaster(clim, "data/clim.tif", overwrite = TRUE)



# # bb <- ext(-125, -66, 24, 50)
#
# scale_rast <- function(x) (x - mean(values(x), na.rm = T)) / sd(values(x), na.rm = T)
#
#
# t1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio1_1981-2010_V.2.1.tif") %>% crop(bb)
# t2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
# t2 <- t2[grepl("2071", t2) & grepl("ssp370_tas_", t2)] %>% rast() %>% crop(bb) %>% mean()
# t <- c(t1, t2) %>% scale_rast()
#
# p1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio12_1981-2010_V.2.1.tif") %>% crop(bb)
# p2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
# p2 <- p2[grepl("2071", p2) & grepl("ssp370_pr_", p2)] %>% rast() %>% crop(bb) %>% mean() %>% "*"(12)
# p <- c(p1, p2) %>% log() %>% scale_rast()
#
# clim1 <- c(t[[1]], p[[1]])
# clim2 <- c(t[[2]], p[[2]])
#
# clim <- c(clim1, clim2) %>% setNames(c("t1", "p1", "t2", "p2"))
#
# writeRaster(clim, "data/clim.tif", overwrite = TRUE)
