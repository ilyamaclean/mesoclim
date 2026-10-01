# =============================================================================
# ukcp_sample_data.R
#
# Creates package datasets ukcpinput (May 2018), ukcpfuture (May 2030) and
# bcmodel_list (bias-correction models fitted for May 2018) for Southwest UK on
# a 4 x 3 grid of 12 km UKCP18 RCM cells (x 144000-192000, y 0-36000, OSGB),
# which covers inst/extdata/dtms/lizard50m.tif and dtmf.tif.
#
# Run from the package root after devtools::load_all():
#   source("data-raw/ukcp_sample_data.R")
#
# Requires the full UKCP18 RCM daily files (land-rcm, uk, rcp85, member 01) for
# 20101201-20201130 and 20201201-20301130 in dir_ukcp.
# =============================================================================
library(terra)

dir_ukcp <- "/Users/jonathanmosedale/Data/mesoclim_inputs"
grid_ext <- ext(144000, 192000, 0, 36000)

ukdtmc <- rast(system.file("extdata/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc", package = "mesoclim"))
crs(ukdtmc) <- "EPSG:27700"
dtmc <- crop(ukdtmc, grid_ext)

# ukcpinput: May 2018, climate variables stored as arrays
ukcpinput <- ukcp18toclimarray(dir_ukcp, dtmc, as.POSIXlt("2018/05/01", tz = "UTC"), as.POSIXlt("2018/05/31", tz = "UTC"),
                               collection = "land-rcm", domain = "uk", member = "01")
ukcpinput$dtm <- wrap(ukcpinput$dtm)
usethis::use_data(ukcpinput, overwrite = TRUE)

# ukcpfuture: May 2030, climate variables stored as PackedSpatRasters
ukcpfuture <- ukcp18toclimarray(dir_ukcp, dtmc, as.POSIXlt("2030/05/01", tz = "UTC"), as.POSIXlt("2030/05/31", tz = "UTC"),
                                collection = "land-rcm", domain = "uk", member = "01", toArrays = FALSE)
ukcpfuture <- lapply(ukcpfuture, function(x) if (inherits(x, "SpatRaster")) wrap(x) else x)
usethis::use_data(ukcpfuture, overwrite = TRUE)

# bcmodel_list: tmax, tmin and prec models fitted to HadUK-Grid 1 km data
# (aggregated to 12 km) for May 2018, as in vignette 4
modeldata <- read_climdata(ukcpinput)
lsm12km <- crop(rast(system.file("extdata/biascorrect/uk_seamask_12km.tif", package = "mesoclim")), modeldata$dtm)
prep_obs <- function(f) {
  r <- rast(system.file(file.path("extdata/haduk", f), package = "mesoclim"))
  r <- crop(extend(r, ext(lsm12km)), ext(lsm12km))
  r <- aggregate(r, 12, fun = "mean", na.rm = TRUE)
  mask(.spatinterp(r), lsm12km, maskvalue = 0)
}
obs <- list(prec = prep_obs("rainfall1km.tif"), tmax = prep_obs("tasmax1km.tif"), tmin = prep_obs("tasmin1km.tif"))
msk <- ifel(lsm12km == 0, NA, lsm12km)
days <- nlyr(obs$tmax)
model_list <- list()
for (v in c("tmax", "tmin")) {
  model_list[[v]] <- biascorrect(mask(obs[[v]], msk), hist_mod = mask(modeldata[[v]], msk), fut_mod = NA,
                                 mod_out = TRUE, rangelims = NA, samplenum = days)
}
model_list$prec <- precipcorrect(mask(obs$prec, msk), hist_mod = mask(modeldata$prec, msk), mod_out = TRUE)
save(model_list, file = "data/bcmodel_list.rda", compress = "bzip2")
