landsea_template <- function(xmin, xmax, ymin, ymax, res = 100, land = NULL) {
  r <- terra::rast(terra::ext(xmin, xmax, ymin, ymax), res = res, crs = "EPSG:27700")
  terra::values(r) <- NA
  if (!is.null(land)) r[terra::cells(r, land)] <- 1
  r
}

test_that("coastalexposure returns land fraction with sea as NA", {
  r <- landsea_template(0, 20000, 0, 20000, land = terra::ext(8000, 12000, 8000, 12000))
  ce <- coastalexposure(r, 270)
  v <- terra::values(ce)
  expect_true(all(is.na(v[is.na(terra::values(r))])))
  expect_true(all(v >= 0 & v <= 1, na.rm = TRUE))
  expect_lt(max(v, na.rm = TRUE), 1)
})

test_that("coastalexposure increases downwind of the coast", {
  r <- landsea_template(0, 20000, 0, 10000, land = terra::ext(10000, 20000, 0, 10000))
  ce <- coastalexposure(r, 270)
  near <- terra::extract(ce, cbind(10550, 5050))[1, 1]
  far <- terra::extract(ce, cbind(19450, 5050))[1, 1]
  expect_lt(near, far)
  expect_equal(unique(stats::na.omit(terra::values(coastalexposure(r, 90))[, 1])), 1)
})

test_that("coastalexposure uses coarse rasters beyond landsea and is reproducible with jitter", {
  fine <- landsea_template(10000, 12000, 0, 2000, land = terra::ext(10000, 12000, 0, 2000))
  coarse <- landsea_template(0, 20000, -10000, 12000, res = 1000, land = terra::ext(9000, 20000, -10000, 12000))
  expect_equal(unique(terra::values(coastalexposure(fine, 270))[, 1]), 1)
  ce <- coastalexposure(fine, 270, coarse = coarse)
  expect_lt(max(terra::values(ce), na.rm = TRUE), 1)
  expect_identical(terra::values(ce), terra::values(coastalexposure(fine, 270, coarse = coarse)))
})

test_that("coastalexposure requires a projected coordinate reference system", {
  r <- terra::rast(terra::ext(-5, -4, 50, 51), res = 0.1, crs = "EPSG:4326")
  terra::values(r) <- 1
  expect_snapshot(coastalexposure(r, 270), error = TRUE)
})

test_that("calculate_coastalexposure returns directional and mean layers matching dtmf", {
  dtmf <- terra::rast(system.file("extdata/dtms/dtmf.tif", package = "mesoclim"))
  dtmm <- terra::rast(system.file("extdata/dtms/dtmm.tif", package = "mesoclim"))
  cex <- calculate_coastalexposure(dtmf, dtmm, ndir = 8)
  expect_equal(terra::nlyr(cex), 9)
  expect_equal(names(cex), c(paste0("dir_", seq(0, 315, 45)), "all"))
  expect_true(terra::compareGeom(cex, dtmf))
  expect_true(all(terra::values(cex) >= 0 & terra::values(cex) <= 1, na.rm = TRUE))
})

test_that(".check_cex crops larger exposure to dtmf and rejects mismatched grids", {
  r <- landsea_template(0, 4000, 0, 4000, land = terra::ext(1000, 4000, 0, 4000))
  cex <- calculate_coastalexposure(r, ndir = 4, smooth = 3)
  sub <- terra::crop(r, terra::ext(2000, 4000, 0, 2000))
  expect_true(terra::compareGeom(.check_cex(cex, sub), sub))
  expect_snapshot(.check_cex(terra::aggregate(cex, 2), sub), error = TRUE)
})

test_that(".tempcoastal works for a single timestep and moves temperature towards SST", {
  r <- landsea_template(0, 4000, 0, 4000, land = terra::ext(1000, 4000, 0, 4000))
  cex <- calculate_coastalexposure(r, ndir = 8, smooth = 3)
  dtmc <- terra::rast(terra::ext(0, 4000, 0, 4000), res = 2000, crs = "EPSG:27700")
  terra::values(dtmc) <- 50
  wdir <- terra::rast(dtmc)
  terra::values(wdir) <- 270
  tc <- r * 0 + 5
  out <- .tempcoastal(tc, sstf = r * 0 + 12, u2 = r * 0 + 5, wdir = wdir, dtmc = dtmc, cex = cex)
  v <- terra::values(out - tc)
  expect_equal(terra::nlyr(out), 1)
  expect_true(all(v >= 0, na.rm = TRUE))
  expect_gt(max(v, na.rm = TRUE), 0)
})

test_that("tempdaily_downscale gives the same result with supplied or internally calculated exposure", {
  climdata <- read_climdata(mesoclim::ukcpinput)
  dtmf <- terra::rast(system.file("extdata/dtms/dtmf.tif", package = "mesoclim"))
  dtmm <- terra::rast(system.file("extdata/dtms/dtmm.tif", package = "mesoclim"))
  sst <- terra::unwrap(mesoclim::ukcp18sst)
  climdata <- subset_climdata(climdata, as.POSIXlt("2018-05-01", tz = "UTC"), as.POSIXlt("2018-05-03", tz = "UTC"))
  wca <- calculate_windcoeffs(climdata$dtm, dtmm, dtmf, zo = 2)
  uzf <- winddownscale(climdata$windspeed, climdata$winddir, dtmf, dtmm, climdata$dtm, wca, zi = climdata$windheight_m, zo = 2)
  cex <- calculate_coastalexposure(dtmf, dtmm)
  a <- tempdaily_downscale(climdata, NA, sst, dtmf, dtmm, NA, uzf, cad = FALSE, coastal = TRUE)
  b <- tempdaily_downscale(climdata, NA, sst, dtmf, dtmm, NA, uzf, cad = FALSE, coastal = TRUE, cex = cex)
  expect_equal(terra::values(a$tmax), terra::values(b$tmax))
  expect_equal(terra::values(a$tmin), terra::values(b$tmin))
})
