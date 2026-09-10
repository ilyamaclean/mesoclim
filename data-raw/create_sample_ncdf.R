# =============================================================================
# create_sample_ncdf.R
#
# Creates small sample NetCDF files from large ERA5 and UKCP18RCM source files
# for use in package examples and vignettes. Output files are kept under 1 MB
# each and preserve the full variable/attribute structure of the originals.
#
# Run from the package root directory:
#   source("data-raw/create_sample_ncdf.R")
#
# Dependencies: ncdf4
# =============================================================================

library(ncdf4)

# =============================================================================
# USER CONFIGURATION — adjust paths and extents as needed
# =============================================================================

# Source files (relative to package root, or use absolute paths)
era5_file <- "raw_data_files/era5/era5_surface_ukeire_2020.nc"
ukcp_dir  <- "raw_data_files/ukcprcm"

# Output directories (created if absent)
era5_out  <- "inst/extdata/era5raw"
ukcp_out  <- "inst/extdata/ukcprcm_sample"

# ERA5: spatial bounding box (WGS84 lon/lat) and number of hourly time steps
# Default: SW Scotland, 1 month = 744 hours
# At 0.25° resolution: 7 lon x 5 lat x 744 h x 12 vars x 2 bytes ≈ 310 KB
era5_lon_range  <- c(-6, -4.5) # longitude (degrees East)  — 7 grid cells
era5_lat_range  <- c(55, 56)   # latitude  (degrees North) — 5 grid cells
era5_n_time     <- 744         # hourly steps (744 = 1 month)

# UKCP18: spatial bounding box in British National Grid (metres) and number of
# daily time steps.  Default: Scottish Highlands, 2 months (60 days of 360-day
# calendar).  At 12 km resolution 11 x 11 cells × 60 days × 11 files ≈ 600 KB
ukcp_x_range <- c(200000, 320000)  # BNG easting  (metres)
ukcp_y_range <- c(750000, 870000)  # BNG northing (metres)
ukcp_n_time  <- 61                 # daily steps = 2 months Dec to Jan

# ERA5 Southwest
era5_lon_range  <- c(-5.9, -4) # longitude (degrees East)  — 7 grid cells
era5_lat_range  <- c(49.8, 51.1)   # latitude  (degrees North) — 5 grid cells

# UKCP Southwest
ukcp_x_range <- c(130000, 250000)  # BNG easting  (metres)
ukcp_y_range <- c(10000, 130000)  # BNG northing (metres)

# SST: spatial bounding box (WGS84 lon/lat) and output directory.
# Files contain 12 monthly time steps per year — all kept by default (Inf).
# Grid resolution: 0.111° lon x 0.067° lat.
# Default: SW England / Celtic Sea, ~36 lon x 45 lat x 12 months x 10 vars ≈ 760 KB
sst_dir       <- "raw_data_files/sst"
sst_out       <- "inst/extdata/sst_sample"
sst_lon_range <- c(-7, -3)    # longitude range (degrees East)
sst_lat_range <- c(48, 51)    # latitude  range (degrees North) — 48 cells ≈ 760 KB
sst_n_time    <- Inf          # keep all monthly steps (12 per annual file)

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

# Copy all non-FillValue attributes from one nc object to another
.copy_atts <- function(nc_in, nc_out, varid) {
  atts <- ncatt_get(nc_in, varid)
  for (nm in names(atts)) {
    if (nm != "_FillValue")  # must be set at var creation; skip here
      ncatt_put(nc_out, varid, nm, atts[[nm]])
  }
}

# Map ncdf4 precision strings to ncvar_def prec argument
.nc_prec <- function(p) {
  switch(p,
    float  = "float",
    double = "double",
    int    = "integer",
    short  = "short",
    char   = "char",
    byte   = "byte",
    "float"   # safe fallback
  )
}

# Return a type-safe missing value for ncvar_def.
# ncdf4 cannot write an R NA as _FillValue for integer/short/byte types;
# use CF-convention defaults for those types instead.
.nc_missval <- function(prec, fill_att) {
  if (fill_att$hasatt) return(fill_att$value)
  switch(prec,
    float  = NA,
    double = NA,
    int    = 2147483647L,
    short  = 32767L,
    byte   = 127L,
    char   = NULL,
    NA
  )
}

# =============================================================================
# 1.  SUBSET AN ERA5 FILE
# =============================================================================
# ERA5 structure (this file): all variables share dims [longitude x latitude x time],
# all stored as 'short' with scale/offset attributes applied automatically by ncdf4.

subset_era5 <- function(infile, outfile, lon_range, lat_range, n_time) {

  cat("\n=== ERA5 ===\n")
  nc <- nc_open(infile)
  on.exit(nc_close(nc), add = TRUE)

  lons  <- nc$dim$longitude$vals
  lats  <- nc$dim$latitude$vals
  times <- nc$dim$time$vals

  lon_idx  <- which(lons  >= lon_range[1] & lons  <= lon_range[2])
  lat_idx  <- which(lats  >= lat_range[1] & lats  <= lat_range[2])
  time_idx <- seq_len(min(n_time, length(times)))

  if (length(lon_idx) == 0) stop("No ERA5 longitudes in requested range")
  if (length(lat_idx) == 0) stop("No ERA5 latitudes in requested range")

  cat("  Lons :", length(lon_idx), "cells [", lons[range(lon_idx)], "]\n")
  cat("  Lats :", length(lat_idx), "cells [", lats[range(lat_idx)], "]\n")
  cat("  Times:", length(time_idx), "steps\n")

  # --- New dimension objects ---
  dim_lon  <- ncdim_def("longitude", nc$dim$longitude$units, lons[lon_idx])
  dim_lat  <- ncdim_def("latitude",  nc$dim$latitude$units,  lats[lat_idx])
  dim_time <- ncdim_def("time",      nc$dim$time$units,      times[time_idx],
                         unlim = TRUE)

  dim_map <- list(longitude = dim_lon, latitude = dim_lat, time = dim_time)

  # --- Variable definitions ---
  var_defs <- lapply(names(nc$var), function(vname) {
    v     <- nc$var[[vname]]
    vdims <- lapply(sapply(v$dim, `[[`, "name"), function(dn) dim_map[[dn]])
    fv    <- ncatt_get(nc, vname, "_FillValue")
    ncvar_def(
      name    = vname,
      units   = if (nchar(v$units) > 0) v$units else " ",
      dim     = vdims,
      missval = .nc_missval(v$prec, fv),
      prec    = .nc_prec(v$prec)
    )
  })
  names(var_defs) <- names(nc$var)

  # --- Create file and copy global attributes ---
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  nc_out <- nc_create(outfile, vars = var_defs)
  on.exit(nc_close(nc_out), add = TRUE)

  for (nm in names(ncatt_get(nc, 0)))
    ncatt_put(nc_out, 0, nm, ncatt_get(nc, 0)[[nm]])

  for (dn in c("longitude", "latitude", "time"))
    .copy_atts(nc, nc_out, dn)

  # --- Write subsetted data and copy variable attributes ---
  for (vname in names(nc$var)) {
    v      <- nc$var[[vname]]
    dnames <- sapply(v$dim, `[[`, "name")
    start  <- rep(1L, length(dnames))
    count  <- sapply(v$dim, `[[`, "len")

    for (i in seq_along(dnames)) {
      if (dnames[i] == "longitude") { start[i] <- min(lon_idx);  count[i] <- length(lon_idx) }
      if (dnames[i] == "latitude")  { start[i] <- min(lat_idx);  count[i] <- length(lat_idx) }
      if (dnames[i] == "time")      { start[i] <- 1L;            count[i] <- length(time_idx) }
    }

    dat <- ncvar_get(nc, vname, start = start, count = count, raw_datavals = TRUE)
    ncvar_put(nc_out, vname, dat)
    .copy_atts(nc, nc_out, vname)
  }

  cat("  Written:", outfile, sprintf("(%.0f KB)\n", file.size(outfile) / 1024))
}

# =============================================================================
# 2.  SUBSET A SINGLE UKCP18 RCM FILE
# =============================================================================
# UKCP18 RCM structure (12 km, daily, BNG):
#   Main var   : [projection_x_coordinate x projection_y_coordinate x time x ensemble_member]
#   Aux 2-D    : grid_latitude, grid_longitude [x x y]
#   Bounds     : time_bnds [bnds x time], *_bnds [bnds x x/y]
#   Char vars  : ensemble_member_id [string27 x ensemble_member]
#                yyyymmdd           [string64 x time]
#   Scalars    : transverse_mercator (no dims), month_number [time], year [time]

subset_ukcp18 <- function(infile, outfile, x_range, y_range, n_time) {

  cat("\n=== UKCP18:", basename(infile), "===\n")
  nc <- nc_open(infile)
  on.exit(nc_close(nc), add = TRUE)

  xvals <- nc$dim$projection_x_coordinate$vals
  yvals <- nc$dim$projection_y_coordinate$vals
  tvals <- nc$dim$time$vals

  x_idx   <- which(xvals >= x_range[1] & xvals <= x_range[2])
  y_idx   <- which(yvals >= y_range[1] & yvals <= y_range[2])
  t_idx   <- seq_len(min(n_time, length(tvals)))

  if (length(x_idx) == 0) stop("No UKCP18 x-coords in requested range")
  if (length(y_idx) == 0) stop("No UKCP18 y-coords in requested range")

  cat("  X :", length(x_idx), "cells [", xvals[range(x_idx)], "] m\n")
  cat("  Y :", length(y_idx), "cells [", yvals[range(y_idx)], "] m\n")
  cat("  T :", length(t_idx), "steps\n")

  # --- Dimension objects (only what exists in this file) ---
  mk_dim <- function(dname) {
    d <- nc$dim[[dname]]
    if (is.null(d)) return(NULL)
    switch(dname,
      projection_x_coordinate = ncdim_def(dname, d$units, xvals[x_idx]),
      projection_y_coordinate = ncdim_def(dname, d$units, yvals[y_idx]),
      time                    = ncdim_def(dname, d$units, tvals[t_idx], unlim = TRUE),
      ensemble_member         = ncdim_def(dname, d$units, d$vals),
      bnds                    = ncdim_def(dname, "",       1:2),
      string27                = ncdim_def(dname, "",       1:27),
      string64                = ncdim_def(dname, "",       1:64),
      ncdim_def(dname, d$units, d$vals)   # passthrough for any others
    )
  }

  dim_map <- Filter(Negate(is.null),
                    setNames(lapply(names(nc$dim), mk_dim), names(nc$dim)))

  # --- Variable definitions ---
  # transverse_mercator is a scalar (0-dim) variable used only for CF grid_mapping
  var_defs <- list()
  for (vname in names(nc$var)) {
    v      <- nc$var[[vname]]
    dnames <- sapply(v$dim, `[[`, "name")

    if (length(dnames) == 0) {
      # Scalar variable (transverse_mercator)
      vdims <- list()
    } else {
      vdims <- lapply(dnames, function(dn) dim_map[[dn]])
      if (any(sapply(vdims, is.null))) {
        cat("  Skipping", vname, "(unrecognised dimension)\n")
        next
      }
    }

    fv <- ncatt_get(nc, vname, "_FillValue")
    var_defs[[vname]] <- ncvar_def(
      name    = vname,
      units   = if (!is.null(v$units) && v$units != "0") v$units else " ",
      dim     = vdims,
      missval = .nc_missval(v$prec, fv),
      prec    = .nc_prec(v$prec)
    )
  }

  # --- Create file and copy global attributes ---
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  # force_v4 = TRUE required because UKCP18 uses an unlimited time dimension
  # that is not in the last position — a NetCDF-4 (HDF5) feature.
  nc_out <- nc_create(outfile, vars = var_defs, force_v4 = TRUE)
  on.exit(nc_close(nc_out), add = TRUE)

  for (nm in names(ncatt_get(nc, 0)))
    ncatt_put(nc_out, 0, nm, ncatt_get(nc, 0)[[nm]])

  for (dn in names(dim_map))
    if (dn %in% names(nc$dim)) .copy_atts(nc, nc_out, dn)

  # --- Write subsetted data and copy variable attributes ---
  for (vname in names(var_defs)) {
    v      <- nc$var[[vname]]
    dnames <- sapply(v$dim, `[[`, "name")

    if (length(dnames) == 0) {
      # Scalar: read and write without start/count
      dat <- ncvar_get(nc, vname)
      ncvar_put(nc_out, vname, dat)
    } else {
      start <- rep(1L, length(dnames))
      count <- sapply(v$dim, `[[`, "len")

      for (i in seq_along(dnames)) {
        if (dnames[i] == "projection_x_coordinate") { start[i] <- min(x_idx); count[i] <- length(x_idx) }
        if (dnames[i] == "projection_y_coordinate") { start[i] <- min(y_idx); count[i] <- length(y_idx) }
        if (dnames[i] == "time")                    { start[i] <- 1L;         count[i] <- length(t_idx) }
        # bnds, ensemble_member, string* dims: full extent kept (count unchanged)
      }

      tryCatch({
        dat <- ncvar_get(nc, vname, start = start, count = count)
        ncvar_put(nc_out, vname, dat)
      }, error = function(e) {
        cat("  Warning: could not write", vname, "—", conditionMessage(e), "\n")
      })
    }

    .copy_atts(nc, nc_out, vname)
  }

  cat("  Written:", outfile, sprintf("(%.0f KB)\n", file.size(outfile) / 1024))
}

# =============================================================================
# 3.  SUBSET AN SST FILE
# =============================================================================
# SST structure (NWSClim NWS PPE, monthly):
#   dims : lon (297, 0.111°), lat (375, 0.067°), time (12 monthly), bnds (2)
#   vars : SST, NBT, DFT, SSS, NBS, DFS, SSH, MLD, PEA, DMUV  [lon x lat x time]
#          time_bounds                                           [bnds x time]
# All main vars stored as float with 1e+20 fill value.
# time units: seconds since 1950-01-01.  No special calendar.

subset_sst <- function(infile, outfile, lon_range, lat_range, n_time = Inf) {

  cat("\n=== SST:", basename(infile), "===\n")
  nc <- nc_open(infile)
  on.exit(nc_close(nc), add = TRUE)

  lons  <- nc$dim$lon$vals
  lats  <- nc$dim$lat$vals
  times <- nc$dim$time$vals

  lon_idx  <- which(lons  >= lon_range[1] & lons  <= lon_range[2])
  lat_idx  <- which(lats  >= lat_range[1] & lats  <= lat_range[2])
  time_idx <- seq_len(min(n_time, length(times)))

  if (length(lon_idx) == 0) stop("No SST longitudes in requested range")
  if (length(lat_idx) == 0) stop("No SST latitudes in requested range")

  cat("  Lon :", length(lon_idx), "cells [", lons[range(lon_idx)], "]\n")
  cat("  Lat :", length(lat_idx), "cells [", lats[range(lat_idx)], "]\n")
  cat("  Time:", length(time_idx), "monthly steps\n")

  # --- New dimension objects ---
  dim_lon  <- ncdim_def("lon",  nc$dim$lon$units,  lons[lon_idx])
  dim_lat  <- ncdim_def("lat",  nc$dim$lat$units,  lats[lat_idx])
  dim_time <- ncdim_def("time", nc$dim$time$units, times[time_idx], unlim = TRUE)
  dim_bnds <- ncdim_def("bnds", "", 1:2)

  dim_map  <- list(lon = dim_lon, lat = dim_lat, time = dim_time, bnds = dim_bnds)

  # --- Variable definitions ---
  var_defs <- list()
  for (vname in names(nc$var)) {
    v      <- nc$var[[vname]]
    dnames <- sapply(v$dim, `[[`, "name")
    vdims  <- lapply(dnames, function(dn) dim_map[[dn]])

    if (any(sapply(vdims, is.null))) {
      cat("  Skipping", vname, "(unrecognised dimension)\n")
      next
    }

    fv <- ncatt_get(nc, vname, "_FillValue")
    var_defs[[vname]] <- ncvar_def(
      name    = vname,
      units   = if (nchar(v$units) > 0) v$units else " ",
      dim     = vdims,
      missval = .nc_missval(v$prec, fv),
      prec    = .nc_prec(v$prec)
    )
  }

  # --- Create file and copy global attributes ---
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  nc_out <- nc_create(outfile, vars = var_defs)
  on.exit(nc_close(nc_out), add = TRUE)

  for (nm in names(ncatt_get(nc, 0)))
    ncatt_put(nc_out, 0, nm, ncatt_get(nc, 0)[[nm]])

  for (dn in c("lon", "lat", "time"))
    .copy_atts(nc, nc_out, dn)

  # --- Write subsetted data and copy variable attributes ---
  for (vname in names(var_defs)) {
    v      <- nc$var[[vname]]
    dnames <- sapply(v$dim, `[[`, "name")
    start  <- rep(1L, length(dnames))
    count  <- sapply(v$dim, `[[`, "len")

    for (i in seq_along(dnames)) {
      if (dnames[i] == "lon")  { start[i] <- min(lon_idx);  count[i] <- length(lon_idx) }
      if (dnames[i] == "lat")  { start[i] <- min(lat_idx);  count[i] <- length(lat_idx) }
      if (dnames[i] == "time") { start[i] <- 1L;            count[i] <- length(time_idx) }
      # bnds: full extent kept
    }

    dat <- ncvar_get(nc, vname, start = start, count = count)
    ncvar_put(nc_out, vname, dat)
    .copy_atts(nc, nc_out, vname)
  }

  cat("  Written:", outfile, sprintf("(%.0f KB)\n", file.size(outfile) / 1024))
}

# =============================================================================
# RUN
# =============================================================================

cat("Creating sample NetCDF files...\n")
cat("ERA5 source    :", era5_file, "\n")
cat("UKCP18 source  :", ukcp_dir, "\n")

# ERA5
subset_era5(
  infile    = era5_file,
  outfile   = file.path(era5_out, "era5_sample.nc"),
  lon_range = era5_lon_range,
  lat_range = era5_lat_range,
  n_time    = era5_n_time
)

# UKCP18 — process every .nc file in the directory
ukcp_files <- list.files(ukcp_dir, pattern = "\\.nc$", full.names = TRUE)
if (length(ukcp_files) == 0) warning("No .nc files found in: ", ukcp_dir)

for (f in ukcp_files) {
  # Shorten the decade suffix in the filename (e.g. _20201201-20301130 → _sample)
  outname <- sub("_[0-9]{8}-[0-9]{8}(\\.nc)$", "_sample\\1", basename(f))
  # OR just add sample but include decade
  #outname <- gsub("(_[0-9]{8}-[0-9]{8})", "\\1_sample", basename(f))
  outname <- basename(f)
  subset_ukcp18(
    infile   = f,
    outfile  = file.path(ukcp_out, outname),
    x_range  = ukcp_x_range,
    y_range  = ukcp_y_range,
    n_time   = ukcp_n_time
  )
}

# SST — process every annual gridT file in the directory
sst_files <- list.files(sst_dir, pattern = "_gridT\\.nc$", full.names = TRUE)
if (length(sst_files) == 0) warning("No gridT .nc files found in: ", sst_dir)

for (f in sst_files) {
  subset_sst(
    infile    = f,
    outfile   = file.path(sst_out, basename(f)),
    lon_range = sst_lon_range,
    lat_range = sst_lat_range,
    n_time    = sst_n_time
  )
}

cat("\nDone.\n")
cat("ERA5 output  :", era5_out, "\n")
cat("UKCP18 output:", ukcp_out, "\n")
cat("SST output   :", sst_out, "\n")
