# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`mesoclim` is an R package for mechanistic downscaling of climate data from coarse (10s km) to fine (<1 km) spatial resolution and from daily to hourly temporal resolution. It handles ERA5 reanalysis data and UKCP18 climate model projections.

## Standard workflow and use
1 Prepare coarse resolution climate data (eg ERA5 or UKCP18RCM) to standard format (climdata)
2 Prepare other data required for downscaling including digital elevation maps (dtms) at different resolutions and extents
3 Bias correction against observational data at coarse resolution
4 Temporal and/or Spatial Downscaling 

The Analyses has many examples of use of the package

## Development Commands

```r
# Install dependencies and build package
devtools::install_deps()
devtools::build()

# Load package for interactive development (faster than install)
devtools::load_all()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-workerfunctions.R")

# Run a specific test by name
testthat::test_file("tests/testthat/test-workerfunctions.R", filter = "test name pattern")

# Generate documentation from roxygen2 comments
devtools::document()

# Check package for CRAN compliance
devtools::check()

# Compile C++ code (Rcpp)
Rcpp::compileAttributes()
devtools::build(compile_attributes = TRUE)
```

After modifying C++ files in `src/`, run `Rcpp::compileAttributes()` to regenerate `R/RcppExports.R` and `src/RcppExports.cpp` before testing.

## Architecture

### Standard "climdata" Format

All climate data is standardized into a named list with the following structure:
- Named numeric arrays with dimensions `[x, y, time]` for variables: `tmin`, `tmax`, `temp`, `relhum`, `pres`, `swrad`, `lwrad`, `windspeed`, `winddir`, `prec`
- `dtm`: a `SpatRaster` of the elevation grid
- `dtmc`: optional coarse-resolution DTM
- `times`: vector of POSIXct timestamps

Functions in `data_functions.R` handle I/O (`read_climdata`, `write_climdata`) and validation (`checkinputs`).

### Processing Pipeline

1. **Input conversion** (`era5_functions.R`, `ukcp_functions.R`): NetCDF → climdata list
2. **Spatial downscaling** (`spatialdownscale.R`, `workerfunctions.R`): Applies elevation, coastal, and cold-air drainage effects
3. **Temporal downscaling** (`temporaldownscale.R`): Daily → hourly via physically-based diurnal curves
4. **Bias correction** (`biascorrect_functions.R`): Statistical correction against observations
5. **Large-domain tiling** (`tile_functions.R`): Overlapping tiles to handle memory constraints

### File Responsibilities

| File | Purpose |
|------|---------|
| `spatialdownscale.R` | Public spatial downscaling functions (`spatialdownscale`, `winddownscale`, `precipdownscale`, `tempdownscale`) |
| `temporaldownscale.R` | Daily→hourly conversion for all variables (`*_dailytohourly`, `subdailyrain`, `temporaldownscale`) |
| `climatefunctions.R` | Physical process models: radiation, humidity, lapse rate, coastal exposure, wind elevation |
| `workerfunctions.R` | ~50 internal helpers (dot-prefixed): type conversions, spatial interpolation, terrain analysis, cold-air drainage (`basindelin`, `flowacc`, `.tempcad`, `.tempcoastal`, `.tempelev`) |
| `biascorrect_functions.R` | Bias correction model fitting and application |
| `era5_functions.R` | ERA5 NetCDF → climdata conversion |
| `ukcp_functions.R` | UKCP18 NetCDF → climdata conversion |
| `data_functions.R` | I/O, validation, plotting |
| `tile_functions.R` | Tiled processing for large domains |
| `data.R` | Dataset documentation only |

### C++ Acceleration (src/)

Performance-critical calculations in `mesoclimCpp.cpp` are called via Rcpp wrappers in `RcppExports.R`:
- `hourlytempv` / `hourlytempm`: Daily→hourly temperature interpolation (vector/matrix variants)
- `clearskyradmCpp`: Clear-sky solar radiation for a matrix of lat/lon/times
- `difpropmCpp`: Diffuse/direct radiation partitioning
- `rainadjustv` / `rainadjustm`: Precipitation temporal disaggregation
- `basinCpp`: Basin delineation (flow-direction kernel)
- `hourtodayCpp`: Hourly→daily aggregation

### Internal Conventions

- Internal functions are dot-prefixed (e.g., `.resample()`, `.spatinterp()`, `.tempelev()`)
- Data type conversion helpers: `.is()` (SpatRaster→array), `.rta()` (raster→array), `.vta()` (vector→array), `.rast()` (array→SpatRaster)
- Functions accept both `terra::SpatRaster` and legacy `raster` objects; internally converted as needed
- `PackedSpatRaster` is used for serialization (see `write_climdata`/`read_climdata`)

### Test Data

Example datasets are in `inst/extdata/`:
- `dtms/`: Fine (dtmf) and medium (dtmm) resolution DTMs for test area
- `era5/`: Sample ERA5 NetCDF files (instant, accum, avg stream types)
- `haduk/`: Met Office HadUK gridded observations
- `ukcp18rcm/`: UKCP18 regional climate model orography
- `sst/`: Sea surface temperature from NWS PPE

Package-level datasets in `data/` (loaded via `data(name)`): `climdata`, `era5input`, `daily100m`, `era5sst`, `ukcp18sst`, `landsea`, `ukcphighland`, `ukcpinput`, etc.
