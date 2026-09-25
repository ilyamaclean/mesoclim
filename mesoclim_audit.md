# mesoclim audit

Last updated 25 September 2026 (`main` @ `4791e6a`). This file lists open actions only. What has been done is recorded in the revision history at the end.

---

## 1. Setup and structure

**Package:** `mesoclim` 0.1.1 (GPL-2). It downscales coarse climate data (UKCP18 RCM 12 km, ERA5 0.25°) to sub-km resolution and from daily to hourly, with bias correction against HadUK-Grid. The target is GitHub plus a pkgdown site, not CRAN.

**Repository:** `https://github.com/ilyamaclean/mesoclim`.

| Branch | Role |
|---|---|
| `main` | Active development |
| `gh-pages` | pkgdown site output, written by `pkgdown::deploy_to_branch()` |
| `old_ver` | Old `main` before 10 Sept 2026 |
| `dev` | Superseded; everything is merged into `main` |

**Local only (gitignored):**
- `CLAUDE.md`, which holds the working guidance for Claude Code sessions.
- `docs/`, the pkgdown build output, now deleted locally.
- `analyses/`, `R_new/`, `raw_data_files/` and `extradocs/`. Don't save new files under `analyses/`.

This audit is kept in the repository root and committed so it can be shared. It's excluded from the package build via `.Rbuildignore`.

**Environment:**
- R 4.3.3 (arm64). Some CRAN packages now build from source; Rust is installed via Homebrew for `yaml12`, which `btw` needs.
- pandoc 3.10 comes from Positron, which `rmarkdown` picks up. The system pandoc at `/usr/local/bin` is 1.13.1 and too old for pkgdown.
- The `r-btw` MCP server gives Claude Code access to R documentation.

**Code layout:**
- `R/` holds 11 source files; see `CLAUDE.md` for what each file does.
- `src/mesoclimCpp.cpp` holds the Rcpp routines.
- `data/` holds 7 datasets, stored as wrapped SpatRasters that `read_climdata()` unwraps.
- `inst/extdata/` holds the sample NetCDF files, DTMs, HadUK observations and data-preparation scripts.
- `vignettes/` holds 4 vignettes, which are also the site articles. All run on package data at build time, except the downloads chunk in vignette 1 and the tiles example in vignette 2.

**Tests:** `devtools::test()` gives 62 passes and 1 skip. Every test is in `test-workerfunctions.R`.

---

## 2. Actions

### P1 — correctness and methods

1. **Precipitation temporal downscaling.** `temporaldownscale()` → `prec_dailytohourly()` spreads each day's total evenly over 24 hours (constant rate, with days below `noraincut` set to 0). `subdailyrain()` (Bartlett-Lewis / HyetosMinute approach) isn't wired in. Add it as a method, or decide how it should be used. `HyetosMinute` isn't on CRAN or declared, and must be attached with `library()`. Vignette 3's precipitation section is still marked "UNDER DEVELOPMENT".
2. **Tiling and coastal effect.**
   - Tile-boundary lines in `spatialdownscale_tiles()` output are probably caused by the coastal effect. Change the coastal effect method.
   - Get vignette 2's tiles example working at an acceptable run time; it currently takes 2.5 minutes and has `eval=FALSE`.
   - Fix the floating-point noise in `inst/extdata/dtms/lizard50m.tif`'s extent (xmin `160000.000192`), which triggers the "tile size NOT divisible by resolution" warning.
3. **Longwave radiation.** Check `lw_dailytohourly()` against the open questions in vignette 3: sky emissivity formulation, and whether to derive downward longwave from net longwave by removing an upward longwave term based on hourly temperature.
4. **Tmean.** Decide when daily mean temperature can be `(tmax+tmin)/2` and when it must come from hourly temperatures. Vignette 2's humidity and longwave examples now use `(tmax+tmin)/2`, while the internal code uses `.hourtoday()`.
5. **Issue #6:** `.resample()` (`R/workerfunctions.R`) calls `terra::resample()` without the BIGTIFF option, so large areas fail with "TIFFAppendToStrip: Maximum TIFF file size exceeded". Add `gdal = c("BIGTIFF=TRUE")`, or write to a temporary file with BIGTIFF.
6. **Issue #7 / PR #8:** `subset_climdata()` returns an unbound `newdata` when the input is neither a SpatRaster nor a list. Review PR #8. The SpatRaster branch also pads the selection by ±1 month, which the list branch doesn't do; check that this is intended.
7. **Tests.** 7 of the 8 test files are empty stubs. Priorities are regression tests for #6 and #7, `era5toclimarray()` (with and without `lsm`/`aoi`), `biascorrect_climdata()`, and a small `spatialdownscale()` / `temporaldownscale()` run on `ukcpinput`.
8. **Wind and precipitation fixes from `dev`.** These are merged, but `analyses/new_windz.R` and `analyses/new_precipdownscale.R` haven't been reviewed for further changes, and there are no tests for either.

### P1 — data and vignettes

9. **Data documentation and references.** Check that every dataset in `R/data.R` and every `inst/extdata` source has a correct description and citation (for example, `ukcp18sst` → Tinker et al. 2024).
10. **Geopotential → elevation.** Add a function to replace the inline conversion in vignette 1's ERA5 chunk. Check whether any other ancillary sources need it.
11. **Vignette 1:**
    - Where the UKCP `dtmc` comes from: it's cropped from `extdata/ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc` (UKCP orography). Explain this in the text.
    - `dtmm`'s southern edge (y = 10,000) lies 2 km outside the snapped `dtmc` (y = 12,000).
    - Add an NCEP example.
    - Add humidity documentation.
    - Fill the draft gaps ("???", "link", "MORE", "`????`").
12. **Remaining vignette text:**
    - Vignette 2: the humidity section text is copied from the pressure section.
    - Vignette 3: "**CHECK**" notes on humidity and wind.
    - Vignette 4: statistics plots are labelled "obs−bc_model" but show raw future minus corrected data.
    - pkgdown warns about missing figure alt-text in all 4 vignettes.
13. **Broken example paths:**
    - `extdata/sst/…2018_gridT.nc` in the `.sea_to_coast` and `create_ukcpsst_data` examples should be `extdata/sst_sample/` (2020/2021).
    - `extdata/preprepdata/ukcp18rcm.Rds` (`R/ukcp_functions.R`) doesn't exist.
    - The `.cad_conditions` example uses the deleted `mesoclim::climdata` and a non-existent `$dtmc`.
14. **`inst/extdata/data_scripts/pkg_data_prep.R:33`** writes `dtmc` over `dtmf.tif`. It should write to `dtmc.tif`. Check whether the shipped `dtmf.tif` is correct.

### P2 — packaging and tidying

15. **`man/` is out of sync with roxygen.** `devtools::document()` changes 13 `.Rd` files and deletes 3 (`blend_tile_lists`, `climdata`, `mosaicblend`), which are currently published on the site as stale pages. Run it, commit, and redeploy.
16. **DESCRIPTION:**
    - Remove `abind` (unused).
    - Drop `raster` (one `brick()` call in `R/workerfunctions.R`). Its full import masks `lubridate::union/intersect/origin` and `magrittr::extract`.
    - Add `grDevices`, `graphics` and `stats` to `Imports` (they're used via `importFrom`).
    - Decide on `curl`/`elevatr`, which are only needed for `get_dem()` (see Q1).
17. **Exported dot-functions:** `.resample`, `.sea_to_coast`, `.spatinterp`, `.tempcad`, `.tempcoastal`, `.tmeinterp`. Either remove the dot or remove `@export`. Vignette 4 uses `.spatinterp()`.
18. **Apparently unused code** (verify before deleting): `.writenc`, `savenc`, `plotrain`, `.tempcad`, `rainadjustv` (Rcpp); possibly `.cropnc`, `.clearskyraddaily`, `era5todaily`, `wet`.
19. **Unused data:**
    - `era5sst` (dataset).
    - `inst/extdata`: `era5/data_stream-*.nc` (0.8 MB), `dtms/era5lsm.tif` (1.6 MB), `dtms/era5dtmc.tif`, `dtms/shap_dtm1km.tif`, `dtms/altnaharra_1km.tif`, `haduk/altnaharra_wind_2020.tif` and `haduk/shap_rain_*`.
    - `dtms/era5dtm.tif` is now used only in examples.
20. **Missing `\examples`** in 27 `.Rd` files, including all bias-correction functions, `subset_climdata`, `landfill_climdata` and `calculate_windcoeffs`.

### Infrastructure

21. **Turn on GitHub Pages** (repo admin needed): Settings → Pages → Deploy from branch `gh-pages` / root. The site is deployed but returns 404 (`has_pages: false`).
22. **CI:** `usethis::use_github_action("check-standard")` and `use_github_action("pkgdown")`.
23. **Repository size:** `.git` is about 458 MB of history (old `climdata.rda`, `docs/`, `.o` files). Consider `git filter-repo` before any publication.
24. **Housekeeping:**
    - Delete the `dev` branch (local and remote) once it's no longer needed.
    - Review and drop the 3 stashes.
    - `analyses/issues_code.R` depends on deleted datasets.

### Open questions

- **Q1:** Keep `get_dem()` (and `curl`/`elevatr`), or delete it in favour of `microclimdata::dem_download()`?
- **Q2:** Keep or delete `.tempcad`, `era5todaily` and `landfill_climdata`? The last is exported and documented but has no example or caller.
- **Q3:** Should any `analyses/` scripts become articles or tests?
- **Q4:** Who owns the repo and site going forward (`ilyamaclean` vs `jrmosedale`)? This affects Pages activation and CI.

---

## Revision history

| Date | Summary |
|---|---|
| 2026-09-09 | Initial static audit of `dev` @ `eea5324`. |
| 2026-09-10 | Work on `claude/funny-goodall`, then fast-forwarded to `main`: basin delineation rewritten in C++ (O(N log N)); redundant helpers removed; `climdata` (39 MB) replaced by `bcmodel_list`; vignettes moved to `vignettes/`; dot-function exports cleaned; `.Rbuildignore`/`.gitignore` updated; `dev` merged (wind and precipitation fixes, `landfill_climdata`); pkgdown config rewritten and site deployed to `gh-pages`; README fixed; `old_ver` created. |
| 2026-09-24 | Local checkout moved from `dev` to `main`. `btw` MCP set up (Rust installed). `era5toclimarray()` fixed: missing `lsm`, cropping without `aoi`, CRS84/EPSG:4326 resampling shift (`6af9311`). All four vignettes made to run on package data (`4791e6a`). Site rebuilt and pushed to `gh-pages`. |
| 2026-09-25 | Audit restructured into setup / open actions / history; to-dos from `to_dos.md` merged in; audit moved to repository root. Tinker et al. (2024) SST reference and typo fixes added to vignettes. `claude/funny-goodall` worktree and branch deleted. |
