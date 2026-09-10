# mesoclim package audit — 9 September 2026

Static review of `~/OneDrive-UniversityofExeter/Rprojects/mesoclim` (branch `dev`, HEAD `eea5324` "Remove defunct pkg data and tests", in sync with `origin/dev`).

---

## Revision history

| Date | Branch | Changes |
|---|---|---|
| 2026-09-09 | `dev` @ `eea5324` | Initial audit |
| 2026-09-10 | `claude/funny-goodall` (from `main` @ `b264eea`) | See below |

### Changes made 10 September 2026 (branch `claude/funny-goodall`)

**Completed since audit — code changes (commits `49a69fc`, `680ff18`):**

- ✅ **Section 3 dead helpers** — `.edge()`, `.edgec()` deleted from `workerfunctions.R` (removed as part of basin merge rewrite)
- ✅ **Section 3 dead helper** — `.nctoarray_old()` deleted (was already broken and self-annotated)
- ✅ **Redundant wrapper removed** — `.clearskyrad()` deleted; all call sites in `data_functions.R` now call `clearskyrad()` directly with inline `jd`/`lt` extraction
- ✅ **Redundant wrapper removed** — `.skyview()` deleted; `calculate_terrain_shading()` gains `skyview_only=FALSE` flag, `steps` default changed 24→36, both call sites in `spatialdownscale.R` updated
- ✅ **Basin functions replaced** — `basinCpp` rewritten with min-heap priority queue (O(N log N), was O(N²)); old C++ helpers `whichmin`, `whichmin2`, `assignhigher`, `sel3n`, `sel3i`, `slotin` deleted; `.basinmerge()` + `.edge()` + `.edgec()` + `.asign3()` replaced by C++ `basinmerge_cpp()` (pour-point + union-find, single-pass transitive merge); public `basindelin()` gains `method="any"|"steepest"` parameter; all 63 tests pass
- ✅ **Sample NetCDF script added** — `data-raw/create_sample_ncdf.R` subsets ERA5, UKCP18RCM and SST files to <1 MB; outputs to `inst/extdata/era5raw/`, `inst/extdata/ukcprcm_sample/`, `inst/extdata/sst_sample/` (partial progress on P1 item 12 — sample data now exists, path references in examples still need updating)
- ✅ **`CLAUDE.md` added** — documents build/test/document/check/Rcpp commands and package architecture

**Completed since audit — P0 fixes (commit `9fe9103`):**

- ✅ **P0 item 1: `data/climdata.rda` removed** — 39 MB dataset deleted; replaced by `data/bcmodel_list.rda` (100 kB, tmin/tmax/prec bias-correction models for UKCP18 RCM member 01, Cornwall area); documented in `R/data.R`
- ✅ **P0 item 2 & 3: Vignettes fixed** — all four `.Rmd` files moved from `vignettes/articles/` to `vignettes/` with proper `%\VignetteIndexEntry{}` / `%\VignetteEngine{knitr::rmarkdown}` / `output: rmarkdown::html_vignette` metadata; global `eval=FALSE`; `library(mesoclim)` replaces `devtools::load_all()`; hard-coded OneDrive bc model path in vignette 2 replaced with `data(bcmodel_list, package="mesoclim")`
- ✅ **P0 item 5: Missing `.Rd` files fixed** — `@export` removed from 7 dot-functions in `ukcp_functions.R` that had `@noRd` (`.get_ukcp18_dates`, `.correct_ukcp_dates`, `.fill_calendar_data`, `.change_rast_units`, `.find_ukcp_decade`, `.lwup`, `.swdown`); all `#'` documentation blocks preserved inside the function files; `get_dem` fixed separately — `@noRd` removed, proper roxygen docs added with `@examples \dontrun{}` block; `devtools::document()` re-run; `get_dem.Rd` and `model_list.Rd` generated, NAMESPACE cleaned (7 dot-functions removed, `get_dem` retained)
- ✅ **P0 item 6: `.Rbuildignore` updated** — added `^docs$`, `^pkgdown$`, `^cds\.txt$`, `^README\.Rmd$`, `^\.github$`, `^analyses$`, `^R_new$`
- ✅ **P0 item 7: Compiled artefacts removed** — added `src/*.o`, `src/*.so`, `src/*.dll` to `.gitignore`; `git rm --cached` run on `src/RcppExports.o`, `src/mesoclimCpp.o`, `src/mesoclim.so`
- ✅ **`biascorrect_climdata()` partial model list** — changed hard `stop()` to `warning()` + `vars <- intersect(vars, model_vars)` so the function applies correction only to variables present in the supplied model list and warns about those missing (relevant because `bcmodel_list` covers only tmin/tmax/prec)

**Still present / not yet addressed:**
- ❌ **P1 items 8–13** — all still open except partial progress on item 12
- ❌ **P2 items 14–20** — all still open
- ❌ **Remaining dead helpers**: `.cropnc`, `.writenc`, `.clearskyraddaily` still in `workerfunctions.R`; `rainadjustv` Rcpp export still has no R caller
- ❌ **`dmean`, `.htd`** — not found on either `dev` or `claude/funny-goodall`; likely already deleted before the audit

**Branch note:** `claude/funny-goodall` was branched from `main` @ `b264eea`, not from `dev` @ `eea5324`. The changes above need merging into `dev` before the rest of the P0/P1 list is addressed there.

---

**Caveat:** R is not installed in the environment I can run code in, so nothing below comes from `R CMD check`, `devtools::check()` or executing any example. Everything is static analysis of source, `NAMESPACE`, `man/`, vignettes and git metadata. Anything marked *verify* needs an actual run.

---

## 1. What the package is

`mesoclim` 0.1.1 (GPL-2; Maclean `cre`, Gardner + Mosedale `aut`). Mechanistic downscaling of coarse climate data (UKCP18 RCM 12 km, ERA5 ~23 km) to sub-kilometre resolution, plus bias correction against observations.

**Scale:** 148 R functions in 11 files (~275 kB), 12 Rcpp routines in `src/mesoclimCpp.cpp`, 65 exports, 71 `.Rd`, 7 package datasets, ~4.6 MB of `inst/extdata` samples.

**Pipeline as coded:**

| Stage | Entry points |
|---|---|
| Preprocess | `ukcp18toclimarray`, `era5toclimarray`, `create_ukcpsst_data` → `checkinputs` → `write_climdata` / `read_climdata` |
| Bias correct | `biascorrect` / `precipcorrect` (fit) → `biascorrect_apply` / `precipcorrect_apply`, or `biascorrect_climdata` (whole climdata list) |
| Spatial downscale | `spatialdownscale` (wrapper), `spatialdownscale_tiles` (large areas), or per-variable: `tempdaily_downscale`, `temphrly_downscale`, `presdownscale`, `relhumdownscale`, `lwdownscale`, `swdownscale`, `winddownscale`, `precipdownscale` |
| Temporal downscale | `temporaldownscale` (wrapper) or `*_dailytohourly` per variable |
| Process/terrain helpers | `calculate_windcoeffs`, `basindelin`, `flowacc`, `coastalexposure`, `calculate_terrain_shading`, `lapserate`, `windelev`, `Tpsdownscale` |

**Data contract:** a `climdata` list — `dtm` (wrapped SpatRaster), `tme` (POSIXlt), `windheight_m`, `tempheight_m`, then `tmin/tmax/relhum/spechum/pres/swrad/lwrad/windspeed/winddir/prec` as SpatRaster or array. Consistent across preprocessing, checking, bias correction and downscaling — this is the strongest part of the design.

**Mesoclimate mechanisms implemented:** elevation lapse rate, cold-air drainage (basin delineation + flow accumulation, Rcpp-accelerated), coastal exposure with SST, wind-shelter coefficients by direction, terrain shading / skyview.

---

## 2. Prioritised actions

### P0 — blocks anyone else building or using the package

1. ✅ **`data/climdata.rda` removed** (2026-09-10) — 39 MB dataset deleted; replaced by `data/bcmodel_list.rda` (100 kB sample bias-correction models for tmin/tmax/prec). Note: `.cad_conditions` example still uses `mesoclim::climdata` — see P1 item 13.
2. ✅ **Vignettes fixed** (2026-09-10) — all four `.Rmd` files moved from `vignettes/articles/` to `vignettes/` with proper `%\VignetteIndexEntry{}` / `%\VignetteEngine{knitr::rmarkdown}` / `output: rmarkdown::html_vignette` metadata; `eval=FALSE` globally.
3. ✅ **Hard-coded OneDrive bc model path replaced** (2026-09-10) — vignette 2 now calls `data(bcmodel_list, package="mesoclim")`; `biascorrect_climdata()` updated to warn (not stop) for variables not covered by the model list.
4. ✅ **Broken chunk fence resolved** (2026-09-10) — the `{r era5check}` bare header without backticks was in the original `vignettes/articles/` version; static analysis of the current `vignettes/mesoclim_1_preparedata.Rmd` confirms all 7 chunk fences open and close correctly (final depth 0).
5. ✅ **Missing `.Rd` files fixed** (2026-09-10) — `@export` removed from 7 dot-functions in `ukcp_functions.R` (keeping all `#'` doc blocks); `get_dem` docs added with `\dontrun{}` example; `devtools::document()` re-run; NAMESPACE cleaned.
6. ✅ **`.Rbuildignore` updated** (2026-09-10) — added `^docs$`, `^pkgdown$`, `^cds\.txt$`, `^README\.Rmd$`, `^\.github$`, `^analyses$`, `^R_new$`.
7. ✅ **Compiled artefacts removed from git** (2026-09-10) — `src/*.o`, `src/*.so`, `src/*.dll` added to `.gitignore`; `git rm --cached` run on tracked `.o`/`.so` files.

### P1 — correctness

8. **Open issue #6 (unfixed in `dev`):** `.resample` (`R/workerfunctions.R:122–138`) calls `terra::resample()` without BIGTIFF, so `spatialdownscale` → `winddownscale` → `.resample` fails on large areas with *"TIFFAppendToStrip: Maximum TIFF file size exceeded"*. Reporter suggests `gdal = c("BIGTIFF=TRUE")`.
9. **Open issue #7 / PR #8 (unfixed in `dev`):** `subset_climdata` (`R/data_functions.R:370`) returns `newdata`, which is never bound when `climdata` is neither a `SpatRaster` nor a `list` (e.g. `NA`) → "object 'newdata' not found". Note the reporter's *second* claim — `else` on a new line at lines 379–381 — is not a bug: that code is inside a `{}` function body, where R parses it fine. Fix the first, close out the second explicitly.
10. **Two documented, reproducible science problems** already written up in `analyses/issues_code.R` but with no home in the package:
    - **A.** Downscaled wind speed at high elevation is far lower than HadUK monthly `sfcWind` (Altnaharra case). `analyses/new_windz.R` (28 kB, Nov 2025) looks like the candidate replacement for `.windz` / `windelev`.
    - **B.** Precipitation bias correction gives different answers one-step (`precipcorrect(..., mod_out=FALSE)`) vs two-step (`mod_out=TRUE` then `precipcorrect_apply`). It should not. `analyses/new_precipdownscale.R` is the candidate replacement.
    `issues_code.R` no longer runs — it needs `mesoclim::altnaharra_era5daily_2020`, `mesoclim::shap_era5daily_2020`, `extdata/dtms/altnaharra_5km.tif`, `extdata/haduk/altnaharra_rain_12km_2020.tif`, `extdata/haduk/altnaharra_rain_2020.tif`, all of which the `eea5324` commit removed.
11. **7 of 8 test files are 64-byte stubs.** Only `test-workerfunctions.R` (7.3 kB) has content. `test-spatialdownscale.R`, `test-temporaldownscale.R`, `test-biascorrect_functions.R`, `test-ukcp_functions.R`, `test-era5_functions.R`, `test-climatefunctions.R`, `test-data_functions.R` are empty. Issues #6 and #7 are exactly the kind of thing a regression test would have caught.
12. **Broken `system.file()` paths in shipped code and examples** (all resolve to nothing):
    - `extdata/sst/NWSClim_NWSPPE_r001i1p00000_2018_gridT.nc` — in `R/workerfunctions.R` (`.sea_to_coast` example) and `man/dot-sea_to_coast.Rd`. Actual directory is `extdata/sst_sample/`, years 2020 and 2021.
    - `extdata/sst` — in `R/ukcp_functions.R` and `man/create_ukcpsst_data.Rd`.
    - `extdata/preprepdata/ukcp18rcm.Rds` — in `R/ukcp_functions.R`.
13. **`.cad_conditions` example is doubly broken**: uses `mesoclim::climdata` (the 41 MB dataset, item 1) and passes `climdata$dtmc`, which is not a field of the climdata list (it is `dtm`).

### P2 — tidying

14. **`README.md`/`README.Rmd` reference five vignettes that do not exist**: `mesoclim_1_download`, `mesoclim_2_preparedata`, `mesoclim_3_spdscale`, `mesoclim_4_tmpdscale`, `mesoclim_5_biascorrect`. Actual files are `mesoclim_1_preparedata`, `_2_spdscale`, `_3_tmpdscale`, `_4_biascorrect`. This is open issue #1 (MarieFinocchiaro, Jan 2025) — the download vignette was dropped in favour of pointing at `microclimdata`, but the README never followed.
15. **`docs/` is committed and badly stale** — 272 tracked files, last built 30 April 2026, with 72 reference pages for functions that no longer exist (`SSTinterpolate`, `blendtemp_hadukera5`, `blend_tile_lists`, `create_parcel_list`, `download_ukcp18`, `hourlytemp`, `tempdownscale`, `write_parcels`, …). Also contains an `articles/articles/` directory using the *old* vignette names.
16. **`_pkgdown.yml` points at the wrong site**: `url: https://jrmosedale.github.io/mesoclim`, but the repo is `ilyamaclean/mesoclim` and `origin/gh-pages` exists there.
17. **`DESCRIPTION` dependency hygiene**:
    - `abind` — zero uses anywhere. Remove.
    - `curl`, `elevatr` — used only by `get_dem`, which is exported, undocumented and never called. Decide `get_dem`'s fate first (see Q4).
    - `raster` — `import(raster)`, one `brick(` call. Masks terra generics; drop it.
    - `gstat` (`gstat`/`idw`), `mgcv` (`gam`), `fields` (`Tps`), `fmsb` (`radarchart`), `zoo` (`na.approx`), `sf`, `units`, `ncdf4`, `lubridate` are all genuinely used.
    - `grDevices`, `graphics`, `stats` are in `NAMESPACE` `importFrom()` but missing from `DESCRIPTION` `Imports:` → `R CMD check` WARNING.
18. **Seven `.Rd` files document dot-functions or non-exported functions.** `dot-flowdir.Rd` lacks `\keyword{internal}` (the others have it). Separately, seven `.` functions are `@export`ed (`.resample`, `.spatinterp`, `.tmeinterp`, `.sea_to_coast`, `.tempcad`, `.tempcoastal`, `.lwup`, `.swdown`, `.change_rast_units`, `.correct_ukcp_dates`, `.fill_calendar_data`, `.find_ukcp_decade`, `.get_ukcp18_dates`) — either they are API and should lose the dot, or they are internal and should lose the `@export`.
19. **16 `.Rd` files have no `\examples{}`**, including all four core bias-correction functions: `biascorrect`, `biascorrect_apply`, `biascorrect_climdata`, `precipcorrect`, `precipcorrect_apply`, plus `calculate_windcoeffs`, `get_nearest_val`, `invls_calc`, `landfill_climdata`, `plot_timestats_r`.
20. **`inst/extdata/data_scripts/pkg_data_prep.R` has a copy-paste bug at line ~33**: `writeRaster(dtmc, "inst/extdata/dtms/dtmf.tif", overwrite=TRUE)` — writes the *coarse* DTM over `dtmf.tif`. Should be `dtmc.tif`. Worth checking whether the shipped `dtmf.tif` is correct.

---

## 3. Uncalled or non-working code

### Defined but never called anywhere (R/, man/, vignettes/, tests/, inst scripts)

| Function | File | Note |
|---|---|---|
| `get_dem` | workerfunctions.R:1620 | **exported**, no `.Rd`; sole reason `curl` + `elevatr` are dependencies |
| `landfill_climdata` | data_functions.R:431 | **exported**, has `.Rd`, no example |
| `.tempcad` | workerfunctions.R | **exported**, has `.Rd`, no example — superseded by `.apply_cad`/`.cad_conditions`? |
| `era5todaily` | biascorrect_functions.R:386 | not exported, has `.Rd` with `\keyword{internal}`; large function, looks orphaned |
| `savenc` | biascorrect_functions.R:624 | helper for `era5todaily` |
| `plotrain` | temporaldownscale.R:814 | plotting helper |
| `wet` | temporaldownscale.R | |
| `dmean`, `.htd` | workerfunctions.R | dead helpers — **not found on either branch; likely already deleted** |
| ~~`.edge`, `.edgec`~~ | workerfunctions.R | ✅ **deleted 2026-09-10** (removed with basin merge rewrite) |
| `.cropnc`, `.writenc` | workerfunctions.R | dead helpers — still present |
| `rainadjustv` | RcppExports.R | Rcpp routine with no R caller (`rainadjustm` is used) — still present |

`biascorrect_climdata` has no caller inside `R/` but *is* used in vignette 2 and `calc_biascorrect_models.R` — live, just untested.

### Broken

- `subset_climdata` — unbound `newdata` (issue #7). ❌ still open
- `.resample` — BIGTIFF failure on large areas (issue #6). ❌ still open
- `.cad_conditions` example — `mesoclim::climdata` + non-existent `$dtmc`. ❌ still open (blocked by P0 item 1)
- `.sea_to_coast`, `create_ukcpsst_data` examples — missing `extdata/sst/…` paths. ⚠️ sample SST files now exist at `extdata/sst_sample/` (added 2026-09-10) but example paths still reference old `extdata/sst/` directory
- `R/ukcp_functions.R` — missing `extdata/preprepdata/ukcp18rcm.Rds`. ❌ still open
- Vignette 1 chunk fence; vignette 2 hard-coded path. ❌ still open
- `analyses/issues_code.R` — depends on four deleted datasets and three deleted extdata files. ❌ still open
- Vignette 3 §"Precipitation downscaling" is marked **UNDER DEVELOPMENT** with the only call commented out. ❌ still open

No undefined internal (`.`-prefixed) function calls were found — the internal call graph is coherent.

---

## 4. Package data

### `data/` — ~2.1 MB total (was 41.1 MB)

| Dataset | Size | Explicit `mesoclim::` / `data()` refs | Verdict |
|---|---|---|---|
| ~~`climdata.rda`~~ | ~~41.0 MB~~ | — | ✅ **deleted 2026-09-10** |
| `bcmodel_list.rda` | 100 kB | vignette 2 (`data(bcmodel_list)`) | ✅ added 2026-09-10; covers tmin/tmax/prec only |
| `ukcp18sst.rda` | 54 kB | 11 | keep |
| `ukcpinput.rda` | 17 kB | 54 | keep — the workhorse |
| `era5sst.rda` | 24 kB | **0** | unused |
| `ukcpfuture.rda` | 18 kB | 2 (vignette 4) | keep |
| `landsea.rda` | 2.5 kB | **0** (2 examples use it as a variable name) | unused |
| `ukcp18lookup.rda` | 0.6 kB | **0** | unused, but documented and plausibly useful |

CRAN's limit is 5 MB for the whole tarball. Present `data/` alone is 8× that; `data/` + `inst/extdata` + `docs/` puts the tarball near 55 MB.

### `inst/extdata/` — 4.6 MB, mostly fine

- Referenced and present: `dtms/dtmf.tif`, `dtms/dtmm.tif`, `dtms/lizard50m.tif`, `dtms/altnaharra_1km.tif`, `dtms/era5dtm.tif`, `ukcp18rcm/orog_land-rcm_uk_12km_osgb.nc`, `ukcprcm_sample/*` (12 nc), `sst_sample/*` (2 nc, 1.5 MB), `era5raw/*`, `haduk/*`, `biascorrect/*`.
- Referenced but **missing**: `dtms/altnaharra_5km.tif`, `haduk/altnaharra_rain_12km_2020.tif`, `haduk/altnaharra_rain_2020.tif`, `sst/`, `preprepdata/ukcp18rcm.Rds`, `dtms/dtmf_inland.tif`.
- Present but **unreferenced**: `dtms/era5lsm.tif` (1.6 MB — the largest single unused file), `dtms/era5dtmc.tif`, `era5/data_stream-oper_stepType-{accum,avg,instant}.nc` (0.8 MB combined; the vignette uses `era5raw/` instead).

### Repository weight

`.git` pack is **458 MB** for 453 tracked files. That is not explained by the current tree (~55 MB) — it is history: earlier versions of `climdata.rda`, `docs/`, `src/*.o` and previously-committed sample NetCDFs. Fresh `git clone` and `install_github` are both slow as a result. If the repo is going public with a paper, a history rewrite (`git filter-repo`) or a fresh-start repo is worth considering.

`raw_data_files/` (2.0 GB) and `analyses/` are correctly excluded from both git and the build.

---

## 5. GitHub setup

| Item | State |
|---|---|
| Remote | `https://github.com/ilyamaclean/mesoclim.git` (HTTPS) |
| Local branch | `dev` @ `eea5324`, **identical to `origin/dev`** |
| `main` | local and `origin/main` both `b264eea` |
| `origin/gh-pages` | exists (`b193ea2`) |
| Last fetch | 30 April 2026 |
| GitHub Actions | **none** — no `.github/` directory locally or workflows in the repo |
| GitHub Pages | **not linked** in the repo's About; no website field set |
| Open issues | #1 (missing download vignette), #6 (large areas / BIGTIFF), #7 (subset_climdata) |
| Open PRs | #8 "Fix: syntax & handle unexpected input in subset_climdata" (timothee-bacri, 30 Oct 2025) — **unreviewed for 10 months** |
| `.git/config` | contains 43 duplicated `github-pr-owner-number` lines for PR #5 (dev→main) — harmless cruft from the VS Code GitHub extension |

### Local repository debris

- **Uncommitted deletions**: `vignettes/articles/mesoclim_1_preparedata.Rmd`, `_2_spdscale.Rmd`, `_3_tmpdscale.Rmd`, `_4_biascorrect.Rmd` and `vignettes/vignette_1_preparedata.Rmd` are in the index but gone from disk.
- **Zero-byte ref** `.git/refs/heads/claude` — a broken branch pointer.
- **Orphan worktree** `.claude/worktrees/funny-goodall/` with an old `man/` tree (including `blendtemp_hadukera5.Rd` for a deleted function), plus a registered worktree in `.git/worktrees/funny-goodall`. Left over from an earlier agent session; `git worktree prune` + delete.
- **A stash** (`refs/stash` → `740185c`) of unknown vintage.
- Git operations intermittently fail on this OneDrive mount ("Resource deadlock avoided", bus errors) because files are cloud-only placeholders. Anything doing real git work needs the repo fully hydrated, or moved off OneDrive.

### Recommended GitHub actions

1. Review and merge/close PR #8; fix #6 and #7 on `dev`; close #1 with the README fix.
2. Add `.github/workflows/R-CMD-check.yaml` and `pkgdown.yaml` (`usethis::use_github_action()`); once the pkgdown action publishes to `gh-pages`, delete `docs/` from the repo.
3. Fix `_pkgdown.yml` `url:` to match wherever the site actually lives, set the repo's website field, and fix `BugReports:` in `DESCRIPTION` — it is currently `http://github.com/mesoclim/issues`, which is not a valid URL. Should be `https://github.com/ilyamaclean/mesoclim/issues`.
4. Decide the `dev` → `main` merge cadence (PR #5 is the standing route).

---

## 6. Questions

1. **`climdata` (41 MB):** it and `ukcpinput` are documented identically. What is actually in it, and does anything still need it? If the answer is "nothing", deleting it fixes the single biggest packaging problem outright.
2. **Which functions do you want replaced?** My read of `analyses/` says `.windz`/`windelev`/`winddownscale` (from `new_windz.R`, and the high-elevation wind bias in `issues_code.R`) and `precipdownscale`/`precipcorrect` (from `new_precipdownscale.R`, and the 1-step vs 2-step inconsistency). Confirm, and say whether the new versions in `analyses/` are ready to promote or still exploratory.
3. **Vignettes: true vignettes or pkgdown articles?** True vignettes must build from package data alone on a clean machine — which rules out the OneDrive bias-correction models and pushes towards small shipped model objects. Articles can assume more, but then `VignetteBuilder: knitr` should go.
4. **`get_dem`, `landfill_climdata`, `.tempcad`, `era5todaily`:** keep and document, or delete? Vignette 1 explicitly says the package does not source climate data and points to `microclimdata`, which argues for deleting `get_dem` (and with it `curl` and `elevatr`).
5. **Target: CRAN, or GitHub-only?** CRAN forces the 5 MB limit, `R CMD check --as-cran` cleanliness, and no dot-function exports. GitHub-only is much more forgiving. This decision drives roughly half the P0/P2 list.
6. **`analyses/` (34 scripts, ~340 kB, gitignored):** any of these destined to become vignettes, articles, or tests — or are they permanently scratch?
7. **Who owns the repo going forward?** `_pkgdown.yml` says `jrmosedale`, `DESCRIPTION` and the remote say `ilyamaclean`. Worth settling before wiring up CI and Pages.
