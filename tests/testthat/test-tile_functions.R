tile_template <- function(w, h, res = 50) {
  r <- terra::rast(terra::ext(0, w, 0, h), res = res, crs = "EPSG:27700")
  terra::values(r) <- 1
  r
}

tile_coverage <- function(r, tiles) {
  cov <- terra::rast(r)
  terra::values(cov) <- 0
  for (e in tiles$tile_extents) {
    idx <- terra::cells(r, e)
    cov[idx] <- cov[idx] + 1
  }
  terra::values(cov)
}

test_that("create_overlapping_tiles covers area with required overlap", {
  for (a in list(c(22000, 20000, 10000, 500), c(10000, 10000, 2000, 200),
                 c(47000, 33000, 5000, 1000), c(15500, 10000, 2000, 0))) {
    r <- tile_template(a[1], a[2])
    tiles <- suppressWarnings(create_overlapping_tiles(r, overlap = a[4], sz = a[3]))
    e <- t(sapply(tiles$tile_extents, as.vector))
    expect_true(all(tile_coverage(r, tiles) > 0))
    expect_true(all(e[, 1] >= 0 & e[, 2] <= a[1] & e[, 3] >= 0 & e[, 4] <= a[2]))
    xs <- sort(unique(e[, 1])); xe <- sort(unique(e[, 2]))
    if (length(xs) > 1) expect_gte(min(head(xe, -1) - tail(xs, -1)), a[4])
    expect_lte(max(e[, 2] - e[, 1]), 1.5 * a[3])
    expect_length(tiles$tile_land, length(tiles$tile_extents))
  }
})

test_that("create_overlapping_tiles handles tile size larger than area in one direction", {
  r <- tile_template(22000, 5000)
  tiles <- create_overlapping_tiles(r, overlap = 500, sz = 10000)
  expect_length(tiles$tile_extents, 2)
  expect_true(all(sapply(tiles$tile_extents, function(e) e$ymax - e$ymin) == 5000))
})

test_that("create_overlapping_tiles identifies tiles without land", {
  r <- tile_template(20000, 20000)
  terra::values(r) <- NA
  r[terra::cells(r, terra::ext(0, 4000, 0, 4000))] <- 1
  tiles <- create_overlapping_tiles(r, overlap = 200, sz = 5000)
  expect_equal(sum(tiles$tile_land == "y"), 1)
})

test_that("create_overlapping_tiles warns and errors on invalid inputs", {
  r <- tile_template(8000, 8000)
  expect_snapshot(tiles <- create_overlapping_tiles(r, overlap = 200, sz = 10000))
  expect_snapshot(create_overlapping_tiles(r, overlap = 5000, sz = 5000), error = TRUE)
  expect_snapshot(tiles <- create_overlapping_tiles(r, overlap = 200, sz = 2025))
})
