# coastalexposure requires a projected coordinate reference system

    Code
      coastalexposure(r, 270)
    Condition
      Error in `coastalexposure()`:
      ! landsea must have a projected coordinate reference system!!

# .check_cex crops larger exposure to dtmf and rejects mismatched grids

    Code
      .check_cex(terra::aggregate(cex, 2), sub)
    Condition
      Error in `.check_cex()`:
      ! cex must match the resolution and coordinate reference system of dtmf!!

